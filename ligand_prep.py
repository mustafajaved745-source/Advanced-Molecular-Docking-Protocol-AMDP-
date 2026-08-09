"""
Molecular Docking Pipeline
=====================================
ligand_prep.py — Ligand + receptor PDBQT preparation via Meeko

Responsibilities:
  1. Convert ligand SDF files  → PDBQT  using mk_prepare_ligand
  2. Convert receptor PDB file → PDBQT  using mk_prepare_receptor
  3. Handle .py scripts, bare executables, and PATH-resident commands
  4. Surface clear, actionable error messages when Meeko fails

Public API
----------
  prepare_ligand_pdbqt(sdf_file, output_dir, cmd, prefix, log_cb)
      -> Path   (output .pdbqt)

  prepare_receptor_pdbqt(pdb_file, output_dir, cmd, log_cb)
      -> Path   (output .pdbqt)
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, List, Optional


# ════════════════════════════════════════════════════════════════════════════
#  Ligand preparation
# ════════════════════════════════════════════════════════════════════════════

def prepare_ligand_pdbqt(
    sdf_file: Path,
    output_dir: Path,
    cmd: str,
    prefix: Optional[str] = None,
    log_cb: Callable[[str, str], None] | None = None,
) -> Path:
    """
    Convert a single SDF (or PDB fallback) ligand file to PDBQT using
    Meeko's mk_prepare_ligand tool.

    Parameters
    ----------
    sdf_file   : Input SDF file path.
    output_dir : Directory where the .pdbqt will be written.
    cmd        : Path to the mk_prepare_ligand executable or .py script.
    prefix     : Stem name for the output file. Defaults to sdf_file.stem.
    log_cb     : Optional logging callback log_cb(message, level).

    Returns
    -------
    Path to the produced .pdbqt file.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    sdf_file = Path(sdf_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = prefix or sdf_file.stem
    out_pdbqt = output_dir / f"{stem}.pdbqt"

    # Pre-process ligand: add explicit hydrogens & isolate primary fragment
    input_file_for_meeko = _ensure_explicit_hs(sdf_file, output_dir, cmd, log)

    # Build and execute preparation command
    run_cmd = _build_cmd(cmd, ["-i", str(input_file_for_meeko), "-o", str(out_pdbqt)])
    log(f"  cmd : {' '.join(str(x) for x in run_cmd)}", "dim")

    result = _run(run_cmd, timeout=120)
    _relay(result, log)

    # ── Validate output ───────────────────────────────────────────────────
    if result.returncode != 0:
        hint = _meeko_hint(result.stderr, input_file_for_meeko)
        raise RuntimeError(
            f"mk_prepare_ligand failed for '{input_file_for_meeko.name}' "
            f"(exit {result.returncode}).\n"
            f"stderr: {result.stderr[:600]}\n"
            f"{hint}"
        )

    if not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
        raise RuntimeError(
            f"mk_prepare_ligand produced no output for '{input_file_for_meeko.name}'.\n"
            "Possible causes:\n"
            "  - The input file contains no valid 3D molecules\n"
            "  - Meeko could not perceive bond orders\n"
            "  - The output directory path was rejected"
        )

    log(f"  -> {out_pdbqt.name}  ({out_pdbqt.stat().st_size} bytes)", "dim")
    return out_pdbqt


# ════════════════════════════════════════════════════════════════════════════
#  Receptor preparation
# ════════════════════════════════════════════════════════════════════════════

def prepare_receptor_pdbqt(
    pdb_file: Path,
    output_dir: Path,
    cmd: str,
    target_chain: str = "",
    grid: dict[str, float] | None = None,
    log_cb: Callable[[str, str], None] | None = None,
) -> Path:
    """
    Convert a cleaned receptor PDB to PDBQT using Meeko's
    mk_prepare_receptor tool.

    Parameters
    ----------
    pdb_file   : Cleaned receptor PDB (output of receptor_prep.py).
    output_dir : Directory where the .pdbqt will be written.
    cmd        : Path to the mk_prepare_receptor executable or .py script.
    log_cb     : Optional logging callback.

    Returns
    -------
    Path to the produced .pdbqt file.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    pdb_file = Path(pdb_file).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    out_pdbqt = output_dir / f"{pdb_file.stem}.pdbqt"
    def _try_meeko(
        input_pdb: Path,
        timeout_s: int,
        args_05: list[str],
        args_04: list[str],
    ) -> subprocess.CompletedProcess | None:
        """Try Meeko 0.5+ then 0.4.x fallback. Returns result or None."""
        def _run_one(args: list[str], label: str) -> subprocess.CompletedProcess | None:
            rc = _build_cmd(cmd, args)
            log(f"  {label} : {' '.join(str(x) for x in rc)}", "dim")
            try:
                r = _run(rc, timeout=timeout_s)
            except RuntimeError as e:
                if "timed out" in str(e):
                    log(f"  {label} timed out after {timeout_s}s", "dim")
                    return None
                raise
            _relay(r, log)
            return r

        # ── Meeko >= 0.5 (basename + -p flag) ──────────────────────────
        result = _run_one(args_05, "cmd")
        if result and result.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            return result

        # ── Meeko 0.4.x (full .pdbqt path, no -p) ─────────────────────
        result = _run_one(args_04, "retry (0.4.x)")
        return result

    # ── Attempt 1: try with HETATMs (fast path, 60 s timeout) ────────────
    out_basename = str(output_dir / pdb_file.stem)
    result = _try_meeko(
        pdb_file, 60,
        ["-i", str(pdb_file), "-o", out_basename, "-p",
         "--default_altloc", "A", "--allow_bad_res"],
        ["-i", str(pdb_file), "-o", str(out_pdbqt),
         "--default_altloc", "A", "--allow_bad_res"],
    )

    # ── Validate output ───────────────────────────────────────────────────
    if result is None:
        raise RuntimeError(
            "TIMEOUT: mk_prepare_receptor exceeded the preparation timeout; "
            "the receptor model was not changed or retried without cofactors."
        )
    if result and result.returncode != 0:
        hint = _receptor_hint(result.stderr)
        raise RuntimeError(
            f"mk_prepare_receptor failed (exit {result.returncode}).\n"
            f"stderr: {result.stderr[:600]}\n"
            f"{hint}"
        )

    if not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
        raise RuntimeError(
            "INVALID_RECEPTOR_PDBQT: mk_prepare_receptor produced no output file.\n"
            "Possible causes:\n"
            "  - The cleaned PDB has no ATOM records\n"
            "  - Meeko could not assign atom types to non-standard residues\n"
            "  - Output file path was blocked or invalid"
        )

    integrity = validate_receptor_pdbqt_integrity(
        pdb_file, out_pdbqt, target_chain=target_chain, grid=grid
    )
    log(
        "  Receptor integrity: "
        f"{integrity['output_heavy_atoms']}/{integrity['input_heavy_atoms']} heavy atoms, "
        f"{integrity['output_residues']}/{integrity['input_residues']} residues retained",
        "dim",
    )

    log(
        f"  -> {out_pdbqt.name}  ({out_pdbqt.stat().st_size // 1024} KB)",
        "dim",
    )
    return out_pdbqt


def _structure_inventory(path: Path) -> dict:
    atoms: list[tuple[tuple[str, str, str, str], tuple[float, float, float]]] = []
    chains: set[str] = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            atom_name = line[12:16].strip()
            ad_type = line[77:79].strip().upper() if len(line) >= 79 else ""
            element = line[76:78].strip().upper() if path.suffix.lower() == ".pdb" else ad_type
            if not element:
                element = "".join(c for c in atom_name if c.isalpha())[:1].upper()
            if element in ("H", "HD", "HS", "D"):
                continue
            chain = line[21].strip()
            residue = (chain, line[17:20].strip(), line[22:26].strip(), line[26].strip())
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms.append((residue, xyz))
            chains.add(chain)
    return {"atoms": atoms, "residues": {a[0] for a in atoms}, "chains": chains}


def validate_receptor_pdbqt_integrity(
    cleaned_pdb: Path,
    receptor_pdbqt: Path,
    target_chain: str = "",
    grid: dict[str, float] | None = None,
    min_atom_fraction: float = 0.80,
    min_residue_fraction: float = 0.90,
    pocket_shell: float = 2.0,
) -> dict:
    """Fail closed if Meeko silently drops the receptor or target pocket."""
    source = _structure_inventory(Path(cleaned_pdb))
    output = _structure_inventory(Path(receptor_pdbqt))
    input_atoms = len(source["atoms"])
    output_atoms = len(output["atoms"])
    input_residues = len(source["residues"])
    output_residues = len(output["residues"])
    if input_atoms == 0 or input_residues == 0:
        raise RuntimeError("INVALID_RECEPTOR_INPUT: cleaned receptor has no heavy atoms/residues.")
    if output_atoms < input_atoms * min_atom_fraction or output_residues < input_residues * min_residue_fraction:
        raise RuntimeError(
            "INVALID_RECEPTOR_PDBQT: Meeko retained only "
            f"{output_atoms}/{input_atoms} heavy atoms and "
            f"{output_residues}/{input_residues} residues."
        )
    if target_chain and target_chain not in output["chains"]:
        raise RuntimeError(
            f"INVALID_RECEPTOR_PDBQT: selected target chain '{target_chain}' is absent."
        )
    if grid:
        bounds = [
            (
                grid[f"center_{axis}"] - grid[f"size_{axis}"] / 2.0 - pocket_shell,
                grid[f"center_{axis}"] + grid[f"size_{axis}"] / 2.0 + pocket_shell,
            )
            for axis in "xyz"
        ]
        pocket_residues = {
            residue for residue, xyz in source["atoms"]
            if all(bounds[i][0] <= xyz[i] <= bounds[i][1] for i in range(3))
        }
        missing = sorted(pocket_residues - output["residues"])
        if missing:
            preview = ", ".join(f"{r[0]}:{r[1]}{r[2]}{r[3]}" for r in missing[:8])
            raise RuntimeError(
                "INVALID_RECEPTOR_PDBQT: receptor residues inside/near the docking "
                f"box were omitted: {preview}"
            )
    return {
        "input_heavy_atoms": input_atoms,
        "output_heavy_atoms": output_atoms,
        "input_residues": input_residues,
        "output_residues": output_residues,
    }


# ════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ════════════════════════════════════════════════════════════════════════════

def _ensure_explicit_hs(
    sdf_file: Path,
    output_dir: Path,
    cmd: str,
    log: Callable[[str, str], None],
) -> Path:
    """
    Standardize ligand by adding explicit hydrogens and keeping the primary
    largest fragment via RDKit inside the target environment.
    """
    hs_sdf = output_dir / f"{sdf_file.stem}_hs.sdf"
    interpreter = _build_cmd(cmd, [])[0]

    # Clean inline script for execution
    script = (
        "import sys\n"
        "from rdkit import Chem\n"
        f"sdf_path = r'{sdf_file}'\n"
        f"out_path = r'{hs_sdf}'\n"
        "mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)\n"
        "if mol is not None:\n"
        "    try:\n"
        "        Chem.SanitizeMol(mol)\n"
        "    except Exception:\n"
        "        pass\n"
        "    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)\n"
        "    if len(frags) > 1:\n"
        "        print(f'WARNING: SDF contains {len(frags)} fragments; keeping the largest.')\n"
        "        mol = max(frags, key=lambda f: f.GetNumAtoms())\n"
        "    try:\n"
        "        mol = Chem.AddHs(mol, addCoords=True)\n"
        "    except Exception:\n"
        "        mol = Chem.AddHs(mol)\n"
        "    writer = Chem.SDWriter(out_path)\n"
        "    writer.write(mol)\n"
        "    writer.close()\n"
        "    print('SUCCESS')\n"
        "else:\n"
        "    print('FAILED_LOAD')\n"
    )

    try:
        result = subprocess.run(
            [interpreter, "-c", script],
            capture_output=True,
            text=True,
            timeout=40,
        )

        if result.returncode == 0 and "SUCCESS" in result.stdout and hs_sdf.exists():
            log(f"  Added explicit Hs → {hs_sdf.name}", "dim")
            return hs_sdf
    except Exception as e:
        log(f"  RDKit protonation check skipped: {e}", "dim")

    # Fallback: OpenBabel if present on the system
    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        log("  RDKit load failed; attempting OpenBabel fallback (-h)...", "warning")
        ob_cmd = [obabel_bin, str(sdf_file), "-O", str(hs_sdf), "-h"]
        ob_res = subprocess.run(ob_cmd, capture_output=True, text=True, timeout=60)
        if ob_res.returncode == 0 and hs_sdf.exists():
            log(f"  Added explicit Hs via OpenBabel → {hs_sdf.name}", "dim")
            return hs_sdf

    log("  Warning: Proceeding with original SDF (hydrogens could not be modified).", "warning")
    return sdf_file


# ════════════════════════════════════════════════════════════════════════════
#  2D → 3D Embedding
# ════════════════════════════════════════════════════════════════════════════

def ensure_3d_sdf(
    input_sdf: Path,
    output_sdf: Path,
    log: Callable[[str, str], None] | None = None,
) -> bool:
    """
    Check whether an SDF contains 2D (flat) coordinates and, if so,
    embed proper 3D coordinates using RDKit ETKDGv3 + MMFF94 optimization.

    Returns True if *output_sdf* was written successfully, False on failure.
    If the molecule already has valid 3D coordinates, copies it unchanged.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    def _log(msg: str, lvl: str = "info") -> None:
        if log:
            log(msg, lvl)

    input_sdf = Path(input_sdf)
    output_sdf = Path(output_sdf)

    if not input_sdf.exists():
        _log(f"  ensure_3d_sdf: input not found — {input_sdf}", "error")
        return False

    try:
        suppl = Chem.SDMolSupplier(str(input_sdf), removeHs=False, sanitize=False)
        mol = next(suppl, None)
    except Exception as e:
        _log(f"  ensure_3d_sdf: failed to read SDF — {e}", "error")
        return False

    if mol is None:
        _log("  ensure_3d_sdf: no valid molecule in SDF.", "error")
        return False

    # Sanitize the molecule (best-effort)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    mol = Chem.AddHs(mol, addCoords=True)

    # ── Detect 2D / flat coordinates ───────────────────────────────────
    is_3d = False
    try:
        conf = mol.GetConformer()
        # A conformer is considered 3D only if Is3D() is set AND at least
        # one atom has a non-negligible Z offset
        if conf.Is3D():
            for i in range(mol.GetNumAtoms()):
                if abs(conf.GetAtomPosition(i).z) > 0.5:
                    is_3d = True
                    break
    except Exception:
        pass

    if is_3d:
        _log(f"  SDF already has 3D coordinates — copying to {output_sdf.name}", "dim")
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(mol)
        writer.close()
        return True

    # ── Embed 3D coordinates ───────────────────────────────────────────
    _log("  SDF has 2D / flat coordinates — embedding 3D via ETKDGv3…", "info")
    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    embed_ok = AllChem.EmbedMolecule(mol, params)
    if embed_ok != 0:
        # Retry with ETKDGv2
        embed_ok = AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())
    if embed_ok != 0:
        # Retry with random coordinates
        try:
            embed_ok = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)
        except Exception:
            embed_ok = 1
    if embed_ok != 0:
        # Retry with distance geometry (fallback for rigid macrocycles)
        try:
            from rdkit.Chem import rdDistGeom
            params_dg = rdDistGeom.EmbedParameters()
            params_dg.randomSeed = 42
            embed_ok = rdDistGeom.EmbedMolecule(mol, params_dg)
        except (ImportError, Exception):
            embed_ok = 1
    if embed_ok != 0:
        _log("  ensure_3d_sdf: RDKit 3D embedding failed (all methods).", "error")
        return False

    # ── Geometry optimization (MMFF94) ─────────────────────────────────
    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94")
    except Exception as e:
        _log(f"  MMFF optimization skipped: {e}", "warning")

    # Write the 3D-embedded molecule
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()

    _log(f"  3D embedding complete → {output_sdf.name}", "success")
    return True


def _build_cmd(cmd_str: str, extra_args: List[str]) -> List[str]:
    """Resolve Python scripts or standalone system binaries to execution arrays.

    Cross-platform logic:
      - ``.py`` scripts → resolved via the sibling Python interpreter or
        ``sys.executable``; ``py`` (Windows launcher) as last fallback
      - Absolute / explicit paths → used directly (handles ``C:\\…`` on Windows)
      - Bare command names → resolved via ``shutil.which``
    """
    cmd_str = cmd_str.strip()
    p = Path(cmd_str)

    if cmd_str.lower().endswith(".py"):
        venv_python = p.parent / ("python.exe" if os.name == "nt" else "python")
        if venv_python.exists():
            interpreter = str(venv_python)
        else:
            venv_python3 = p.parent / ("python3.exe" if os.name == "nt" else "python3")
            if venv_python3.exists():
                interpreter = str(venv_python3)
            elif os.name == "nt" and shutil.which("py"):
                interpreter = "py"
            else:
                interpreter = sys.executable
        return [interpreter, "-s", str(p)] + extra_args

    if p.is_absolute():
        # Absolute path — use directly (works for both POSIX and Windows)
        return [str(p)] + extra_args

    # Relative path containing a separator — resolve by joining with CWD
    if "/" in cmd_str or "\\" in cmd_str:
        abs_path = (Path.cwd() / cmd_str).resolve()
        if abs_path.is_file():
            return [str(abs_path)] + extra_args
        return [cmd_str] + extra_args

    # Bare command — search PATH (shutil.which handles .exe appending on Windows)
    resolved = shutil.which(cmd_str)
    if resolved:
        return [resolved] + extra_args

    # Last attempt on Windows: try with .exe appended
    if os.name == "nt" and not cmd_str.lower().endswith(".exe"):
        resolved = shutil.which(cmd_str + ".exe")
        if resolved:
            return [resolved] + extra_args

    return [cmd_str] + extra_args


def _run(
    cmd: List[str],
    timeout: int = 120,
) -> subprocess.CompletedProcess:
    """Safely run a subprocess, capturing output."""
    try:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"Executable or interpreter not found: {cmd[0]}\n"
            "Please check tool paths in Settings."
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"Command timed out after {timeout} seconds:\n  {' '.join(str(c) for c in cmd)}"
        )


def _relay(
    result: subprocess.CompletedProcess,
    log: Callable[[str, str], None],
) -> None:
    """Relay stdout and stderr from process executions back to the log pump."""
    # Warnings that are harmless and handled automatically by Meeko
    # (e.g. unknown cofactor residues like NAD — Meeko builds templates on the fly)
    _SUPPRESS_WARNING_SUBSTRINGS = (
        "NOT IN RESIDUE_TEMPLATES",
        "TRYING TO RESOLVE UNKNOWN RESIDUES",
    )
    _ZERO_ERROR_SUMMARIES = {
        "PDBQT FILES NOT WRITTEN DUE TO ERROR: 0",
        "INPUT MOLECULES WITH ERRORS: 0",
    }

    combined = (result.stdout or "") + (result.stderr or "")
    for raw in combined.splitlines():
        line = raw.strip()
        if not line:
            continue
        upper = line.upper()
        if upper in _ZERO_ERROR_SUMMARIES:
            # Successful Meeko summary; keep it in the full log only.
            log(f"    Meeko > {line}", "dim")
        elif "ERROR" in upper or "TRACEBACK" in upper:
            log(f"    Meeko > {line}", "error")
        elif "WARNING" in upper or "WARN" in upper:
            # Downgrade known-harmless Meeko warnings to dim
            if any(sub in upper for sub in _SUPPRESS_WARNING_SUBSTRINGS):
                log(f"    Meeko > {line}", "dim")
            else:
                log(f"    Meeko > {line}", "warning")
        else:
            log(f"    Meeko > {line}", "dim")


def _meeko_hint(stderr: str, sdf_file: Path) -> str:
    """Return an actionable troubleshooting hint based on stderr messages."""
    s = (stderr or "").lower()

    if "sanitize" in s or "valence" in s or "kekulize" in s:
        return (
            "Hint: RDKit/Meeko could not sanitize the ligand.\n"
            "  - Inspect structure for invalid valences or bond types\n"
            "  - Pre-process with OpenBabel: obabel input.sdf -O fixed.sdf -h"
        )

    if "no such file" in s or "not found" in s:
        return f"Hint: File not found or unreadable. Path: {sdf_file}"

    if "importerror" in s or "modulenotfounderror" in s:
        return (
            "Hint: Missing required modules in Meeko environment.\n"
            "  Run: pip install meeko rdkit-pypi"
        )

    return (
        "Hint: Verify Meeko version (>= 0.4) and ensure valid 3D coordinates\n"
        "and explicit bond configurations exist in the SDF."
    )


def _receptor_hint(stderr: str) -> str:
    """Return an actionable hint for receptor conversion errors."""
    s = (stderr or "").lower()

    if "importerror" in s or "modulenotfounderror" in s:
        return "Hint: Ensure Meeko is installed (pip install meeko)."

    if "non-standard" in s or "unknown residue" in s:
        return (
            "Hint: Unrecognized non-standard amino acid/heteroatom found.\n"
            "  Clean or remove non-standard residues prior to preparation."
        )

    return (
        "Hint: Verify mk_prepare_receptor is from Meeko (0.4+ or 0.5+) and that\n"
        "the protein structure contains proper 3D ATOM entries."
    )
