"""
AI-Guided Molecular Docking Pipeline
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

    # Meeko >= 0.5 expects output basename (without .pdbqt) + flag '-p'
    out_basename = str(output_dir / pdb_file.stem)

    run_cmd = _build_cmd(
        cmd,
        [
            "-i", str(pdb_file),
            "-o", out_basename,
            "-p",
            "--default_altloc", "A",
            "--allow_bad_res",
        ],
    )
    log(f"  cmd : {' '.join(str(x) for x in run_cmd)}", "dim")

    result = _run(run_cmd, timeout=300)
    _relay(result, log)

    # ── Validate output ───────────────────────────────────────────────────
    if result.returncode != 0:
        hint = _receptor_hint(result.stderr)
        raise RuntimeError(
            f"mk_prepare_receptor failed (exit {result.returncode}).\n"
            f"stderr: {result.stderr[:600]}\n"
            f"{hint}"
        )

    if not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
        raise RuntimeError(
            "mk_prepare_receptor produced no output PDBQT file.\n"
            "Possible causes:\n"
            "  - The cleaned PDB has no ATOM records\n"
            "  - Meeko could not assign atom types to non-standard residues\n"
            "  - Output file path was blocked or invalid"
        )

    log(
        f"  -> {out_pdbqt.name}  ({out_pdbqt.stat().st_size // 1024} KB)",
        "dim",
    )
    return out_pdbqt


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


def _build_cmd(cmd_str: str, extra_args: List[str]) -> List[str]:
    """Resolve Python scripts or standalone system binaries to complete execution arrays."""
    cmd_str = cmd_str.strip()
    p = Path(cmd_str)

    if cmd_str.lower().endswith(".py"):
        venv_python = p.parent / ("python.exe" if os.name == "nt" else "python")
        if venv_python.exists():
            interpreter = str(venv_python)
        else:
            venv_python3 = p.parent / ("python3.exe" if os.name == "nt" else "python3")
            interpreter = str(venv_python3) if venv_python3.exists() else sys.executable
        return [interpreter, "-s", str(p)] + extra_args

    if p.is_absolute() or "/" in cmd_str or "\\" in cmd_str:
        return [str(p)] + extra_args

    resolved = shutil.which(cmd_str)
    return [resolved if resolved else cmd_str] + extra_args


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
    combined = (result.stdout or "") + (result.stderr or "")
    for raw in combined.splitlines():
        line = raw.strip()
        if not line:
            continue
        upper = line.upper()
        if "ERROR" in upper or "TRACEBACK" in upper:
            log(f"    Meeko > {line}", "error")
        elif "WARNING" in upper or "WARN" in upper:
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
        "Hint: Verify mk_prepare_receptor is from Meeko >= 0.5 and that\n"
        "the protein structure contains proper 3D ATOM entries."
    )