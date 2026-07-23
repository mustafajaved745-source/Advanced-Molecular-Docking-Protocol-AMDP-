"""
AI-Guided Molecular Docking Pipeline
=====================================
receptor_prep.py — Receptor cleaning (PyMOL) + docking grid generation

Responsibilities:
  1. Generate and run a PyMOL Python script that:
       - Loads the raw receptor PDB
       - Removes solvent (HOH/WAT) and non-significant heteroatoms
       - Isolates a SINGLE target chain containing the ligand (drops extra chains)
       - Retains essential cofactors (e.g., NAD/NADH) required for binding
       - Extracts the native ligand as both PDB and SDF from that single chain
       - Saves the cleaned receptor (protein + cofactor) as PDB
  2. Parse the native ligand PDB to compute the docking grid box
       - Centre  = centroid of ligand heavy atoms
       - Size    = bounding box + 2 × padding (default 8 Å per side)
       - Minimum box size enforced at 20 Å per axis

Public API
----------
  clean_receptor_with_pymol(receptor_pdb, ligand_code,
                             output_dir, pymol_exe, log_cb)
      -> (cleaned_pdb: Path, native_lig_pdb: Path, native_lig_sdf: Path)

  calculate_grid_box(native_ligand_pdb, padding=8.0)
      -> dict with keys center_x/y/z and size_x/y/z
"""

from __future__ import annotations

import json
import subprocess
import textwrap
from pathlib import Path
from typing import Callable, Dict, List, Tuple

# ════════════════════════════════════════════════════════════════════════════
#  Non-significant heteroatom residue names & cofactors
# ════════════════════════════════════════════════════════════════════════════

_CLEANUP_RESNS: set[str] = {
    # Water / deuterium water
    "HOH", "WAT", "DOD", "H2O",
    # Alkali / alkaline-earth / transition-metal ions
    "NA", "K", "LI", "RB", "CS", "MG", "CA", "SR", "BA",
    "ZN", "CU", "FE", "MN", "CO", "NI", "CD", "HG", "PB",
    "CR", "MO",
    # Halide / noble-gas ions
    "CL", "BR", "IOD", "F", "XE", "KR",
    # Inorganic polyatomic ions
    "SO4", "SUL", "PO4", "PHO", "NO3", "CO3", "BO4", "SCN", "NCO",
    # Cryoprotectants
    "GOL", "EDO", "EG", "PEG", "PG4", "PE4", "PE5", "P6G",
    "MPD", "MRD", "DMS", "DMSO",
    # Crystallisation buffers / additives
    "ACT", "ACY", "ACE", "FMT", "TRS", "MES", "HEP", "EPE",
    "PIP", "IMD", "DTT", "DTE", "BME", "MSE",
    "TAR", "TLA", "SUC", "TBU", "NH4", "NH2",
    "EOH", "MOH", "IPH", "PGE", "PG",
    "BU1", "BU2", "BTB", "MLA", "MLI", "MOL",
    "CIT", "ISO",
    # Detergents / lipids
    "OLA", "PLM", "SDS", "LMT", "C8E",
    # Amino-acid capping groups
    "NME",
}

_PRESERVED_COFACTORS: set[str] = {
    "NAD", "NDP", "NAP", "FAD", "FMN", "HEM", "PLP", "CoA"
}


# ════════════════════════════════════════════════════════════════════════════
#  Main Receptor Preparation API
# ════════════════════════════════════════════════════════════════════════════

def clean_receptor_with_pymol(
    receptor_pdb: str | Path,
    ligand_code: str,
    output_dir: Path,
    pymol_exe: str,
    log_cb: Callable[[str, str], None] | None = None,
) -> Tuple[Path, Path, Path]:
    """
    Run PyMOL in headless mode to clean the receptor and extract
    the native co-crystallised ligand from a single isolated chain.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    output_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = Path(receptor_pdb).resolve()
    cleaned_pdb = output_dir / "cleaned_receptor.pdb"
    native_lig_pdb = output_dir / "native_ligand.pdb"
    native_lig_sdf = output_dir / "native_ligand.sdf"
    script_path = output_dir / "_pymol_cleanup.py"

    # Build residue selection filters
    target_ligand = ligand_code.upper()
    excl_resns = [r for r in _CLEANUP_RESNS if r != target_ligand]
    removal_sel = " or ".join([f"resn {r}" for r in excl_resns])
    cofactor_sel = " or ".join([f"resn {c}" for c in _PRESERVED_COFACTORS])

    # Safely escape file paths into Python string representations
    script_content = textwrap.dedent(f"""\
        # Auto-generated PyMOL cleanup script
        import sys
        import shutil
        from pymol import cmd

        RECEPTOR_PDB   = {json.dumps(str(receptor_path))}
        CLEANED_PDB    = {json.dumps(str(cleaned_pdb))}
        NATIVE_LIG_PDB = {json.dumps(str(native_lig_pdb))}
        NATIVE_LIG_SDF = {json.dumps(str(native_lig_sdf))}
        LIGAND_CODE    = {json.dumps(target_ligand)}
        REMOVAL_SEL    = {json.dumps(removal_sel)}
        COFACTOR_SEL   = {json.dumps(cofactor_sel)}

        # 1. Load Structure & Remove Solvent
        cmd.load(RECEPTOR_PDB, "receptor")
        cmd.remove("solvent")
        
        # 2. Identify Target Chain
        target_chains = [c for c in cmd.get_chains(f"resn {{LIGAND_CODE}}") if c]
        if not target_chains:
            print(f"ERROR: Ligand '{{LIGAND_CODE}}' not found in structure.")
            cmd.quit()
            sys.exit(1)

        chosen_chain = "A" if "A" in target_chains else target_chains[0]
        print(f"Selected Chain '{{chosen_chain}}' containing ligand {{LIGAND_CODE}}")

        # 3. Isolate Target Chain & Clean Non-significant Heteroatoms
        cmd.remove(f"not chain {{chosen_chain}}")

        if REMOVAL_SEL:
            removed = cmd.count_atoms(f"hetatm and ({{REMOVAL_SEL}})")
            cmd.remove(f"hetatm and ({{REMOVAL_SEL}})")
            print(f"Removed {{removed}} non-significant heteroatom(s)")

        # 4. Save Native Ligand
        lig_sel = f"resn {{LIGAND_CODE}} and chain {{chosen_chain}}"
        lig_count = cmd.count_atoms(lig_sel)
        print(f"Target ligand {{LIGAND_CODE}} count: {{lig_count}} atom(s)")

        cmd.save(NATIVE_LIG_PDB, lig_sel)

        try:
            cmd.save(NATIVE_LIG_SDF, lig_sel)
        except Exception as err:
            print(f"WARNING: SDF export failed ({{err}}). Falling back to copying PDB.")
            shutil.copy(NATIVE_LIG_PDB, NATIVE_LIG_SDF)

        # 5. Prepare & Save Cleaned Receptor (Protein + Essential Cofactors)
        cmd.remove(f"resn {{LIGAND_CODE}}")

        catch_all_sel = f"hetatm and not ({{COFACTOR_SEL}})" if COFACTOR_SEL else "hetatm"
        rem_het = cmd.count_atoms(catch_all_sel)
        if rem_het > 0:
            print(f"Removing {{rem_het}} non-cofactor heteroatom(s)")
            cmd.remove(catch_all_sel)

        protein_atoms = cmd.count_atoms("all")
        if protein_atoms == 0:
            print("ERROR: No atoms remaining after cleanup.")
            cmd.quit()
            sys.exit(1)

        cmd.save(CLEANED_PDB, "all")
        print(f"Saved cleaned receptor: {{protein_atoms}} atoms -> {{CLEANED_PDB}}")

        cmd.quit()
    """)

    script_path.write_text(script_content, encoding="utf-8")
    log(f"PyMOL script generated: {script_path.name}", "dim")

    try:
        log(f"Running PyMOL engine ({Path(pymol_exe).name})...", "info")
        result = _run_pymol(pymol_exe, script_path, timeout=180)
        _relay_pymol_output(result.stdout + result.stderr, log)

        if result.returncode not in (0, 1):
            raise RuntimeError(
                f"PyMOL execution failed with exit code {result.returncode}.\n"
                f"Stderr:\n{result.stderr[:800]}"
            )

        _verify_outputs(
            [
                (cleaned_pdb, "cleaned receptor PDB"),
                (native_lig_pdb, "native ligand PDB"),
                (native_lig_sdf, "native ligand SDF"),
            ],
            target_ligand,
        )
    finally:
        if script_path.exists():
            script_path.unlink(missing_ok=True)

    return cleaned_pdb, native_lig_pdb, native_lig_sdf


# ════════════════════════════════════════════════════════════════════════════
#  Grid Box Calculation API
# ════════════════════════════════════════════════════════════════════════════

def calculate_grid_box(
    native_ligand_pdb: Path,
    padding: float = 8.0,
    min_size: float = 20.0,
) -> Dict[str, float]:
    """Compute the AutoDock Vina docking grid box from a native ligand PDB."""
    coords = _parse_heavy_atom_coords(native_ligand_pdb)

    if not coords:
        raise ValueError(
            f"No valid heavy-atom coordinates found in {native_ligand_pdb}.\n"
            "Cannot calculate docking grid box."
        )

    xs, ys, zs = zip(*coords)

    # Compute bounding box centroid
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    cx = (max_x + min_x) / 2.0
    cy = (max_y + min_y) / 2.0
    cz = (max_z + min_z) / 2.0

    # Dimension = span + padding on both sides, capped at minimum size
    sx = max(max_x - min_x + 2.0 * padding, min_size)
    sy = max(max_y - min_y + 2.0 * padding, min_size)
    sz = max(max_z - min_z + 2.0 * padding, min_size)

    return {
        "center_x": round(cx, 4),
        "center_y": round(cy, 4),
        "center_z": round(cz, 4),
        "size_x": round(sx, 2),
        "size_y": round(sy, 2),
        "size_z": round(sz, 2),
    }


# ════════════════════════════════════════════════════════════════════════════
#  Internal Helper Functions
# ════════════════════════════════════════════════════════════════════════════

def _run_pymol(
    pymol_exe: str,
    script_path: Path,
    timeout: int = 180,
) -> subprocess.CompletedProcess[str]:
    """Execute PyMOL headlessly using compatible CLI flags."""
    for flag_group in (["-cq"], ["-c"], []):
        cmd = [pymol_exe, *flag_group, str(script_path)]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            if "unknown option" not in (result.stderr or "").lower():
                return result
        except FileNotFoundError:
            raise RuntimeError(
                f"PyMOL executable not found: {pymol_exe}\n"
                "Please verify the installation path in settings."
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"PyMOL execution timed out ({timeout} s).\n"
                "The receptor file may be abnormally large, or PyMOL hung."
            )

    # Final fallback attempt
    return subprocess.run(
        [pymol_exe, str(script_path)],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _relay_pymol_output(text: str, log_cb: Callable[[str, str], None]) -> None:
    """Parse PyMOL stdout/stderr and output structured logs."""
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        upper = line.upper()
        if "ERROR" in upper:
            log_cb(f"  PyMOL > {line}", "error")
        elif "WARNING" in upper:
            log_cb(f"  PyMOL > {line}", "warning")
        elif any(k in upper for k in ("SAVED", "ATOMS", "LOADED")):
            log_cb(f"  PyMOL > {line}", "info")
        else:
            log_cb(f"  PyMOL > {line}", "dim")


def _verify_outputs(
    files: List[Tuple[Path, str]],
    ligand_code: str,
) -> None:
    """Validate that required files exist and are not empty."""
    for path, label in files:
        if not path.exists():
            raise RuntimeError(
                f"PyMOL failed to create {label}: {path}\n"
                f"Ensure ligand code '{ligand_code}' exists in the receptor PDB."
            )
        if path.stat().st_size == 0:
            raise RuntimeError(
                f"PyMOL produced an empty file for {label}: {path}"
            )


def _parse_heavy_atom_coords(pdb_file: Path) -> List[Tuple[float, float, float]]:
    """Extract Cartesian coordinates for non-hydrogen atoms from a PDB file."""
    coords: List[Tuple[float, float, float]] = []

    with open(pdb_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                # Element symbol parsing based on standard PDB formatting (cols 77-78)
                element = line[76:78].strip().upper() if len(line) >= 78 else ""

                # Fallback to atom name column (cols 13-16) if element is absent
                if not element:
                    atom_name = line[12:16].strip()
                    element = "".join(c for c in atom_name if c.isalpha())[:1].upper()

                if element == "H":
                    continue

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append((x, y, z))

            except (ValueError, IndexError):
                continue

    return coords