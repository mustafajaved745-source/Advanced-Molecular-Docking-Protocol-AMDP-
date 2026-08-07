"""
Molecular Docking Pipeline
=====================================
receptor_prep.py — Receptor cleaning (PyMOL) + docking grid generation

Responsibilities:
  1. Generate and run a PyMOL Python script that:
       - Loads the raw receptor structure (PDB or mmCIF; PyMOL auto-detects)
       - Selects the user-chosen target chain
       - Removes all other chains
       - Removes all heteroatoms (ligands, waters, cryoprotectants, lipids, …)
       - Keeps only the selected protein chain and saves the cleaned receptor
       - The cleaned receptor is always written as PDB, so downstream
         Meeko / Vina steps never need to handle mmCIF themselves
  2. Extract the selected native ligand's crystal coordinates directly
     from the original receptor structure via pure Python (no PyMOL,
     no intermediate PDB file) for grid-box calculation and RMSD reference
  3. Compute the docking grid box from those coordinates
       - Centre  = centroid of ligand heavy atoms
       - Size    = bounding box + 2 × padding (default 4 Å per side)
       - Minimum box size enforced at 12 Å per axis

Public API
----------
  clean_receptor_with_pymol(receptor_path, ligand_code, target_chain,
                             output_dir, pymol_exe, log_cb)
      -> cleaned_pdb: Path

  extract_native_ligand_coords(receptor_path, ligand_code, target_chain)
      -> List[Tuple[float, float, float]]  (heavy-atom crystal coordinates)

  calculate_grid_box(coords, padding=4.0)
      -> dict with keys center_x/y/z and size_x/y/z
"""

from __future__ import annotations

import json
import logging
import subprocess
import textwrap
from pathlib import Path
from typing import Callable, Dict, List, Tuple

from hetatm_parser import iter_hetatm_atoms


# ════════════════════════════════════════════════════════════════════════════
#  Main Receptor Preparation API
# ════════════════════════════════════════════════════════════════════════════

def clean_receptor_with_pymol(
    receptor_path: str | Path,
    ligand_code: str,
    target_chain: str,
    output_dir: Path,
    pymol_exe: str,
    preserved_cofactor: str | None = None,
    log_cb: Callable[[str, str], None] | None = None,
) -> Path:
    """
    Run PyMOL in headless mode to clean the receptor:

      1. Select the user-chosen target chain
      2. Remove all other chains except the selected chain
      3. Deselect the current selection
      4. Remove all heteroatoms (ligands, water, cryoprotectants, lipids, …)
         except any preserved cofactor
      5. Keep only the selected protein chain and save the cleaned receptor

    Native-ligand extraction is intentionally removed because PyMOL-based
    PDB/SDF export often produces broken fragment files.  Grid-box origin
    is now derived externally from the input ligand file instead.

    Parameters
    ----------
    target_chain : str
        Single-character chain ID to isolate (e.g. "A", "B").
        Chosen explicitly by the user in the Active Site Selection panel —
        there is no automatic chain detection.
    ligand_code : str
        Accepted for API compatibility; no longer used in PyMOL script.
    preserved_cofactor : str or None
        3-letter residue name of an essential cofactor (e.g. "NAD") to
        preserve during heteroatom removal, or None to remove everything.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    output_dir.mkdir(parents=True, exist_ok=True)
    receptor_path = Path(receptor_path).resolve()
    cleaned_pdb = output_dir / "cleaned_receptor.pdb"
    script_path = output_dir / "_pymol_cleanup.py"

    target_chain_str = target_chain.upper()
    cofactor_code = preserved_cofactor.upper() if preserved_cofactor else None

    script_content = textwrap.dedent(f"""\
        # Auto-generated PyMOL cleanup script
        import sys
        from pymol import cmd

        RECEPTOR_PDB = {json.dumps(str(receptor_path))}
        CLEANED_PDB  = {json.dumps(str(cleaned_pdb))}
        TARGET_CHAIN = {json.dumps(target_chain_str)}
        LIGAND_CODE  = {json.dumps(ligand_code.upper())}
        COFACTOR     = {json.dumps(cofactor_code) if cofactor_code else 'None'}

        # 1. Select the user-chosen target chain
        cmd.load(RECEPTOR_PDB, "receptor")
        cmd.select("selected_chain", f"chain {{TARGET_CHAIN}}")

        if cmd.count_atoms("selected_chain") == 0:
            print(f"ERROR: Chain '{{TARGET_CHAIN}}' not found in structure.")
            cmd.quit()
            sys.exit(1)

        print(f"Selected chain '{{TARGET_CHAIN}}'")

        # 2. Include cofactor in the preserved selection (may be on a
        #    different chain), then remove everything else.
        if COFACTOR:
            print(f"Including cofactor {{COFACTOR}} in preserved selection")
            cmd.select("selected_chain", f"selected_chain or (resn {{COFACTOR}})")

        cmd.remove("not selected_chain")

        # 3. Deselect the current selection
        cmd.deselect()

        # 4. Remove the identified ligand explicitly (by residue name) so
        #    it never contaminates the cleaned receptor regardless of whether
        #    the PDB uses ATOM or HETATM records for the ligand.
        cmd.remove(f"resn {{LIGAND_CODE}}")

        # 5. Remove all remaining heteroatoms except the preserved cofactor
        if COFACTOR:
            print(f"Preserving cofactor {{COFACTOR}} during cleanup")
            cmd.remove(f"hetatm and not resn {{COFACTOR}}")
        else:
            cmd.remove("hetatm")

        # 6. Save cleaned receptor — target-chain protein + preserved cofactor
        if COFACTOR:
            cmd.select(
                "cleaned_receptor",
                f"chain {{TARGET_CHAIN}} and (polymer.protein or resn {{COFACTOR}})",
            )
        else:
            cmd.select(
                "cleaned_receptor",
                f"chain {{TARGET_CHAIN}} and polymer.protein",
            )

        if cmd.count_atoms("cleaned_receptor") == 0:
            print(
                f"ERROR: No protein atoms remain for chain '{{TARGET_CHAIN}}' "
                f"after cleaning."
            )
            cmd.quit()
            sys.exit(1)

        cmd.save(CLEANED_PDB, "cleaned_receptor")
        print(
            f"Saved cleaned receptor: "
            f"{{cmd.count_atoms('cleaned_receptor')}} atoms -> {{CLEANED_PDB}}"
        )

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
            ],
            ligand_code.upper(),
        )
    finally:
        if script_path.exists():
            script_path.unlink(missing_ok=True)

    return cleaned_pdb


# ════════════════════════════════════════════════════════════════════════════
#  Native Ligand Crystal Coordinates (Pure Python — no PyMOL)
# ════════════════════════════════════════════════════════════════════════════

def extract_native_ligand_coords(
    receptor_file: str | Path,
    ligand_code: str,
    target_chain: str,
) -> List[Tuple[float, float, float]]:
    """
    Return the heavy-atom crystal coordinates of the selected native ligand,
    parsed directly from the receptor structure (PDB or mmCIF).

    No intermediate PDB/SDF file is created: the grid box and the redock
    RMSD reference both use these in-memory coordinates.  This is the
    reliable replacement for the old PyMOL/PDB-based ligand extraction,
    which produced broken fragment files.

    Parameters
    ----------
    receptor_file : Path to the original (uncleaned) receptor file.
    ligand_code   : Three-letter residue name of the native ligand (e.g. "TCU").
    target_chain  : Chain ID to filter on (e.g. "A").

    Returns
    -------
    List of (x, y, z) tuples for heavy atoms of the matching HETATM ligand.
    Heavy atoms only (elements H/D excluded) — cofactor and waters are not
    included, so the grid stays centred on the ligand.

    Raises
    ------
    RuntimeError
        If no matching HETATM atoms are found in the receptor structure.
    """
    receptor_path = Path(receptor_file).resolve()

    target_ligand = ligand_code.upper()
    target_chain_str = target_chain.upper()

    # Group matching HETATM atoms by chain:  {chain: [(x, y, z), ...]}
    by_chain: Dict[str, List[Tuple[float, float, float]]] = {}
    for atom in iter_hetatm_atoms(receptor_path):
        if atom["resn"].upper() != target_ligand:
            continue
        if atom["element"] in ("H", "D"):
            continue
        by_chain.setdefault(atom["chain"] or "A", []).append(
            (atom["x"], atom["y"], atom["z"])
        )

    if not by_chain:
        raise RuntimeError(
            f"No HETATM atoms found for ligand '{target_ligand}' "
            f"on any chain in {receptor_path}.\n"
            "Cannot compute the docking grid box."
        )

    if target_chain_str in by_chain:
        coords = by_chain[target_chain_str]
    else:
        # Fallback: use the chain with the most matching atoms.
        best_chain = max(by_chain, key=lambda k: len(by_chain[k]))
        coords = by_chain[best_chain]
        logging.warning(
            f"Ligand '{target_ligand}' not found on chain '{target_chain_str}' "
            f"but found on chain(s) {list(by_chain.keys())}. "
            f"Using chain '{best_chain}' ({len(coords)} heavy atoms)."
        )

    return coords


# ════════════════════════════════════════════════════════════════════════════
#  Grid Box Calculation API
# ════════════════════════════════════════════════════════════════════════════

def calculate_grid_box(
    coords: List[Tuple[float, float, float]],
    padding: float = 4.0,
    min_size: float = 12.0,
    max_dimension: float = 25.0,
) -> Dict[str, float]:
    """
    Compute the AutoDock Vina docking grid box from native-ligand
    heavy-atom coordinates.

    Parameters
    ----------
    coords        : Native ligand heavy-atom coordinates (x, y, z) in Å,
                    e.g. from extract_native_ligand_coords().
    padding       : Extra padding per side in Å (default 4.0).
    min_size      : Minimum grid dimension in Å (default 12.0).
    max_dimension : Maximum grid dimension in Å (default 25.0).
                    Keeps volume ≤ 25³ = 15,625 Å³, well under
                    Vina's recommended 27,000 Å³ limit.

    Returns
    -------
    Dict with center_x/y/z, size_x/y/z, and a boolean '_clamped'
    flag indicating whether any axis was truncated.
    """
    if not coords:
        raise ValueError(
            "No native-ligand heavy-atom coordinates provided.\n"
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

    # Dimension = span + padding on both sides, clamped [min_size, max_dimension]
    sx = max(min(max_x - min_x + 2.0 * padding, max_dimension), min_size)
    sy = max(min(max_y - min_y + 2.0 * padding, max_dimension), min_size)
    sz = max(min(max_z - min_z + 2.0 * padding, max_dimension), min_size)

    clamped = (
        (max_x - min_x + 2.0 * padding > max_dimension)
        or (max_y - min_y + 2.0 * padding > max_dimension)
        or (max_z - min_z + 2.0 * padding > max_dimension)
    )

    return {
        "center_x": round(cx, 4),
        "center_y": round(cy, 4),
        "center_z": round(cz, 4),
        "size_x": round(sx, 2),
        "size_y": round(sy, 2),
        "size_z": round(sz, 2),
        "_clamped": clamped,
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
                f"Ensure ligand code '{ligand_code}' exists in the receptor structure."
            )
        if path.stat().st_size == 0:
            raise RuntimeError(
                f"PyMOL produced an empty file for {label}: {path}"
            )