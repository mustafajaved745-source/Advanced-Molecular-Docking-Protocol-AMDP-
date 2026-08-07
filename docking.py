"""
Molecular Docking Pipeline
=====================================
docking.py — AutoDock Vina execution, RMSD validation, summary report

Responsibilities:
  1. Run AutoDock Vina for any receptor + ligand PDBQT pair
  2. Parse binding affinity scores from Vina stdout / log file
  3. Calculate redocking RMSD between crystal pose and top docked pose
       - Heavy atoms only
       - Hungarian-algorithm optimal atom matching (scipy) with
         direct RMSD fallback when scipy is unavailable
  4. Write a formatted plain-text summary report

Public API
----------
  run_vina_docking(receptor_pdbqt, ligand_pdbqt, grid, output_dir,
                   prefix, vina_exe, exhaustiveness, num_modes, cpu, log_cb)
      -> dict  {best_score, all_scores, output_pdbqt}

  calculate_redock_rmsd(native_coords, redocked_pdbqt)
      -> float  (Angstroms; 999.0 on failure)

  generate_summary_report(output_dir, receptor_file, ligand_files,
                           active_site_code, selection_reason, grid,
                           redock_rmsd, redock_score, docking_results)
      -> Path  (report file)
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple


# ════════════════════════════════════════════════════════════════════════════
#  Vina docking
# ════════════════════════════════════════════════════════════════════════════

def run_vina_docking(
    receptor_pdbqt: Path,
    ligand_pdbqt:   Path,
    grid:           Dict[str, float],
    output_dir:     Path,
    prefix:         str,
    vina_exe:       str,
    exhaustiveness: int = 8,
    num_modes:      int = 9,
    cpu:            int = 0,
    log_cb:         Callable[[str, str], None] | None = None,
) -> Dict:
    """
    Run AutoDock Vina for one receptor–ligand pair.

    Parameters
    ----------
    receptor_pdbqt : Prepared receptor PDBQT.
    ligand_pdbqt   : Prepared ligand PDBQT.
    grid           : Dict with center_x/y/z and size_x/y/z (Angstroms).
    output_dir     : Folder where poses and log are written.
    prefix         : Filename stem for output files.
    vina_exe       : Path to the vina binary or PATH command name.
    exhaustiveness : Search exhaustiveness (default 8).
    num_modes      : Maximum poses to generate (default 9).
    cpu            : CPU cores to use; 0 = Vina auto-detect.
    log_cb         : Optional GUI log callback.

    Returns
    -------
    dict with keys:
        best_score   : float  — top binding affinity (kcal/mol)
        all_scores   : list   — all pose scores in rank order
        output_pdbqt : Path   — multi-model PDBQT with all poses
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    receptor_pdbqt = Path(receptor_pdbqt).resolve()
    ligand_pdbqt   = Path(ligand_pdbqt).resolve()
    output_dir     = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_pdbqt = output_dir / f"{prefix}_poses.pdbqt"
    vina_log     = output_dir / f"{prefix}_vina.log"

    # Resolve Vina executable binary across platforms
    exe_path = _resolve_executable(vina_exe)

    # ── Build Vina command ────────────────────────────────────────────────
    cmd = [
        str(exe_path),
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--center_x", f"{grid['center_x']:.4f}",
        "--center_y", f"{grid['center_y']:.4f}",
        "--center_z", f"{grid['center_z']:.4f}",
        "--size_x",   f"{grid['size_x']:.2f}",
        "--size_y",   f"{grid['size_y']:.2f}",
        "--size_z",   f"{grid['size_z']:.2f}",
        "--out",      str(output_pdbqt),
        "--log",      str(vina_log),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(num_modes),
    ]
    if cpu > 0:
        cmd += ["--cpu", str(cpu)]

    log(f"  Vina: {exe_path.name}  exhaustiveness={exhaustiveness}"
        f"  modes={num_modes}  cpu={'auto' if cpu == 0 else cpu}", "dim")

    # ── Execute ───────────────────────────────────────────────────────────
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=900,   # 15-minute ceiling per ligand
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"AutoDock Vina executable not found at: {vina_exe}\n"
            "Checked PATH and common locations (conda, /usr/bin, /usr/local/bin).\n"
            "Install: conda install -c conda-forge autoDock-vina\n"
            "Then set the path in Settings (⚙)."
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"Vina timed out (>15 min) on: {ligand_pdbqt.name}\n"
            "Try reducing exhaustiveness, or check that the grid box is reasonable."
        )

    # ── Relay Vina output to log ──────────────────────────────────────────
    combined = (result.stdout or "") + (result.stderr or "")
    _relay_vina_output(combined, log)

    if result.returncode != 0:
        raise RuntimeError(
            f"AutoDock Vina failed (exit {result.returncode}) "
            f"on '{ligand_pdbqt.name}'.\n"
            f"stderr: {result.stderr[:600]}\n"
            f"{_vina_hint(result.stderr)}"
        )

    # ── Parse scores ──────────────────────────────────────────────────────
    scores = _parse_vina_scores(combined)

    # Fallback to reading the dedicated log file if stdout output lacked the table
    if not scores and vina_log.exists():
        scores = _parse_vina_scores(vina_log.read_text(encoding="utf-8", errors="replace"))

    if not scores:
        log("  Warning: could not parse binding scores from Vina output.", "warning")
        scores = [0.0]

    # ── Verify output file ────────────────────────────────────────────────
    if not output_pdbqt.exists() or output_pdbqt.stat().st_size == 0:
        raise RuntimeError(
            f"Vina ran but produced no output PDBQT for '{ligand_pdbqt.name}'.\n"
            "The ligand may have too many rotatable bonds, or the grid box\n"
            "may be too small to accommodate it."
        )

    return {
        "best_score":   scores[0],
        "all_scores":   scores,
        "output_pdbqt": output_pdbqt,
    }


# ════════════════════════════════════════════════════════════════════════════
#  Redocking RMSD
# ════════════════════════════════════════════════════════════════════════════

def calculate_redock_rmsd(
    native_coords: List[Tuple[float, float, float]],
    redocked_pdbqt: Path,
) -> float:
    """
    Calculate the RMSD between the crystal-pose native ligand and the
    top-ranked redocked pose. Heavy atoms only.

    ``native_coords`` are the native ligand's crystal heavy-atom
    coordinates (from receptor_prep.extract_native_ligand_coords), so no
    ligand file needs to be created from the receptor.

    Returns 999.0 if the calculation fails for any reason (missing file,
    parse error, atom count mismatch).
    """
    try:
        docked_coords = _parse_pdbqt_model1_heavy_coords(Path(redocked_pdbqt))

        if not native_coords or not docked_coords:
            return 999.0

        return _optimal_rmsd(native_coords, docked_coords)

    except Exception:
        return 999.0


# ════════════════════════════════════════════════════════════════════════════
#  Summary report
# ════════════════════════════════════════════════════════════════════════════

def generate_summary_report(
    output_dir:      Path,
    receptor_file:   str,
    ligand_files:    List[str],
    active_site_code: str,
    target_chain:    str = "A",
    selection_reason: str = "",
    grid:            Dict[str, float] = None,
    redock_rmsd:     Optional[float] = None,
    redock_score:    float = 0.0,
    docking_results: List[Dict] = None,
    grid_padding:    float = 8.0,
) -> Path:
    """
    Write a plain-text summary report to output_dir/docking_summary.txt.

    Returns the Path to the written file.
    """
    output_dir  = Path(output_dir)
    report_path = output_dir / "docking_summary.txt"
    sorted_r    = sorted(docking_results, key=lambda x: x.get("best_score", 0.0))

    # ── RMSD status string ────────────────────────────────────────────────
    if redock_rmsd is None:
        rmsd_status = "SKIPPED  ─  (no valid 3D native ligand available)"
    elif redock_rmsd <= 2.0:
        rmsd_status = f"PASS  ✓  ({redock_rmsd:.2f} Å ≤ 2.0 Å threshold)"
    elif redock_rmsd <= 3.0:
        rmsd_status = f"MARGINAL  ⚠  ({redock_rmsd:.2f} Å — interpret with caution)"
    else:
        rmsd_status = f"FAIL  ✗  ({redock_rmsd:.2f} Å > 3.0 Å — check pocket centring)"

    # ── Column width for ligand name ──────────────────────────────────────
    col = max((len(r.get("name", "Unknown")) for r in sorted_r), default=10) + 2

    lines: List[str] = [
        "=" * 72,
        "  MOLECULAR DOCKING PIPELINE — SUMMARY REPORT",
        f"  Generated : {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}",
        "=" * 72,
        "",
        "━━  INPUT FILES",
        f"  Receptor : {receptor_file}",
        f"  Ligands  : {len(ligand_files)} file(s)",
        *[f"    • {f}" for f in ligand_files],
        "",
        "━━  ACTIVE-SITE SELECTION",
        f"  Identified ligand : {active_site_code}",
        f"  Target chain      : {target_chain}",
        f"  Selection         : {selection_reason}",
        "",
        "━━  DOCKING GRID BOX",
        f"  Centre     : ({grid['center_x']:.3f},  {grid['center_y']:.3f},  "
        f"{grid['center_z']:.3f})  Å",
        f"  Dimensions : {grid['size_x']:.1f}  ×  {grid['size_y']:.1f}  ×  "
        f"{grid['size_z']:.1f}  Å",
        f"  Padding    : {grid_padding:.1f} Å per side ({grid_padding * 2:.1f} Å total per axis)",
        "",
        "━━  VALIDATION — NATIVE LIGAND REDOCKING",
        f"  Redock score : {redock_score:.2f} kcal/mol"
        if redock_score and redock_rmsd is not None
        else "  Redock score : SKIPPED",
        f"  RMSD status  : {rmsd_status}",
        "",
        "━━  PRODUCTION DOCKING RESULTS  (ranked by binding affinity)",
        "",
        f"  {'Rank':<5}  {'Ligand':<{col}}  {'Best (kcal/mol)':>16}  "
        f"{'All poses (kcal/mol)'}",
        f"  {'─'*5}  {'─'*col}  {'─'*16}  {'─'*30}",
    ]

    for rank, r in enumerate(sorted_r, 1):
        name = r.get("name", "Unknown")
        best = r.get("best_score", 0.0)
        all_s = r.get("all_scores", [best])
        scores_str = ",  ".join(f"{s:.2f}" for s in all_s[:5])
        if len(all_s) > 5:
            scores_str += " …"
        lines.append(
            f"  {rank:<5}  {name:<{col}}  {best:>16.2f}  "
            f"[{scores_str}]"
        )

    lines += [
        "",
        "━━  OUTPUT DIRECTORY",
        f"  {output_dir}",
        f"  {'─' * 60}",
        "  cleaned_receptor.pdb       protein-only structure",
        "  cleaned_receptor.pdbqt     Meeko-prepared receptor",
        "  grid_box.txt               Vina grid configuration",
        "  redocking/                 fetched native ligand + redock poses",
        "  ligands/                   prepared ligand PDBQT files",
        "  docking/                   docked poses per ligand",
        "  docking_summary.txt        this report",
        "  MDP-log.txt                full pipeline log",
        "",
        "=" * 72,
    ]

    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


# ════════════════════════════════════════════════════════════════════════════
#  Internal Helpers & Utilities
# ════════════════════════════════════════════════════════════════════════════

def _resolve_executable(exe: str) -> Path:
    """Resolve executable path or look it up in PATH and common locations.

    Works cross-platform:
      - Checks path as-is
      - On Windows, auto-appends ``.exe`` if the file is found that way
      - Searches PATH via ``shutil.which``
      - Falls back to platform-specific common install locations
    """
    p = Path(exe)

    # 1. Check if the path as-is points to an existing file
    if p.is_file():
        return p.resolve()

    # 2. On Windows, try appending .exe when no extension is present
    if os.name == "nt" and not p.suffix:
        p_exe = p.with_suffix(".exe")
        if p_exe.is_file():
            return p_exe.resolve()

    # 3. Look up in PATH (handles both "vina" and "vina.exe")
    found = shutil.which(exe)
    if not found and os.name == "nt" and not exe.lower().endswith(".exe"):
        found = shutil.which(exe + ".exe")
    if found:
        return Path(found).resolve()

    # 4. Platform-specific common install locations
    if os.name == "nt":
        common_paths = [
            Path(r"C:\Program Files\AutoDock Vina\vina.exe"),
            Path(r"C:\Program Files (x86)\AutoDock Vina\vina.exe"),
            Path.home() / "AppData" / "Local" / "Programs" / "AutoDock Vina" / "vina.exe",
            Path.home() / "miniconda3" / "Scripts" / "vina.exe",
            Path.home() / "anaconda3" / "Scripts" / "vina.exe",
            Path.home() / "miniforge3" / "Scripts" / "vina.exe",
            Path.home() / "mambaforge" / "Scripts" / "vina.exe",
        ]
    else:
        common_paths = [
            Path("/usr/bin/vina"),
            Path("/usr/local/bin/vina"),
            Path.home() / "miniconda3" / "bin" / "vina",
            Path.home() / "anaconda3" / "bin" / "vina",
            Path.home() / "miniforge3" / "bin" / "vina",
            Path.home() / "mambaforge" / "bin" / "vina",
        ]

    for candidate in common_paths:
        if candidate.is_file():
            return candidate.resolve()

    return p


def _parse_vina_scores(text: str) -> List[float]:
    """Extract binding affinity scores from Vina output or log files."""
    scores: List[float] = []
    in_table = False

    for line in text.splitlines():
        stripped = line.strip()

        if "-----+------------" in stripped or "mode |" in stripped.lower():
            in_table = True
            continue

        if in_table:
            if not stripped:
                break
            parts = stripped.split()
            if parts and parts[0].isdigit():
                try:
                    scores.append(float(parts[1]))
                except (IndexError, ValueError):
                    pass
            elif scores:
                break

    return scores


def _parse_pdbqt_model1_heavy_coords(
    pdbqt_file: Path,
) -> List[Tuple[float, float, float]]:
    """Parse heavy-atom coordinates from MODEL 1 of a PDBQT structure."""
    coords: List[Tuple[float, float, float]] = []

    if not pdbqt_file.is_file():
        return coords

    try:
        text  = pdbqt_file.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
    except OSError:
        return coords

    has_models   = any(l.startswith("MODEL") for l in lines)
    in_target    = not has_models

    for line in lines:
        if line.startswith("MODEL"):
            model_num = line.split()[-1].strip() if line.split() else ""
            in_target = (model_num == "1")
            continue

        if line.startswith("ENDMDL"):
            if in_target:
                break
            in_target = False
            continue

        if not in_target:
            continue

        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue

        try:
            ad_type   = line[77:79].strip().upper() if len(line) >= 79 else ""
            atom_name = line[12:16].strip()

            if ad_type in ("H", "HD", "HS"):
                continue
            if not ad_type:
                first = "".join(c for c in atom_name if not c.isdigit())[:1].upper()
                if first == "H":
                    continue

            coords.append((
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ))
        except (ValueError, IndexError):
            continue

    return coords


def _optimal_rmsd(
    coords1: List[Tuple[float, float, float]],
    coords2: List[Tuple[float, float, float]],
) -> float:
    """
    Compute heavy-atom RMSD with Hungarian algorithm matching fallback.
    """
    n = min(len(coords1), len(coords2))
    if n == 0:
        return 999.0

    try:
        import numpy as np
        from scipy.spatial.distance import cdist
        from scipy.optimize import linear_sum_assignment

        a1 = np.array(coords1, dtype=float)
        a2 = np.array(coords2, dtype=float)

        D = cdist(a1, a2)
        row_ind, col_ind = linear_sum_assignment(D)

        # Slice to maximum available matched atom pairs
        row_ind = row_ind[:n]
        col_ind = col_ind[:n]

        sq_sum = float(sum(D[r, c] ** 2 for r, c in zip(row_ind, col_ind)))
        return math.sqrt(sq_sum / n)

    except ImportError:
        sq_sum = sum(
            (c1[0] - c2[0]) ** 2 +
            (c1[1] - c2[1]) ** 2 +
            (c1[2] - c2[2]) ** 2
            for c1, c2 in zip(coords1[:n], coords2[:n])
        )
        return math.sqrt(sq_sum / n)


def _relay_vina_output(
    text: str,
    log:  Callable[[str, str], None],
) -> None:
    """Relay stdout and stderr from process executions back to the log pump."""
    skip_prefixes = ("*", "=", "Detected")
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if any(line.startswith(p) for p in skip_prefixes):
            log(f"  Vina > {line}", "dim")
            continue
        upper = line.upper()
        if "ERROR" in upper or "FAILED" in upper:
            log(f"  Vina > {line}", "error")
        elif "WARNING" in upper:
            log(f"  Vina > {line}", "warning")
        elif line.lstrip().split() and line.lstrip().split()[0].isdigit():
            log(f"  Vina > {line}", "success")
        else:
            log(f"  Vina > {line}", "dim")


def _vina_hint(stderr: str) -> str:
    """Return context-aware hints for AutoDock Vina errors."""
    s = (stderr or "").lower()

    if "could not open" in s or "no such file" in s:
        return (
            "Hint: One of the input files could not be opened.\n"
            "  Check that the receptor and ligand PDBQT files were created\n"
            "  successfully in the previous preparation steps."
        )
    if "invalid grid" in s or "grid" in s:
        return (
            "Hint: The grid box parameters may be invalid.\n"
            "  Verify that the target pocket center and size dimensions are non-zero."
        )
    if "too many" in s and "torsion" in s:
        return (
            "Hint: Ligand exceeds maximum rotatable bond limit (>32).\n"
            "  Pre-treat or freeze rigid bonds prior to docking."
        )
    if "out of memory" in s or "bad alloc" in s:
        return (
            "Hint: Vina ran out of memory.\n"
            "  Reduce the grid box volume or lower the exhaustiveness setting."
        )

    return (
        "Hint: Ensure AutoDock Vina >= 1.2 is installed and the executable path\n"
        "is correctly configured in pipeline parameters."
    )