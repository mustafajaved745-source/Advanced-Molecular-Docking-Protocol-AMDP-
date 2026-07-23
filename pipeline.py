"""
AI-Guided Molecular Docking Pipeline
=====================================
pipeline.py — Stage orchestrator

Coordinates every step in order, calling the specialist modules:
  Stage 0 ─ Output directory setup
  Stage 1 ─ AI active-site identification   (ai_identify.py)
  Stage 2 ─ Receptor cleaning + grid box    (receptor_prep.py)
  Stage 3 ─ Receptor PDBQT preparation      (ligand_prep.py)
  Stage 4 ─ Native ligand PDBQT (redock)    (ligand_prep.py)
  Stage 5 ─ Input ligand PDBQT preparation  (ligand_prep.py)
  Stage 6 ─ Validation redocking + RMSD     (docking.py)
  Stage 7 ─ Production docking              (docking.py)
  Stage 8 ─ Summary report                  (docking.py)

All log / progress / status updates are pushed through callbacks so the
GUI thread is never touched directly from here.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from rdkit import Chem

from ai_identify import identify_active_site_ligand, parse_hetatm_records
from docking import calculate_redock_rmsd, generate_summary_report, run_vina_docking
from ligand_prep import prepare_ligand_pdbqt, prepare_receptor_pdbqt
from receptor_prep import calculate_grid_box, clean_receptor_with_pymol

# Standard user-agent header for web API requests
_HTTP_HEADERS: dict[str, str] = {
    "User-Agent": "AIGuidedDockingPipeline/1.0 (Academic Research Pipeline)"
}


# ════════════════════════════════════════════════════════════════════════════
#  Online Metadata & Structure Fetching
# ════════════════════════════════════════════════════════════════════════════

def fetch_native_ligand_sdf(
    ligand_code: str,
    output_path: Path,
    log_cb: Optional[Callable[[str, str], None]] = None,
    timeout: int = 15,
) -> bool:
    """
    Fetch clean 3D SDF for a ligand code from RCSB or PubChem,
    ensure explicit hydrogens, and save to output_path.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    code_upper = ligand_code.upper()
    log(f"Fetching clean 3D SDF for ligand [{code_upper}] from RCSB PDB...", "info")

    tmp_path = output_path.with_suffix(".tmp.sdf")

    # 1. Try RCSB PDB Ligand Expo API first
    rcsb_url = f"https://files.rcsb.org/ligands/download/{code_upper}.sdf"
    if _download_and_process_sdf(rcsb_url, tmp_path, output_path, timeout):
        log(f"Successfully retrieved 3D structure for [{code_upper}] from RCSB.", "success")
        return True
    else:
        log(f"RCSB fetch failed for [{code_upper}].", "dim")

    # 2. Fallback: Query PubChem API by ligand code
    log(f"Attempting PubChem lookup for [{code_upper}]...", "info")
    pubchem_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(code_upper)}/SDF?record_type=3d"
    )
    if _download_and_process_sdf(pubchem_url, tmp_path, output_path, timeout):
        log(f"Successfully retrieved 3D structure for [{code_upper}] from PubChem.", "success")
        return True
    else:
        log(f"PubChem fetch failed for [{code_upper}].", "warning")

    return False


def _download_and_process_sdf(
    url: str,
    tmp_path: Path,
    final_path: Path,
    timeout: int = 15,
) -> bool:
    """Download an SDF from a URL and sanitize explicit 3D hydrogens via RDKit."""
    req = urllib.request.Request(url, headers=_HTTP_HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            tmp_path.write_bytes(response.read())

        mol = Chem.MolFromMolFile(str(tmp_path), removeHs=False)
        if mol is not None:
            mol = Chem.AddHs(mol, addCoords=True)
            with Chem.SDWriter(str(final_path)) as writer:
                writer.write(mol)
            return True
    except (urllib.error.URLError, TimeoutError, ValueError, Exception):
        pass
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    return False


# ════════════════════════════════════════════════════════════════════════════
#  Pipeline Execution Entry Point
# ════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    receptor_file: str,
    ligand_files: List[str],
    config: Dict[str, Any],
    log_cb: Callable[[str, str], None],
    progress_cb: Callable[[float], None],
    status_cb: Callable[[str], None],
    stop_flag: Callable[[], bool],
) -> None:
    """
    Execute the full docking pipeline.

    Parameters
    ----------
    receptor_file : Path to the raw receptor PDB.
    ligand_files  : List of raw ligand SDF paths.
    config        : Settings dict (from settings.load_config()).
    log_cb        : log_cb(message, level) — level in
                    {info, success, warning, error, header, dim, bold}
    progress_cb   : progress_cb(0–100 float)
    status_cb     : status_cb(short string for status display)
    stop_flag     : Returns True when the pipeline execution is cancelled.
    """

    def log(msg: str, lvl: str = "info") -> None:
        log_cb(msg, lvl)

    def banner(title: str) -> None:
        log("─" * 58, "header")
        log(f"  {title}", "header")
        log("─" * 58, "header")

    def check_stop() -> bool:
        if stop_flag():
            log("Pipeline halted by user request.", "warning")
            return True
        return False

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 0 — Setup
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 0  ─  Setup")
    status_cb("Setting up output directory…")

    output_dir = _make_output_dir(
        config.get("output_dir", str(Path.home() / "docking_results"))
    )
    log(f"Output directory : {output_dir}", "info")
    log(f"Receptor         : {Path(receptor_file).name}", "info")
    log(f"Ligands          : {len(ligand_files)} file(s)", "info")
    for lf in ligand_files:
        log(f"   • {Path(lf).name}", "dim")
    progress_cb(3.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 1 — AI Active-Site Identification
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 1  ─  AI Active-Site Identification")
    status_cb("Parsing receptor HETATM records…")

    log("Scanning receptor PDB for HETATM records…", "info")
    hetatms = parse_hetatm_records(receptor_file)

    if not hetatms:
        raise RuntimeError(
            "No HETATM records found in the receptor PDB.\n"
            "The protein must contain a co-crystallised ligand for active-site\n"
            "identification. Waters-only / apo structures are not supported."
        )

    candidate_names = [h["resn"] for h in hetatms if not h["is_known_non_ligand"]]
    log(
        f"Found {len(hetatms)} HETATM group(s); "
        f"{len(candidate_names)} candidate ligand(s): "
        + (", ".join(candidate_names) if candidate_names else "none"),
        "info",
    )

    status_cb("Contacting AI model API…")
    log(
        f"Sending HETATM summary to {config.get('api_provider', 'API')} / "
        f"{config.get('model', 'unknown model')}…",
        "info",
    )

    ligand_code, ai_reason = identify_active_site_ligand(
        hetatms=hetatms,
        api_key=config["api_key"],
        model=config.get("model", "meta-llama/llama-3.1-70b-instruct:free"),
        base_url=config.get("api_base_url", "https://openrouter.ai/api/v1"),
    )

    # Fetch chemical metadata
    ligand_name, pubchem_cid = _fetch_rcsb_chem_metadata(ligand_code)

    log(f"✔  Active-site ligand identified : [{ligand_code}]", "success")
    log(f"   Chemical Name : {ligand_name}", "info")
    log(f"   PubChem CID   : {pubchem_cid}", "info")
    log(f"   Reasoning     : {ai_reason}", "dim")
    progress_cb(15.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 2 — Receptor Cleaning (PyMOL) & Grid Box
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 2  ─  Receptor Cleaning  (PyMOL)")
    status_cb("Cleaning receptor with PyMOL…")

    cleaned_pdb, native_lig_pdb, native_lig_sdf = clean_receptor_with_pymol(
        receptor_pdb=receptor_file,
        ligand_code=ligand_code,
        output_dir=output_dir,
        pymol_exe=config["pymol_exe"],
        log_cb=log,
    )

    log(f"✔  Cleaned receptor  : {cleaned_pdb.name}", "success")
    log(f"✔  Native ligand PDB : {native_lig_pdb.name}", "success")
    log(f"✔  Native ligand SDF : {native_lig_sdf.name}", "success")
    progress_cb(28.0)

    if check_stop():
        return

    banner("STAGE 2b ─  Docking Grid Box Calculation")
    status_cb("Calculating docking grid…")

    grid = calculate_grid_box(native_lig_pdb, padding=8.0)

    log(
        f"Grid centre     : "
        f"({grid['center_x']:.3f},  {grid['center_y']:.3f},  {grid['center_z']:.3f})  Å",
        "info",
    )
    log(
        f"Grid dimensions : "
        f"{grid['size_x']:.1f}  ×  {grid['size_y']:.1f}  ×  {grid['size_z']:.1f}  Å",
        "info",
    )

    grid_cfg_path = output_dir / "grid_box.txt"
    _write_grid_config(grid_cfg_path, cleaned_pdb, grid)
    log(f"Grid config saved : {grid_cfg_path.name}", "dim")
    progress_cb(33.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 3 — Receptor PDBQT Preparation  (Meeko)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 3  ─  Receptor PDBQT Preparation  (Meeko)")
    status_cb("Preparing receptor PDBQT…")

    receptor_pdbqt = prepare_receptor_pdbqt(
        pdb_file=cleaned_pdb,
        output_dir=output_dir,
        cmd=config["mk_prepare_receptor_cmd"],
        log_cb=log,
    )

    log(f"✔  Receptor PDBQT : {receptor_pdbqt.name}", "success")
    progress_cb(40.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 4 — Native Ligand PDBQT Preparation
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 4  ─  Native Ligand PDBQT Preparation")
    log(
        "Native ligand preparation deferred to Stage 6 (online clean 3D structure retrieval).",
        "dim",
    )
    progress_cb(42.0)

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 5 — Input Ligand Preparation  (Meeko)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 5  ─  Input Ligand Preparation  (Meeko)")
    status_cb("Preparing input ligand PDBQT files…")

    ligand_dir = output_dir / "ligands"
    ligand_dir.mkdir(exist_ok=True)

    prepared_ligands: list[tuple[str, Path]] = []
    n_ligs = len(ligand_files)

    for i, lig_file in enumerate(ligand_files):
        if check_stop():
            return
        name = Path(lig_file).stem
        log(f"[{i+1}/{n_ligs}]  Preparing : {name}", "info")
        pdbqt = prepare_ligand_pdbqt(
            sdf_file=Path(lig_file),
            output_dir=ligand_dir,
            cmd=config["mk_prepare_ligand_cmd"],
            prefix=name,
            log_cb=log,
        )
        prepared_ligands.append((name, pdbqt))
        progress_cb(42.0 + (i + 1) * (13.0 / n_ligs))

    log(f"✔  {len(prepared_ligands)} ligand(s) prepared.", "success")

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 6 — Validation Redocking
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 6  ─  Validation Redocking")
    status_cb("Running validation redocking…")

    redock_dir = output_dir / "redocking"
    redock_dir.mkdir(parents=True, exist_ok=True)

    fetched_sdf = redock_dir / f"{ligand_code}_online.sdf"
    native_pdbqt: Optional[Path] = None
    redock_score: float = 0.0
    rmsd: Optional[float] = None

    if fetch_native_ligand_sdf(ligand_code, fetched_sdf, log_cb=log):
        try:
            native_pdbqt = prepare_ligand_pdbqt(
                sdf_file=fetched_sdf,
                output_dir=redock_dir,
                cmd=config["mk_prepare_ligand_cmd"],
                prefix=f"{ligand_code}_native",
                log_cb=log,
            )
        except Exception as e:
            log(f"Failed to prepare online SDF into PDBQT: {e}", "warning")

    if native_pdbqt and native_pdbqt.exists():
        log(f"Re-docking native ligand [{ligand_code}] ({ligand_name}) back into pocket…", "info")

        redock_result = run_vina_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=native_pdbqt,
            grid=grid,
            output_dir=redock_dir,
            prefix="native_redock",
            vina_exe=config["vina_exe"],
            exhaustiveness=int(config.get("exhaustiveness", 8)),
            num_modes=int(config.get("num_modes", 9)),
            cpu=int(config.get("cpu", 0)),
            log_cb=log,
        )

        redock_score = redock_result["best_score"]
        log(f"   Top redock score : {redock_score:.2f} kcal/mol", "info")

        try:
            rmsd = calculate_redock_rmsd(native_lig_pdb, redock_result["output_pdbqt"])
            _report_rmsd(rmsd, log)
        except Exception as e:
            log(f"   RMSD calculation skipped: {e}", "warning")
    else:
        log("Online native ligand preparation failed — skipping validation redocking.", "warning")

    progress_cb(65.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 7 — Production Docking
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 7  ─  Production Docking")
    status_cb("Docking ligands…")

    dock_dir = output_dir / "docking"
    dock_dir.mkdir(exist_ok=True)

    docking_results: list[dict[str, Any]] = []
    n = len(prepared_ligands)

    for i, (name, pdbqt) in enumerate(prepared_ligands):
        if check_stop():
            break

        log("", "dim")
        log(f"[{i+1}/{n}]  Docking : {name}", "bold")

        result = run_vina_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=pdbqt,
            grid=grid,
            output_dir=dock_dir,
            prefix=name,
            vina_exe=config["vina_exe"],
            exhaustiveness=int(config.get("exhaustiveness", 8)),
            num_modes=int(config.get("num_modes", 9)),
            cpu=int(config.get("cpu", 0)),
            log_cb=log,
        )

        docking_results.append(
            {
                "name": name,
                "best_score": result["best_score"],
                "all_scores": result["all_scores"],
                "output_file": result["output_pdbqt"],
            }
        )

        log(f"   ✔  Best score : {result['best_score']:.2f} kcal/mol", "success")
        progress_cb(65.0 + (i + 1) * (30.0 / n))

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 8 — Summary Report
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 8  ─  Summary Report")
    status_cb("Generating summary report…")

    report_path = generate_summary_report(
        output_dir=output_dir,
        receptor_file=receptor_file,
        ligand_files=ligand_files,
        active_site_code=ligand_code,
        ai_reason=ai_reason,
        grid=grid,
        redock_rmsd=rmsd,
        redock_score=redock_score,
        docking_results=docking_results,
    )

    log(f"✔  Report written : {report_path.name}", "success")

    # Final ranked table output
    log("", "dim")
    log("─" * 58, "header")
    log("  RANKED DOCKING RESULTS", "header")
    log("─" * 58, "header")

    sorted_r = sorted(docking_results, key=lambda x: x["best_score"])
    col_w = max((len(r["name"]) for r in sorted_r), default=10) + 2

    log(f"  {'Rank':<5}  {'Ligand':<{col_w}}  {'Best (kcal/mol)':>15}", "dim")
    log(f"  {'─'*5}  {'─'*col_w}  {'─'*15}", "dim")

    for rank, r in enumerate(sorted_r, 1):
        log(
            f"  {rank:<5}  {r['name']:<{col_w}}  {r['best_score']:>15.2f}",
            "success",
        )

    log("", "dim")
    log("✔  All output saved to :", "success")
    log(f"   {output_dir}", "bold")
    progress_cb(100.0)
    status_cb("Pipeline complete ✓")


# ════════════════════════════════════════════════════════════════════════════
#  Internal Helper Functions
# ════════════════════════════════════════════════════════════════════════════

def _fetch_rcsb_chem_metadata(ligand_code: str) -> Tuple[str, str]:
    """Query RCSB Data REST API for official chemical name and PubChem CID."""
    ligand_name = ligand_code
    pubchem_cid = "N/A"
    meta_url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{urllib.parse.quote(ligand_code.upper())}"
    req = urllib.request.Request(meta_url, headers=_HTTP_HEADERS)

    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            chem_comp = data.get("chem_comp", {})
            ligand_name = chem_comp.get("name", ligand_code)

            identifiers = data.get("rcsb_chem_comp_identifiers", {}).get("identifiers", [])
            for ident in identifiers:
                if ident.get("program") == "PubChem":
                    pubchem_cid = str(ident.get("identifier"))
                    break
    except Exception:
        pass

    return ligand_name, pubchem_cid


def _make_output_dir(base: str) -> Path:
    """Create and return a timestamped run directory."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path(base) / f"docking_run_{ts}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "redocking").mkdir(exist_ok=True)
    (out / "ligands").mkdir(exist_ok=True)
    (out / "docking").mkdir(exist_ok=True)
    return out


def _write_grid_config(path: Path, receptor_pdbqt: Path, grid: Dict[str, float]) -> None:
    """Write a Vina-compatible grid box configuration text file."""
    lines = [
        "# AutoDock Vina grid box configuration",
        f"# Generated by AI-Guided Docking Pipeline — {datetime.now():%Y-%m-%d %H:%M:%S}",
        "#",
        f"receptor  = {receptor_pdbqt}",
        f"center_x  = {grid['center_x']:.4f}",
        f"center_y  = {grid['center_y']:.4f}",
        f"center_z  = {grid['center_z']:.4f}",
        f"size_x    = {grid['size_x']:.2f}",
        f"size_y    = {grid['size_y']:.2f}",
        f"size_z    = {grid['size_z']:.2f}",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _report_rmsd(rmsd: float, log: Callable[[str, str], None]) -> None:
    """Log the redocking RMSD with colour-coded threshold evaluation."""
    if rmsd <= 2.0:
        log(f"✔  Redocking RMSD = {rmsd:.2f} Å  ←  PASS  (≤ 2.0 Å threshold)", "success")
        log(
            "   Docking pocket is correctly centred. Production results are reliable.",
            "dim",
        )
    elif rmsd <= 3.0:
        log(f"⚠  Redocking RMSD = {rmsd:.2f} Å  ←  MARGINAL  (2.0 – 3.0 Å)", "warning")
        log(
            "   Results may still be valid but interpret with caution.\n"
            "   Consider checking the identified ligand code or adjusting grid padding.",
            "warning",
        )
    else:
        log(f"✘  Redocking RMSD = {rmsd:.2f} Å  ←  FAIL  (> 3.0 Å)", "error")
        log(
            "   WARNING: The docking pocket may be mis-centred.\n"
            "   Verify the AI-identified ligand code is correct.\n"
            "   You may need to manually inspect the receptor and rerun.",
            "error",
        )