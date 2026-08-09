"""
Molecular Docking Pipeline
=====================================
pipeline.py — Stage orchestrator

Coordinates every step in order, calling the specialist modules:
  Stage 0 ─ Output directory setup
  Stage 1 ─ Active-site selection           (user-selected, deterministic)
  Stage 2 ─ Receptor cleaning + grid box    (receptor_prep.py)
  Stage 3 ─ Receptor PDBQT preparation      (ligand_prep.py)
  Stage 4 ─ Native ligand PDBQT (redock)    (ligand_prep.py)
  Stage 5 ─ Input ligand PDBQT preparation  (ligand_prep.py)
  Stage 6 ─ Validation redocking + RMSD     (docking.py)
  Stage 7 ─ Production docking              (docking.py)
  Stage 8 ─ PLIP interaction analysis       (interaction_analysis.py)
  Stage 9 ─ Summary report                  (docking.py)

All log / progress / status updates are pushed through callbacks so the
GUI thread is never touched directly from here.
"""

from __future__ import annotations

import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger

from hetatm_parser import parse_hetatm_records, build_docking_target_context
from docking import calculate_redock_rmsd, generate_summary_report, run_vina_docking
from interaction_analysis import combine_interaction_csvs, run_plip_interaction_analysis
from ligand_prep import ensure_3d_sdf, prepare_ligand_pdbqt, prepare_receptor_pdbqt
from llm_export import ExportOptions, export_experiment
from receptor_prep import (
    calculate_grid_box,
    clean_receptor_with_pymol,
    extract_native_ligand_coords,
)

# Standard user-agent header for web API requests.
# A browser-style UA is used because RCSB / PubChem / NCI occasionally
# reject or rate-limit bare library agents.
_HTTP_HEADERS: dict[str, str] = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
}


# ════════════════════════════════════════════════════════════════════════════
#  Online Metadata & Structure Fetching
# ════════════════════════════════════════════════════════════════════════════

def fetch_native_ligand_sdf(
    ligand_code: str,
    output_path: Path,
    ligand_name: Optional[str] = None,
    pubchem_cid: Optional[str] = None,
    log_cb: Optional[Callable[[str, str], None]] = None,
    timeout: int = 12,
) -> bool:
    """
    Fetch clean 3D SDF for a ligand code and save it to output_path.

    Identifier priority is CID → chemical name → residue code; source
    priority is PubChem first and RCSB last (RCSB's CCD SDFs are tagged 2D
    even though they carry 3D coordinates, which trips up downstream 3D prep).
    All lookups request explicit 3D coordinates so Meeko/Vina get a usable
    conformer.

    Fetch order:
      1. PubChem by CID (record_type=3d), then by SMILES via CACTUS when the
         ligand is a single complete molecule (disconnected salt/cofactor
         SMILES are skipped).
      2. PubChem by chemical name (record_type=3d), plus name→CID resolution.
      3. PubChem by residue code (record_type=3d).
      4. RCSB PDB — ideal → model → plain (last).
      5. NIH NCI/CADD (CACTUS) — by name → by code, final fallback.

    Returns True if any source produced a valid SDF.
    """
    def log(msg: str, lvl: str = "info") -> None:
        if log_cb:
            log_cb(msg, lvl)

    code_upper = ligand_code.strip().upper()
    tmp_path = output_path.with_suffix(".tmp.sdf")

    # Resolve the PubChem CID up front (param → chemical name → residue code)
    # so the most-specific identifier is always attempted first.
    cid = _normalize_cid(pubchem_cid)
    if not cid and ligand_name and ligand_name != code_upper:
        cid = _pubchem_name_to_cid(ligand_name, timeout=timeout)
    if not cid:
        cid = _pubchem_name_to_cid(code_upper, timeout=timeout)

    # ── 1. By CID (PubChem 3D), then by SMILES via CACTUS ───────────────────
    if cid:
        log(f"Attempting PubChem 3D lookup for CID [{cid}]...", "info")
        pubchem_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{urllib.parse.quote(cid, safe='')}/SDF?record_type=3d"
        )
        if _download_and_process_sdf(pubchem_url, tmp_path, output_path, timeout):
            log(f"Successfully retrieved 3D structure for [{code_upper}] from PubChem by CID.", "success")
            return True
        log(f"PubChem CID fetch failed for [{cid}].", "dim")

        # CID → SMILES → CACTUS 3D, only for a single complete molecule.
        smiles = _pubchem_cid_to_smiles(cid, timeout=timeout)
        if smiles and _is_complete_molecule(smiles):
            log("Attempting NIH NCI CACTUS 3D lookup via SMILES...", "info")
            cactus_url = (
                f"https://cactus.nci.nih.gov/chemical/structure/"
                f"{urllib.parse.quote(smiles, safe='')}/file?format=sdf&get3d=true"
            )
            if _download_and_process_sdf(cactus_url, tmp_path, output_path, timeout):
                log(f"Successfully retrieved 3D structure for [{code_upper}] via SMILES from NIH NCI CACTUS.", "success")
                return True
            log("NIH NCI CACTUS (SMILES) fetch failed.", "dim")

    # ── 2. By chemical name (PubChem) ───────────────────────────────────────
    if ligand_name and ligand_name != code_upper:
        log(f"Attempting PubChem 3D lookup for [{ligand_name}]...", "info")
        pubchem_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{urllib.parse.quote(ligand_name, safe='')}/SDF?record_type=3d"
        )
        if _download_and_process_sdf(pubchem_url, tmp_path, output_path, timeout):
            log(f"Successfully retrieved 3D structure for [{code_upper}] from PubChem by name.", "success")
            return True
        log(f"PubChem name fetch failed for [{ligand_name}].", "dim")

        # The strict name→SDF route often 404s on synonyms / brand names
        # (e.g. "TAXOL"). Resolve the name to a CID, then fetch 3D by CID.
        resolved_cid = _pubchem_name_to_cid(ligand_name, timeout=timeout)
        if resolved_cid and resolved_cid != cid:
            log(f"Attempting PubChem 3D lookup by resolved CID [{resolved_cid}]...", "info")
            pubchem_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{urllib.parse.quote(resolved_cid, safe='')}/SDF?record_type=3d"
            )
            if _download_and_process_sdf(pubchem_url, tmp_path, output_path, timeout):
                log(f"Successfully retrieved 3D structure for [{code_upper}] from PubChem by CID.", "success")
                return True
            log(f"PubChem resolved-CID fetch failed for [{resolved_cid}].", "dim")

    # ── 3. By residue code (PubChem) ────────────────────────────────────────
    log(f"Attempting PubChem 3D lookup for [{code_upper}]...", "info")
    pubchem_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(code_upper, safe='')}/SDF?record_type=3d"
    )
    if _download_and_process_sdf(pubchem_url, tmp_path, output_path, timeout):
        log(f"Successfully retrieved 3D structure for [{code_upper}] from PubChem.", "success")
        return True
    log(f"PubChem fetch failed for [{code_upper}].", "dim")

    # ── 4. RCSB PDB (last) ──────────────────────────────────────────────────
    log(f"Fetching ideal 3D SDF for ligand [{code_upper}] from RCSB PDB...", "info")
    rcsb_urls = [
        f"https://files.rcsb.org/ligands/view/{code_upper}_ideal.sdf",
        f"https://files.rcsb.org/ligands/download/{code_upper}_model.sdf",
        f"https://files.rcsb.org/ligands/download/{code_upper}.sdf",
    ]
    for rcsb_url in rcsb_urls:
        if _download_and_process_sdf(rcsb_url, tmp_path, output_path, timeout):
            log(f"Successfully retrieved 3D structure for [{code_upper}] from RCSB.", "success")
            return True
        log(f"RCSB fetch failed for [{code_upper}]: {rcsb_url}", "dim")

    # ── 5. NIH NCI/CADD (CACTUS) final fallback (name → code) ───────────────
    for name in (ligand_name, code_upper):
        if not name:
            continue
        log(f"Attempting NIH NCI CACTUS 3D lookup for [{name}]...", "info")
        cactus_url = (
            f"https://cactus.nci.nih.gov/chemical/structure/"
            f"{urllib.parse.quote(name, safe='')}/file?format=sdf&get3d=true"
        )
        if _download_and_process_sdf(cactus_url, tmp_path, output_path, timeout):
            log(f"Successfully retrieved 3D structure for [{name}] from NIH NCI CACTUS.", "success")
            return True
        log(f"NIH NCI CACTUS fetch failed for [{name}].", "dim")

    log(
        "⚠  All online 3D SDF fetchers failed (PubChem, RCSB PDB, NIH NCI "
        "CACTUS).\n"
        "   This is usually a temporary network / API issue; Stage 6 "
        "validation redocking will be skipped, production docking is unaffected.",
        "warning",
    )
    return False


def _download_and_process_sdf(
    url: str,
    tmp_path: Path,
    final_path: Path,
    timeout: int = 12,
    retries: int = 3,
) -> bool:
    """Download an SDF from a URL and sanitize explicit 3D hydrogens via RDKit.

    Retries up to *retries* times with exponential back-off (1 s, 2 s, 4 s)
    to tolerate transient failures — network blips and server-busy responses
    (HTTP 429 / 408 / 5xx), which PubChem especially is prone to. Deterministic
    failures — HTTP 4xx (404 / 403 / 410), or a payload that is not a parseable
    SDF — return False immediately without retrying.
    """
    for attempt in range(retries):
        req = urllib.request.Request(url, headers=_HTTP_HEADERS)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                tmp_path.write_bytes(response.read())

            # RCSB CCD SDFs are tagged 2D while carrying 3D coordinates, so
            # RDKit prints a benign "molecule is tagged as 2D… Marking the mol
            # as 3D" warning here. Quiet it — the molecule is written with 3D
            # coordinates and treated as 3D downstream regardless.
            RDLogger.DisableLog("rdApp.warning")
            try:
                mol = Chem.MolFromMolFile(str(tmp_path), removeHs=False)
                if mol is None:
                    # Non-SDF payload (e.g. an HTML error page) — nothing to salvage.
                    return False

                mol = Chem.AddHs(mol, addCoords=True)
                with Chem.SDWriter(str(final_path)) as writer:
                    writer.write(mol)
            finally:
                RDLogger.EnableLog("rdApp.warning")
            return True
        except urllib.error.HTTPError as e:
            # 5xx / 429 / 408 are transient server-side issues — retry with
            # back-off. Other 4xx (404/403/410) are deterministic — give up.
            if e.code >= 500 or e.code in (408, 429):
                pass
            else:
                return False
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            pass  # transient network failure → retry with back-off
        except Exception:
            pass  # RDKit parse/write hiccup → retry is harmless
        finally:
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)

        if attempt < retries - 1:
            time.sleep(1 * (2 ** attempt))  # 1 s, 2 s, 4 s

    return False


def _normalize_cid(value: Optional[str]) -> Optional[str]:
    """Return a clean PubChem CID string, or None for empty / 'N/A' / 'UNKNOWN'."""
    if value is None:
        return None
    s = str(value).strip()
    return None if s in ("N/A", "UNKNOWN", "") else s


def _pubchem_cid_to_smiles(cid: str, timeout: int = 12) -> Optional[str]:
    """Fetch a SMILES string for a PubChem CID via the PUG REST property
    endpoint.

    PubChem returns the requested property under a related key (an
    ``IsomericSMILES`` request comes back as ``SMILES``), so both possible
    keys are accepted. Returns None on any failure.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{urllib.parse.quote(cid, safe='')}/property/IsomericSMILES/JSON"
    )
    req = urllib.request.Request(url, headers=_HTTP_HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
        return props.get("SMILES") or props.get("ConnectivitySMILES")
    except Exception:
        return None


def _is_complete_molecule(smiles: str) -> bool:
    """True if the SMILES parses to a single, connected molecule.

    Disconnected SMILES (dots — e.g. salt / metal-cofactor complexes) are
    rejected so CACTUS is only asked for a complete small-molecule ligand.
    """
    if not smiles or "." in smiles:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return len(Chem.GetMolFrags(mol)) == 1
    except Exception:
        return False


def _pubchem_name_to_cid(name: str, timeout: int = 12) -> Optional[str]:
    """Resolve a chemical name to a PubChem CID via the PUG REST cids/JSON
    endpoint.

    The direct ``name/{name}/SDF`` route is strict about matching the exact
    primary name, whereas ``cids/JSON`` also matches synonyms / brand names
    (e.g. "TAXOL" → 36314). Returns the first CID, or None if unresolvable.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(name, safe='')}/cids/JSON"
    )
    req = urllib.request.Request(url, headers=_HTTP_HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        cids = data.get("IdentifierList", {}).get("CID", [])
        return str(cids[0]) if cids else None
    except Exception:
        return None


# Centroid-to-centroid distance (Å) within which a HETATM group is treated as
# an active-site neighbour of the selected ligand and therefore a candidate
# cofactor to preserve during receptor cleaning.  Chosen generous enough to
# catch elongated cofactors (e.g. NAD) whose centroid sits several Å from the
# ligand while still excluding the bulk of crystallographic waters / ions.
_COFACTOR_NEIGHBOUR_CUTOFF_A: float = 12.0

# Minimum heavy-atom count for a HETATM group to be considered a cofactor.
# Filters out single-atom ions and tiny fragments that happen to be near the
# active site (e.g. a lone crystallographic sulphate) while keeping real
# cofactors such as NAD / FAD / HEM / ATP / PLP.
_COFACTOR_MIN_HEAVY_ATOMS: int = 6


def _detect_preserved_cofactor(
    hetatms: List[Dict],
    ligand_code: str,
    target_chain: str,
) -> Optional[str]:
    """Deterministically pick a cofactor to preserve during receptor cleaning.

    The selected ligand is usually an inhibitor bound next to an essential
    catalytic cofactor (NAD, FAD, HEM, …).  Stripping that cofactor deforms
    the active-site pocket, so redocking validation can no longer reproduce
    the native pose (the TCU/NAD case: 3.84 Å FAIL without NAD vs 0.52 Å PASS
    with it).  This scans the parsed HETATM groups and returns the residue
    code of the most likely cofactor, or None.

    A group qualifies when it is
      * not the selected ligand itself,
      * not a known non-ligand (water / ion / cryoprotectant / buffer),
      * located on the target chain,
      * large enough to be a real cofactor (>= ``_COFACTOR_MIN_HEAVY_ATOMS``
        heavy atoms),
      * within ``_COFACTOR_NEIGHBOUR_CUTOFF_A`` Å of the ligand centroid.

    Among qualifying groups the closest one to the ligand wins.
    """
    ligand_centroid = None
    for h in hetatms:
        if h["resn"] == ligand_code and h["chain"] == target_chain:
            ligand_centroid = h.get("centroid")
            break
    if ligand_centroid is None:
        return None

    best: Optional[Dict] = None
    best_dist = float("inf")
    for h in hetatms:
        if h["resn"] == ligand_code and h["chain"] == target_chain:
            continue  # the ligand itself is removed by name, not preserved
        if h["is_known_non_ligand"]:
            continue
        if h["chain"] != target_chain:
            continue
        if h["heavy_atom_count"] < _COFACTOR_MIN_HEAVY_ATOMS:
            continue
        d = math.dist(ligand_centroid, h.get("centroid", (0.0, 0.0, 0.0)))
        if d > _COFACTOR_NEIGHBOUR_CUTOFF_A:
            continue
        if d < best_dist:
            best = h
            best_dist = d

    return best["resn"] if best else None


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
    write_log_fn: Optional[Callable[[Path], None]] = None,
    selected_chain: str = "",
    selected_ligand_code: str = "",
    ligand_status_cb: Optional[Callable[[str, str, str], None]] = None,
    selected_ligand_instance: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Execute the full docking pipeline.

    Parameters
    ----------
    receptor_file : Path to the raw receptor structure (PDB or mmCIF).
    ligand_files  : List of raw ligand SDF paths.
    config        : Settings dict (from settings.load_config()).
    log_cb        : log_cb(message, level) — level in
                    {info, success, warning, error, header, dim, bold}
    progress_cb   : progress_cb(0–100 float)
    status_cb     : status_cb(short string for status display)
    stop_flag     : Returns True when the pipeline execution is cancelled.
    write_log_fn  : Optional callback(write_path) to persist the full UI log
                    to a file after the pipeline completes.
    selected_chain      : Chain ID chosen by the user in the Active Site
                          Selection panel (e.g. "A").
    selected_ligand_code: 3-letter residue code of the ligand the user
                          chose (e.g. "ATP").  Stage 1 validates that a
                          selection was made; no auto-detection happens.
    ligand_status_cb : Optional callback(ligand_file_path, status, detail)
                       for live per-ligand status updates in the GUI.
                       status ∈ {preparing, prep_failed, docking, ok, failed}.
    selected_ligand_instance: Exact structural identity selected in the GUI.
                              This disambiguates repeated ligand codes within
                              one chain. Older callers may omit it when their
                              chain/code pair resolves to exactly one instance.
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

    # Record pipeline start time for elapsed-time reporting
    _pipeline_start = time.time()
    _experiment_created = datetime.now(timezone.utc).isoformat()

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
    log("", "dim")
    log("Configuration:", "dim")
    log(f"   Grid padding   : {config.get('grid_padding', 8.0)} Å", "dim")
    log(f"   Exhaustiveness : {config.get('exhaustiveness', 8)}", "dim")
    log(f"   Num modes      : {config.get('num_modes', 9)}", "dim")
    log(f"   CPU threads    : {'auto' if config.get('cpu', 0) == 0 else config.get('cpu', 0)}", "dim")
    progress_cb(3.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 1 — Active-Site Selection (user-confirmed, deterministic)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 1  ─  Active-Site Selection")
    status_cb("Verifying active-site selection…")

    if not selected_chain or not selected_ligand_code:
        raise RuntimeError(
            "No active site selected. Load a receptor structure (PDB or mmCIF) "
            "and select a chain and ligand before running."
        )

    ligand_code = selected_ligand_code.strip().upper()
    target_chain = selected_chain.strip().upper()

    log(
        f"Active site chosen by user : [{ligand_code}]  (chain {target_chain})",
        "info",
    )

    # Confirm the chosen chain + ligand really exist as a candidate HETATM
    # group in the receptor PDB (defence in depth — the GUI already guarantees
    # this, but a stale or edited file on disk should fail here instead of
    # three stages later).
    log("Confirming selection against receptor PDB HETATM records…", "dim")
    hetatms = parse_hetatm_records(receptor_file)
    if not hetatms:
        raise RuntimeError(
            "No HETATM records found in the receptor structure.\n"
            "The protein must contain a co-crystallised ligand for active-site "
            "selection."
        )

    if selected_ligand_instance:
        target_context = build_docking_target_context(
            receptor_file, selected_ligand_instance
        )
        selected_key = target_context.instance_key
        if (
            target_context.component_id.upper() != ligand_code
            or selected_key.chain_id.upper() != target_chain
        ):
            raise RuntimeError(
                "INVALID_LIGAND_INSTANCE: selected ligand instance does not match "
                f"[{ligand_code}] on chain [{target_chain}]. Re-load the receptor "
                "and select the active site again."
            )
    else:
        matches = [
            h for h in hetatms
            if h["resn"].upper() == ligand_code
            and h["chain"] == target_chain
            and not h["is_known_non_ligand"]
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"AMBIGUOUS_LIGAND_INSTANCE: [{ligand_code}] on chain [{target_chain}] "
                f"has {len(matches)} candidate instances; select an exact instance."
            )
        target_context = build_docking_target_context(
            receptor_file, matches[0]["instance_key"]
        )

    selection_reason = f"Manually selected by user: {ligand_code} in chain {target_chain}"

    # Detect an essential cofactor bound next to the ligand (NAD, FAD, HEM, …)
    # so it survives receptor cleaning.  Without it the active-site pocket is
    # incomplete and redocking validation cannot reproduce the native pose.
    preserved_cofactor = _detect_preserved_cofactor(hetatms, ligand_code, target_chain)
    if preserved_cofactor:
        log(
            f"✔  Preserving cofactor : {preserved_cofactor} "
            f"(bound beside {ligand_code})",
            "success",
        )
    else:
        log(
            "No active-site cofactor detected next to "
            f"{ligand_code} — cleaning to protein only.",
            "dim",
        )

    # Fetch chemical metadata — a plain lookup by the now-known code
    ligand_name, pubchem_cid = _fetch_rcsb_chem_metadata(ligand_code)

    # Warn about incomplete metadata
    if pubchem_cid in ("N/A", "UNKNOWN", ""):
        log(
            "⚠  Could not retrieve PubChem CID from any source.\n"
            "   Redocking validation may use a suboptimal 3D structure.",
            "warning",
        )
    if ligand_name == ligand_code:
        log(
            "⚠  Chemical name lookup failed; using ligand code as name.",
            "dim",
        )

    log(f"✔  Active-site ligand selected : [{ligand_code}]", "success")
    log(f"   Target Chain  : {target_chain}", "info")
    log(f"   Chemical Name : {ligand_name}", "info")
    log(f"   PubChem CID   : {pubchem_cid}", "info")
    log(f"   Selection     : {selection_reason}", "dim")
    progress_cb(15.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 2 — Receptor Cleaning (PyMOL) & Grid Box
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 2  ─  Receptor Cleaning  (PyMOL)")
    status_cb("Cleaning receptor with PyMOL…")
    log("Stripping water, ions, and non-protein heteroatoms via PyMOL…", "dim")

    # Deterministic: clean exactly the user-selected chain.  No auto-detection
    # and no chain fallback — a wrong selection fails loudly right here rather
    # than silently producing a mis-centred grid three stages later.
    cleaned_pdb = clean_receptor_with_pymol(
        receptor_path=receptor_file,
        ligand_code=ligand_code,
        target_chain=target_chain,
        output_dir=output_dir,
        pymol_exe=config["pymol_exe"],
        preserved_cofactor=preserved_cofactor,
        log_cb=log,
    )

    if not cleaned_pdb or not cleaned_pdb.exists():
        raise RuntimeError(
            "PyMOL failed to produce a valid cleaned receptor PDB.\n"
            "Check that the receptor file is a valid PDB with protein ATOM records."
        )

    log(f"✔  Cleaned receptor  : {cleaned_pdb.name}", "success")

    # Extract the native ligand's crystal coordinates directly from the
    # original receptor structure (pure Python, no PyMOL, no intermediate
    # PDB file) for grid-box calculation and RMSD reference.
    native_coords = extract_native_ligand_coords(
        receptor_file=receptor_file,
        ligand_code=ligand_code,
        target_chain=target_chain,
        ligand_instance=target_context.instance_key,
    )
    log(
        f"✔  Native ligand crystal coords : {len(native_coords)} heavy atoms",
        "success",
    )
    progress_cb(28.0)

    if check_stop():
        return

    banner("STAGE 2b ─  Docking Grid Box Calculation")
    status_cb("Calculating docking grid…")
    log("Computing grid from native ligand heavy-atom centroid…", "dim")

    grid_pad = float(config.get("grid_padding", 4.0))
    log(f"Grid padding : {grid_pad} Å per side (from settings)", "dim")

    grid = calculate_grid_box(native_coords, padding=grid_pad)

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

    # Warn if any axis exceeds Vina's effective operating limit
    oversized_axes = [
        axis for axis, val in
        [("X", grid["size_x"]), ("Y", grid["size_y"]), ("Z", grid["size_z"])]
        if val > 60.0
    ]
    if oversized_axes:
        log(
            f"⚠  Grid axis {', '.join(oversized_axes)} exceeds 60 Å. "
            f"AutoDock Vina may produce unreliable results for boxes this large.\n"
            f"   Consider narrowing the search space or reducing grid padding.",
            "warning",
        )

    # Warn if dimensions were clamped to max_dimension
    if grid.get("_clamped"):
        log(
            f"⚠  Grid dimension(s) clamped to 25.0 Å to keep search volume "
            f"≤ 15,625 Å³ (Vina best-practice limit).\n"
            f"   The binding site may extend beyond the grid — consider "
            f"reducing grid_padding in Settings.",
            "warning",
        )

    grid_cfg_path = output_dir / "grid_box.txt"
    _write_grid_config(grid_cfg_path, cleaned_pdb, grid)
    log(f"Grid config saved : {grid_cfg_path.name}", "dim")

    base_exhaustiveness = int(config.get("exhaustiveness", 8))
    exhaustiveness, grid_volume, pocket_label = _adaptive_exhaustiveness(
        base_exhaustiveness, grid
    )
    log(
        f"Grid volume : {grid_volume:,.0f} Å³  ({pocket_label} pocket) → "
        f"exhaustiveness {base_exhaustiveness} → {exhaustiveness}",
        "info" if exhaustiveness == base_exhaustiveness else "warning",
    )

    progress_cb(33.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 3 — Receptor PDBQT Preparation  (Meeko)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 3  ─  Receptor PDBQT Preparation  (Meeko)")
    status_cb("Preparing receptor PDBQT…")
    log(f"Converting {cleaned_pdb.name} → PDBQT via Meeko…", "dim")

    receptor_pdbqt = prepare_receptor_pdbqt(
        pdb_file=cleaned_pdb,
        output_dir=output_dir,
        cmd=config["mk_prepare_receptor_cmd"],
        target_chain=target_chain,
        grid=grid,
        log_cb=log,
    )

    log(f"✔  Receptor PDBQT : {receptor_pdbqt.name}", "success")
    progress_cb(40.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 4 — Native Ligand 3D Retrieval & PDBQT Preparation
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 4  ─  Native Ligand 3D Retrieval & PDBQT Preparation")
    status_cb("Fetching native ligand 3D structure…")
    log("Fetching clean 3D structure of native ligand from RCSB/PubChem…", "dim")

    redock_dir = output_dir / "redocking"
    redock_dir.mkdir(parents=True, exist_ok=True)

    fetched_sdf = redock_dir / f"{ligand_code}_online.sdf"
    native_pdbqt: Optional[Path] = None

    if fetch_native_ligand_sdf(
        ligand_code, fetched_sdf,
        ligand_name=ligand_name, pubchem_cid=pubchem_cid, log_cb=log,
    ):
        # Ensure 3D coordinates — PubChem may return 2D-only SDFs
        embedded_sdf = redock_dir / f"{ligand_code}_3d.sdf"
        if ensure_3d_sdf(fetched_sdf, embedded_sdf, log=log):
            fetched_sdf = embedded_sdf

        try:
            native_pdbqt = prepare_ligand_pdbqt(
                sdf_file=fetched_sdf,
                output_dir=redock_dir,
                cmd=config["mk_prepare_ligand_cmd"],
                prefix=f"{ligand_code}_native",
                log_cb=log,
            )
            log(f"✔  Native ligand PDBQT : {Path(native_pdbqt).name}", "success")
        except Exception as e:
            log(f"Failed to prepare online SDF into PDBQT: {e}", "warning")
    else:
        log(
            "Online native ligand 3D retrieval failed — "
            "validation redocking will be skipped.",
            "warning",
        )

    progress_cb(42.0)

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 5 — Input Ligand Preparation  (Meeko)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 5  ─  Input Ligand Preparation  (Meeko)")
    status_cb("Preparing input ligand PDBQT files…")
    log(
        "Each ligand will be checked for 3D coordinates (flat/2D SDFs will be",
        "dim",
    )
    log(
        "embedded via ETKDGv3 + MMFF94) then converted to PDBQT by Meeko.",
        "dim",
    )

    ligand_dir = output_dir / "ligands"
    ligand_dir.mkdir(exist_ok=True)

    prepared_ligands: list[tuple[str, Path, str]] = []
    failed_ligands: list[tuple[str, str]] = []
    n_ligs = len(ligand_files)

    for i, lig_file in enumerate(ligand_files):
        if check_stop():
            return
        orig_stem = Path(lig_file).stem
        clean_name = _get_clean_ligand_name(Path(lig_file))
        log(f"[{i+1}/{n_ligs}]  Preparing : {clean_name}", "info")
        if clean_name != orig_stem:
            log(f"   (source: {Path(lig_file).name})", "dim")
        if ligand_status_cb:
            ligand_status_cb(str(lig_file), "preparing", "")

        try:
            # Ensure 3D coordinates — user-supplied SDFs may be 2D
            lig_path = Path(lig_file)
            embedded_path = ligand_dir / f"{orig_stem}_3d.sdf"
            if ensure_3d_sdf(lig_path, embedded_path, log=log):
                lig_path = embedded_path

            pdbqt = prepare_ligand_pdbqt(
                sdf_file=lig_path,
                output_dir=ligand_dir,
                cmd=config["mk_prepare_ligand_cmd"],
                prefix=orig_stem,
                log_cb=log,
            )
            prepared_ligands.append((clean_name, pdbqt, str(lig_file)))
        except Exception as e:
            log(f"  ✘  Preparation failed for {clean_name}: {e}", "error")
            failed_ligands.append((clean_name, str(e)))
            if ligand_status_cb:
                ligand_status_cb(str(lig_file), "prep_failed", str(e))

        progress_cb(42.0 + (i + 1) * (13.0 / n_ligs))

    log(f"✔  {len(prepared_ligands)} ligand(s) prepared.", "success")
    if failed_ligands:
        log(f"⚠  {len(failed_ligands)} ligand(s) failed preparation:", "warning")
        for fail_name, err in failed_ligands:
            log(f"    • {fail_name}: {err[:120]}", "warning")

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 6 — Validation Redocking
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 6  ─  Validation Redocking")
    status_cb("Running validation redocking…")
    log("Re-docking the native ligand to verify the grid box placement.", "dim")
    log("A successful redock (RMSD ≤ 2.0 Å) confirms reliable production results.", "dim")

    redock_score: float = 0.0
    rmsd: Optional[float] = None
    redock_result: Optional[dict[str, Any]] = None

    if not (native_pdbqt and native_pdbqt.exists()):
        raise RuntimeError("VALIDATION_BLOCKED: native ligand preparation failed; production docking is not permitted.")
    if native_pdbqt and native_pdbqt.exists():
        log(f"Re-docking native ligand [{ligand_code}] ({ligand_name}) back into pocket…", "info")

        redock_result = run_vina_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=native_pdbqt,
            grid=grid,
            output_dir=redock_dir,
            prefix="native_redock",
            vina_exe=config["vina_exe"],
            exhaustiveness=exhaustiveness,
            num_modes=int(config.get("num_modes", 9)),
            cpu=int(config.get("cpu", 0)),
            log_cb=log,
        )

        redock_score = redock_result["best_score"]
        log(f"   Top redock score : {redock_score:.2f} kcal/mol", "info")

        try:
            _compare_redock_atom_counts(native_coords, redock_result["output_pdbqt"], log)
            rmsd = calculate_redock_rmsd(native_coords, redock_result["output_pdbqt"])
            _report_rmsd(rmsd, log)
            if rmsd > 2.0:
                raise RuntimeError(f"VALIDATION_BLOCKED: redocking RMSD {rmsd:.2f} Å exceeds 2.0 Å threshold.")
        except Exception as e:
            raise RuntimeError(f"VALIDATION_BLOCKED: {e}") from e

    progress_cb(65.0)

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 7 — Production Docking
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 7  ─  Production Docking")
    status_cb("Docking ligands…")
    log(f"Docking {len(prepared_ligands)} ligand(s) with exhaustiveness={exhaustiveness}…", "dim")

    dock_dir = output_dir / "docking"
    dock_dir.mkdir(exist_ok=True)

    docking_results: list[dict[str, Any]] = []
    n = len(prepared_ligands)

    for i, (name, pdbqt, lig_path) in enumerate(prepared_ligands):
        if check_stop():
            break

        log("", "dim")
        log(f"[{i+1}/{n}]  Docking : {name}", "bold")
        if ligand_status_cb:
            ligand_status_cb(lig_path, "docking", "")

        try:
            result = run_vina_docking(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=pdbqt,
                grid=grid,
                output_dir=dock_dir,
                prefix=name,
                vina_exe=config["vina_exe"],
                exhaustiveness=exhaustiveness,
                num_modes=int(config.get("num_modes", 9)),
                cpu=int(config.get("cpu", 0)),
                log_cb=log,
            )

            docking_results.append(
                {
                    "name": name,
                    "source_file": lig_path,
                    "best_score": result["best_score"],
                    "all_scores": result["all_scores"],
                    "output_file": result["output_pdbqt"],
                }
            )

            log(f"   ✔  Best score : {result['best_score']:.2f} kcal/mol", "success")
            if ligand_status_cb:
                ligand_status_cb(lig_path, "ok", f"{result['best_score']:.2f}")

        except Exception as e:
            log(f"   ✘  Docking failed for {name}: {e}", "error")
            if ligand_status_cb:
                ligand_status_cb(lig_path, "failed", str(e))
            docking_results.append(
                {
                    "name": name,
                    "best_score": float("inf"),
                    "all_scores": [],
                    "output_file": None,
                    "error": str(e),
                }
            )
        progress_cb(65.0 + (i + 1) * (30.0 / n))

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 8 — Validation & Production Interaction Analysis (PLIP)
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 8  ─  Validation & Production Interaction Analysis  (PLIP)")
    status_cb("Analyzing protein-ligand interactions…")

    interaction_dir = output_dir / "interactions"
    interaction_dir.mkdir(exist_ok=True)
    successful_dockings = [
        result for result in docking_results if result.get("output_file")
    ]
    interaction_csvs: list[Path] = []
    total_analyses = len(successful_dockings) + (1 if redock_result else 0)
    completed_analyses = 0

    if redock_result:
        validation_name = f"Validation {ligand_code} ({ligand_name})"
        log("", "dim")
        log(f"[1/{total_analyses}]  PLIP analysis : {validation_name}", "bold")
        try:
            validation_analysis = run_plip_interaction_analysis(
                receptor_pdbqt=receptor_pdbqt,
                docked_poses_pdbqt=redock_result["output_pdbqt"],
                ligand_name=validation_name,
                output_dir=interaction_dir,
                pymol_exe=config["pymol_exe"],
                plip_cmd=config["plip_cmd"],
                docking_type="validation",
                best_docking_score=redock_result["best_score"],
                all_docking_scores=redock_result["all_scores"],
                log_cb=log,
            )
            redock_result["interaction_analysis"] = validation_analysis
            interaction_csvs.append(validation_analysis["csv_report"])
            counts = validation_analysis["interaction_counts"]
            count_text = ", ".join(
                f"{interaction_type}: {count}"
                for interaction_type, count in counts.items()
            ) or "no interactions detected"
            log(
                f"   ✔  {validation_analysis['interaction_count']} interaction(s): "
                f"{count_text}",
                "success",
            )
        except Exception as exc:
            redock_result["interaction_error"] = str(exc)
            log(f"   ✘  PLIP analysis failed for {validation_name}: {exc}", "error")
        completed_analyses += 1
        progress_cb(95.0 + completed_analyses * (4.0 / max(total_analyses, 1)))

    for index, result in enumerate(successful_dockings):
        if check_stop():
            return
        name = result["name"]
        log("", "dim")
        log(
            f"[{completed_analyses + 1}/{total_analyses}]  PLIP analysis : {name}",
            "bold",
        )
        try:
            analysis = run_plip_interaction_analysis(
                receptor_pdbqt=receptor_pdbqt,
                docked_poses_pdbqt=result["output_file"],
                ligand_name=name,
                output_dir=interaction_dir,
                pymol_exe=config["pymol_exe"],
                plip_cmd=config["plip_cmd"],
                docking_type="production",
                best_docking_score=result["best_score"],
                all_docking_scores=result["all_scores"],
                log_cb=log,
            )
            result["interaction_analysis"] = analysis
            interaction_csvs.append(analysis["csv_report"])
            counts = analysis["interaction_counts"]
            count_text = ", ".join(
                f"{interaction_type}: {count}"
                for interaction_type, count in counts.items()
            ) or "no interactions detected"
            log(
                f"   ✔  {analysis['interaction_count']} interaction(s): {count_text}",
                "success",
            )
        except Exception as exc:
            result["interaction_error"] = str(exc)
            log(f"   ✘  PLIP analysis failed for {name}: {exc}", "error")

        completed_analyses += 1
        progress_cb(95.0 + completed_analyses * (4.0 / max(total_analyses, 1)))

    combined_interaction_csv: Path | None = None
    if interaction_csvs:
        combined_interaction_csv = combine_interaction_csvs(
            interaction_csvs,
            interaction_dir / "interaction_summary.csv",
        )
        log(
            f"✔  Combined interaction table : {combined_interaction_csv.name}",
            "success",
        )

    if check_stop():
        return

    # ════════════════════════════════════════════════════════════════════
    #  STAGE 9 — Summary Report
    # ════════════════════════════════════════════════════════════════════
    banner("STAGE 9  ─  Summary Report")
    status_cb("Generating summary report…")

    report_path = generate_summary_report(
        output_dir=output_dir,
        receptor_file=receptor_file,
        ligand_files=ligand_files,
        active_site_code=ligand_code,
        target_chain=target_chain,
        selection_reason=selection_reason,
        grid=grid,
        redock_rmsd=rmsd,
        redock_score=redock_score,
        docking_results=docking_results,
        grid_padding=grid_pad,
        interaction_csv=combined_interaction_csv,
    )

    log(f"✔  Report written : {report_path.name}", "success")

    llm_context = {
        "experiment": {
            "id": output_dir.name,
            "created": _experiment_created,
            "engine": "AutoDock Vina",
            "exhaustiveness": exhaustiveness,
            "num_modes": int(config.get("num_modes", 9)),
            "energy_range": config.get("energy_range"),
            "seed": config.get("seed"),
        },
        "receptor": {
            "id": Path(receptor_file).stem,
            "chains": [target_chain],
            "prepared_from": Path(receptor_file).name,
            "prep": {"waters_removed": True, "hydrogens_added": True},
        },
        "site": {
            "id": 1,
            "method": "reference_ligand",
            "center": [grid["center_x"], grid["center_y"], grid["center_z"]],
            "size": [grid["size_x"], grid["size_y"], grid["size_z"]],
        },
        "reference": {
            "id": ligand_code,
            "name": ligand_name,
            "best_score": redock_result.get("best_score") if redock_result else None,
            "redock_rmsd": rmsd,
            "validation_threshold_A": 2.0,
            "interaction_csv": str(redock_result.get("interaction_analysis", {}).get("csv_report", "")) if redock_result else "",
        },
        "ligands": [
            {
                "id": result["name"],
                "scores": result.get("all_scores", []),
                "source_file": result.get("source_file"),
                "interaction_csv": str(result.get("interaction_analysis", {}).get("csv_report", "")),
            }
            for result in docking_results if result.get("all_scores")
        ],
        "provenance": {
            "receptor_file": str(receptor_pdbqt),
            "interaction_source": str(combined_interaction_csv) if combined_interaction_csv else None,
            "software_version": "2.0.0",
        },
    }
    manifest_path = output_dir / "dockllm-source.json"
    manifest_path.write_text(json.dumps(llm_context, indent=2, default=str) + "\n", encoding="utf-8")
    llm_path = export_experiment(
        llm_context,
        output_dir / "docking_experiment.llm.json",
        ExportOptions(detail=str(config.get("llm_detail", "standard"))),
    )
    log(f"✔  LLM analysis file : {llm_path.name}", "success")

    # Final ranked table output
    log("", "dim")
    log("─" * 58, "header")
    log("  RANKED DOCKING RESULTS", "header")
    log("─" * 58, "header")

    sorted_r = sorted(docking_results, key=lambda x: x.get("best_score", float("inf")))
    successful = [r for r in sorted_r if r.get("best_score", float("inf")) != float("inf")]
    failed = [r for r in sorted_r if r.get("best_score", float("inf")) == float("inf")]
    col_w = max((len(r["name"]) for r in sorted_r), default=10) + 2

    log(f"  {'Rank':<5}  {'Ligand':<{col_w}}  {'Best (kcal/mol)':>15}", "dim")
    log(f"  {'─'*5}  {'─'*col_w}  {'─'*15}", "dim")

    for rank, r in enumerate(successful, 1):
        log(
            f"  {rank:<5}  {r['name']:<{col_w}}  {r['best_score']:>15.2f}",
            "success",
        )

    if failed:
        log(f"  {'—':<5}  FAILED LIGANDS", "warning")
        for r in failed:
            err_short = r.get("error", "unknown error")[:80]
            log(f"  {'  ':<5}  {r['name']:<{col_w}}  {err_short}", "error")

    # Elapsed time summary
    elapsed = time.time() - _pipeline_start
    minutes, seconds = divmod(int(elapsed), 60)
    log("", "dim")
    log(
        f"Total elapsed time : {minutes}m {seconds}s",
        "dim",
    )
    log("", "dim")
    log("✔  All output saved to :", "success")
    log(f"   {output_dir}", "bold")

    # Write full UI log to MDP-log.txt
    log_path = output_dir / "MDP-log.txt"
    if write_log_fn:
        write_log_fn(log_path)
        log(f"✔  UI log saved : MDP-log.txt", "success")

    progress_cb(100.0)
    status_cb("Pipeline complete ✓")


# ════════════════════════════════════════════════════════════════════════════
#  Internal Helper Functions
# ════════════════════════════════════════════════════════════════════════════

def _fetch_rcsb_chem_metadata(ligand_code: str, retries: int = 3) -> Tuple[str, str]:
    """Query RCSB Data REST API for official chemical name and PubChem CID.

    Retries up to *retries* times with exponential back-off (1 s, 2 s, 4 s).
    If RCSB does not expose a PubChem CID (some CCD entries only carry a
    DrugBank / ChEBI id), the official chemical name is resolved against
    PubChem's cids/JSON endpoint as a backfill.
    """
    ligand_name = ligand_code
    pubchem_cid = "N/A"
    meta_url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{urllib.parse.quote(ligand_code.upper())}"

    for attempt in range(retries):
        try:
            req = urllib.request.Request(meta_url, headers=_HTTP_HEADERS)
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                chem_comp = data.get("chem_comp", {})
                ligand_name = chem_comp.get("name", ligand_code)

                identifiers = data.get("rcsb_chem_comp_identifiers", {}).get("identifiers", [])
                for ident in identifiers:
                    if ident.get("program") == "PubChem":
                        pubchem_cid = str(ident.get("identifier"))
                        break
                break
        except Exception:
            if attempt < retries - 1:
                time.sleep(1 * (2 ** attempt))

    # Backfill: RCSB often has no PubChem identifier for a CCD entry.
    if pubchem_cid in ("N/A", "UNKNOWN", "") and ligand_name and ligand_name != ligand_code:
        resolved = _pubchem_name_to_cid(ligand_name)
        if resolved:
            pubchem_cid = resolved

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
        f"# Generated by Molecular Docking Pipeline — {datetime.now():%Y-%m-%d %H:%M:%S}",
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


def _compare_redock_atom_counts(
    native_coords: List[Tuple[float, float, float]],
    redocked_pdbqt: Path,
    log: Callable[[str, str], None],
) -> None:
    """
    Warn when crystal and online SDF ligands have different heavy-atom
    counts — a common source of 999.0 RMSD returns.

    ``native_coords`` are the crystal heavy-atom coordinates from
    extract_native_ligand_coords(); the online redock PDBQT is parsed for
    its own heavy-atom count.
    """
    from docking import _parse_pdbqt_model1_heavy_coords

    native_count = len(native_coords)
    docked_count = len(_parse_pdbqt_model1_heavy_coords(redocked_pdbqt))

    if native_count != docked_count and native_count > 0 and docked_count > 0:
        log(
            f"⚠  Atom count mismatch: crystal native ligand has {native_count} "
            f"heavy atoms, online SDF redock has {docked_count}.\n"
            f"   This may indicate different protonation / tautomer states "
            f"between the crystal and online structure.\n"
            f"   RMSD may be unreliable — check the summary report.",
            "warning",
        )


def _adaptive_exhaustiveness(
    base: int,
    grid: Dict[str, float],
    cap: int = 64,
) -> Tuple[int, float, str]:
    """
    Scale Vina exhaustiveness based on grid box volume.

    Larger search spaces need more sampling to achieve the same coverage.

    Returns
    -------
    (exhaustiveness, volume, label)
        exhaustiveness : base × multiplier, capped at *cap*
        volume         : grid volume in Å³
        label          : human-readable pocket descriptor

    Volume (Å³)       Descriptor     Multiplier
    ─────────────     ───────────     ──────────
    < 10 000          Small/tight     1.0×
    10 000 – 20 000   Standard        1.5×
    20 000 – 35 000   Large           2.0×
    > 35 000          Very large      3.0×
    """
    volume = grid["size_x"] * grid["size_y"] * grid["size_z"]

    if volume > 35_000:
        multiplier, label = 3.0, "Very large"
    elif volume > 20_000:
        multiplier, label = 2.0, "Large"
    elif volume > 10_000:
        multiplier, label = 1.5, "Standard"
    else:
        multiplier, label = 1.0, "Small/tight"

    scaled = int(base * multiplier)
    adjusted = min(scaled, cap)

    return adjusted, volume, label


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
            "   Verify the selected ligand code is correct.\n"
            "   You may need to manually inspect the receptor and rerun.",
            "error",
        )


# ════════════════════════════════════════════════════════════════════════════
#  Ligand Naming Helpers
# ════════════════════════════════════════════════════════════════════════════

_CID_PATTERN = re.compile(r"[_-]CID[_-](\d+)", re.IGNORECASE)


def _extract_cid_from_filename(filename: str) -> str | None:
    """Extract a PubChem CID from filenames like 'Conformer3D_COMPOUND_CID_2244'."""
    m = _CID_PATTERN.search(filename)
    return m.group(1) if m else None


def _extract_sdf_name(sdf_path: Path) -> str | None:
    """Extract the molecule name from the first MOL block in an SDF file.

    Returns the name from the line after the first '$$$$' delimiter
    (or from the MOL block's first line), or None if unavailable.
    """
    try:
        with open(sdf_path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()

        # Method 1: look for the molecule name line after the first $$$$ delimiter
        # SDF format: <name> is on the line immediately after $$$$ (or at file start)
        name_line = None
        for i, line in enumerate(lines):
            if line.strip() == "$$$$":
                if i + 1 < len(lines):
                    candidate = lines[i + 1].strip()
                    if candidate and candidate != "$$$$" and not candidate.startswith(">"):
                        name_line = candidate
                break

        # Method 2: if file starts with a MOL block, first line is the name
        if not name_line and lines:
            first = lines[0].strip()
            if first and not first.startswith("$$$$") and not first.startswith(">"):
                name_line = first

        if name_line:
            # Clean up the name — remove common prefixes and file extensions
            name = name_line.split(".")[0].strip()
            # Remove "Conformer3D_COMPOUND_" type prefixes
            name = re.sub(r"^Conformer3D[_-]COMPOUND[_-]", "", name, flags=re.IGNORECASE)
            # Remove CID prefix if present (we'll use the extracted CID separately)
            name = re.sub(r"^CID[_-]", "", name, flags=re.IGNORECASE)
            if name:
                return name
    except Exception:
        pass
    return None


def _get_clean_ligand_name(lig_file: Path) -> str:
    """Derive a clean, human-readable name for a ligand from its SDF metadata.

    Priority:
      1. PubChem CID extracted from filename (e.g. CID_2244 → 'CID 2244')
      2. Molecule name from SDF metadata
      3. Original file stem (fallback)
    """
    stem = lig_file.stem

    # Try CID extraction from filename first
    cid = _extract_cid_from_filename(stem)
    if cid:
        return f"CID {cid}"

    # Try reading the SDF molecule name
    sdf_name = _extract_sdf_name(lig_file)
    if sdf_name:
        return sdf_name

    return stem
