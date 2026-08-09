"""PLIP analysis for validation and production docking poses.

Builds a protein-ligand PDB complex from receptor/pose PDBQT files, runs
PLIP, and converts its complete XML interaction report into analysis-ready CSV.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import sys
import textwrap
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Callable


_INTERACTION_NAMES = {
    "hydrophobic_interaction": "Hydrophobic contact",
    "hydrogen_bond": "Hydrogen bond",
    "water_bridge": "Water bridge",
    "salt_bridge": "Salt bridge",
    "pi_stack": "Pi-stacking",
    "pi_cation_interaction": "Pi-cation interaction",
    "halogen_bond": "Halogen bond",
    "metal_complex": "Metal complex",
}

_PREFERRED_COLUMNS = [
    "ligand_name",
    "docking_type",
    "best_docking_score_kcal_mol",
    "all_docking_scores_kcal_mol",
    "binding_site_id",
    "interaction_type",
    "ligand_hetid",
    "ligand_chain",
    "ligand_position",
    "restype",
    "resnr",
    "reschain",
    "restype_lig",
    "resnr_lig",
    "reschain_lig",
    "dist",
    "dist_h-a",
    "dist_d-a",
    "dist_a-w",
    "dist_d-w",
    "cent_dist",
    "don_angle",
    "acc_angle",
    "water_angle",
    "angle",
    "offset",
    "donoridx",
    "donortype",
    "acceptoridx",
    "acceptortype",
    "protcarbonidx",
    "ligcarbonidx",
    "metal_idx",
    "metal_type",
    "target_idx",
    "target_type",
    "details_json",
]


def run_plip_interaction_analysis(
    receptor_pdbqt: Path,
    docked_poses_pdbqt: Path,
    ligand_name: str,
    output_dir: Path,
    pymol_exe: str,
    plip_cmd: str,
    docking_type: str = "production",
    best_docking_score: float | None = None,
    all_docking_scores: list[float] | None = None,
    log_cb: Callable[[str, str], None] | None = None,
) -> dict:
    """Analyze a top-ranked docked pose and return generated report paths."""

    def log(message: str, level: str = "info") -> None:
        if log_cb:
            log_cb(message, level)

    safe_name = _safe_filename(ligand_name)
    analysis_dir = Path(output_dir).resolve() / safe_name
    analysis_dir.mkdir(parents=True, exist_ok=True)
    complex_pdb = analysis_dir / f"{safe_name}_complex.pdb"

    _build_complex_pdb(
        receptor_pdbqt=Path(receptor_pdbqt),
        docked_poses_pdbqt=Path(docked_poses_pdbqt),
        output_pdb=complex_pdb,
        pymol_exe=pymol_exe,
        log_cb=log,
    )

    command = [
        *_resolve_plip_command(plip_cmd),
        "-f",
        str(complex_pdb),
        "-o",
        str(analysis_dir),
        "-x",
        "-t",
        "-q",
    ]
    log(f"   PLIP command : {' '.join(command)}", "dim")
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"PLIP executable not found: {plip_cmd}") from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("PLIP interaction analysis timed out after 300 s.") from exc

    _relay_plip_output((result.stdout or "") + (result.stderr or ""), log)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "No diagnostic output")[-1000:]
        raise RuntimeError(f"PLIP failed with exit code {result.returncode}: {detail}")

    xml_report = _find_plip_report(analysis_dir, complex_pdb.stem, "xml")
    text_report = _find_plip_report(analysis_dir, complex_pdb.stem, "txt")
    if not xml_report.is_file() or xml_report.stat().st_size == 0:
        raise RuntimeError("PLIP completed but an XML interaction report was not created.")
    if not text_report.is_file() or text_report.stat().st_size == 0:
        raise RuntimeError("PLIP completed but a text interaction report was not created.")

    csv_report = analysis_dir / "interactions.csv"
    parsed = parse_plip_xml_to_csv(
        xml_report,
        csv_report,
        ligand_name,
        docking_type=docking_type,
        best_docking_score=best_docking_score,
        all_docking_scores=all_docking_scores,
    )
    return {
        "complex_pdb": complex_pdb,
        "xml_report": xml_report,
        "text_report": text_report,
        "csv_report": csv_report,
        **parsed,
    }


def parse_plip_xml_to_csv(
    xml_report: Path,
    csv_report: Path,
    ligand_name: str,
    docking_type: str = "production",
    best_docking_score: float | None = None,
    all_docking_scores: list[float] | None = None,
) -> dict:
    """Flatten every PLIP XML interaction and all atom/geometry fields to CSV."""
    root = ET.parse(xml_report).getroot()
    rows: list[dict[str, str]] = []

    binding_sites = root.findall(".//bindingsite")
    docked_sites = [site for site in binding_sites if _is_docked_ligand_site(site)]
    selected_sites = docked_sites or binding_sites

    for site in selected_sites:
        identifiers = site.find("identifiers")
        ligand_ids = {
            f"ligand_{child.tag}": (child.text or "").strip()
            for child in list(identifiers or [])
            if len(child) == 0
        }
        interactions = site.find("interactions")
        if interactions is None:
            continue
        for group in list(interactions):
            for item in list(group):
                details: dict[str, str] = {}
                _flatten_xml(item, details)
                interaction_tag = _strip_namespace(item.tag)
                row = {
                    "ligand_name": ligand_name,
                    "docking_type": docking_type,
                    "best_docking_score_kcal_mol": _format_score(
                        best_docking_score
                    ),
                    "all_docking_scores_kcal_mol": ";".join(
                        _format_score(score) for score in (all_docking_scores or [])
                    ),
                    "binding_site_id": site.get("id", ""),
                    "interaction_type": _INTERACTION_NAMES.get(
                        interaction_tag,
                        interaction_tag.replace("_", " ").title(),
                    ),
                    **ligand_ids,
                    **details,
                }
                row["details_json"] = json.dumps(details, sort_keys=True)
                rows.append(row)

    _write_rows(csv_report, rows)
    counts = Counter(row["interaction_type"] for row in rows)
    return {
        "interaction_count": len(rows),
        "interaction_counts": dict(sorted(counts.items())),
    }


def combine_interaction_csvs(csv_paths: list[Path], output_csv: Path) -> Path:
    """Combine per-ligand PLIP tables into one analysis table."""
    rows: list[dict[str, str]] = []
    for path in csv_paths:
        if not Path(path).is_file():
            continue
        with Path(path).open(newline="", encoding="utf-8") as handle:
            rows.extend(dict(row) for row in csv.DictReader(handle))
    _write_rows(output_csv, rows)
    return output_csv


def _build_complex_pdb(
    receptor_pdbqt: Path,
    docked_poses_pdbqt: Path,
    output_pdb: Path,
    pymol_exe: str,
    log_cb: Callable[[str, str], None],
) -> None:
    """Use PyMOL to combine receptor with model 1 of docked ligand poses."""
    receptor = Path(receptor_pdbqt).resolve()
    poses = Path(docked_poses_pdbqt).resolve()
    output = Path(output_pdb).resolve()
    if not receptor.is_file() or not poses.is_file():
        raise RuntimeError("Cannot build PLIP complex: receptor or docked poses missing.")

    script_path = output.parent / "_build_plip_complex.py"
    script = textwrap.dedent(
        f"""\
        import sys
        from pymol import cmd

        RECEPTOR = {json.dumps(str(receptor))}
        POSES = {json.dumps(str(poses))}
        OUTPUT = {json.dumps(str(output))}

        cmd.reinitialize()
        cmd.load(RECEPTOR, "receptor")
        cmd.load(POSES, "all_poses")
        cmd.create("ligand_pose", "all_poses", 1, 1)
        cmd.delete("all_poses")
        cmd.alter("ligand_pose", "resn='LIG'; chain='Z'; resi='1'; segi=''; type='HETATM'")
        cmd.sort("ligand_pose")
        if cmd.count_atoms("receptor") == 0 or cmd.count_atoms("ligand_pose") == 0:
            print("ERROR: receptor or top docked pose contains no atoms")
            cmd.quit()
            sys.exit(2)
        cmd.save(OUTPUT, "receptor or ligand_pose", state=1, format="pdb")
        print(f"Saved PLIP complex: {{OUTPUT}}")
        cmd.quit()
        """
    )
    script_path.write_text(script, encoding="utf-8")
    try:
        result = _run_pymol(pymol_exe, script_path)
        _relay_plip_output((result.stdout or "") + (result.stderr or ""), log_cb, "PyMOL")
        if result.returncode not in (0, 1):
            raise RuntimeError(
                f"PyMOL complex assembly failed with exit code {result.returncode}."
            )
    finally:
        script_path.unlink(missing_ok=True)

    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("PyMOL did not create protein-ligand complex PDB.")
    content = output.read_text(encoding="utf-8", errors="replace")
    if "ATOM  " not in content or "HETATM" not in content:
        raise RuntimeError(
            "PLIP complex PDB must contain receptor ATOM and ligand HETATM records."
        )


def _run_pymol(
    pymol_exe: str,
    script_path: Path,
) -> subprocess.CompletedProcess[str]:
    for flags in (["-cq"], ["-c"], []):
        try:
            result = subprocess.run(
                [pymol_exe, *flags, str(script_path)],
                capture_output=True,
                text=True,
                timeout=180,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(f"PyMOL executable not found: {pymol_exe}") from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("PyMOL complex assembly timed out after 180 s.") from exc
        if "unknown option" not in (result.stderr or "").lower():
            return result
    return result


def _resolve_plip_command(plip_cmd: str) -> list[str]:
    command = str(plip_cmd).strip()
    if not command:
        raise RuntimeError("PLIP command is not configured in Settings.")
    path = Path(command).expanduser()
    if command.lower().endswith(".py") and path.is_file():
        return [sys.executable, "-s", str(path.resolve())]
    if path.is_file():
        return [str(path.resolve())]
    resolved = shutil.which(command)
    if resolved:
        return [resolved]
    raise RuntimeError(f"PLIP executable or plipcmd.py not found: {command}")


def _find_plip_report(analysis_dir: Path, input_stem: str, extension: str) -> Path:
    """Find the input-prefixed report emitted by PLIP, with legacy fallback."""
    candidates = (
        analysis_dir / f"{input_stem}_report.{extension}",
        analysis_dir / f"report.{extension}",
    )
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return candidates[0]


def _format_score(score: float | None) -> str:
    return "" if score is None else f"{score:.3f}"


def _is_docked_ligand_site(site: ET.Element) -> bool:
    identifiers = site.find("identifiers")
    if identifiers is None:
        return False
    hetid = (identifiers.findtext("hetid") or "").strip().upper()
    chain = (identifiers.findtext("chain") or "").strip().upper()
    return hetid == "LIG" and chain == "Z"


def _flatten_xml(
    element: ET.Element,
    output: dict[str, str],
    prefix: str = "",
) -> None:
    for key, value in element.attrib.items():
        _append_value(output, f"{prefix}@{key}" if prefix else f"@{key}", value)
    children = list(element)
    if not children:
        key = prefix or _strip_namespace(element.tag)
        _append_value(output, key, (element.text or "").strip())
        return
    for child in children:
        tag = _strip_namespace(child.tag)
        child_prefix = f"{prefix}.{tag}" if prefix else tag
        _flatten_xml(child, output, child_prefix)


def _append_value(output: dict[str, str], key: str, value: str) -> None:
    if key in output and value:
        output[key] = f"{output[key]}|{value}"
    else:
        output[key] = value


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    discovered = {key for row in rows for key in row}
    fieldnames = [key for key in _PREFERRED_COLUMNS if key in discovered]
    fieldnames.extend(sorted(discovered - set(fieldnames)))
    if not fieldnames:
        fieldnames = _PREFERRED_COLUMNS[:3]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _relay_plip_output(
    text: str,
    log_cb: Callable[[str, str], None],
    source: str = "PLIP",
) -> None:
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        upper = line.upper()
        level = "error" if "ERROR" in upper else "warning" if "WARNING" in upper else "dim"
        log_cb(f"   {source} > {line}", level)


def _safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned or "ligand"


def _strip_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]
