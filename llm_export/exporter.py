"""Build compact, self-describing docking experiment JSON files."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA_VERSION = "dockllm-1.0"
RESIDUE_RE = re.compile(r"^[A-Za-z0-9]+:[A-Z0-9]{1,4}-?\d+[A-Za-z]?$")

INTERACTION_CODES = {
    "hydrogen bond": "HB", "hydrophobic contact": "HP",
    "salt bridge": "SB", "pi-stacking": "PI",
    "pi-cation interaction": "PC", "halogen bond": "XB",
    "water bridge": "WB", "metal complex": "MT",
    "steric clash": "CL",
}
CONTACT_FIELDS = {
    "HB": ["type", "residue", "DA_A", "angle_deg", "quality"],
    "HP": ["type", "residue", "distance_A"],
    "SB": ["type", "residue", "distance_A"],
    "PI": ["type", "residue", "centroid_A", "angle_deg"],
    "PC": ["type", "residue", "distance_A"],
    "XB": ["type", "residue", "distance_A", "angle_deg"],
    "WB": ["type", "residue", "distance_A", "angle_deg"],
    "MT": ["type", "residue", "distance_A"],
    "CL": ["type", "residue", "distance_A", "severity"],
}


@dataclass(frozen=True)
class ExportOptions:
    detail: str = "standard"
    pretty: bool = False
    include_coordinates: bool = False
    include_all_pose_interactions: bool = False
    hb_preferred_max_da_A: float = 3.5
    hb_preferred_min_angle_deg: float = 120.0
    score_mismatch_tolerance: float = 0.15


def normalize_residue(chain: Any, name: Any, number: Any) -> str | None:
    """Return CHAIN:RESNAME123, or None when source identity is incomplete."""
    c, n, i = str(chain or "").strip(), str(name or "").strip().upper(), str(number or "").strip()
    if not c or not n or not i:
        return None
    value = f"{c}:{n}{i}"
    return value if RESIDUE_RE.fullmatch(value) else None


def score_stats(scores: Sequence[float]) -> dict[str, float]:
    values = [float(v) for v in scores if math.isfinite(float(v))]
    if not values:
        return {}
    return {"best": _r(min(values), 2), "worst": _r(max(values), 2),
            "range": _r(max(values) - min(values), 2), "mean": _r(sum(values) / len(values), 2)}


def reference_metrics(candidate: Iterable[Sequence[Any]], reference: Iterable[Sequence[Any]]) -> dict[str, Any]:
    candidate_pairs = {(str(x[0]), str(x[1])) for x in candidate if len(x) >= 2}
    reference_pairs = {(str(x[0]), str(x[1])) for x in reference if len(x) >= 2}
    candidate_residues = {x[1] for x in candidate_pairs}
    reference_residues = {x[1] for x in reference_pairs}
    shared = candidate_residues & reference_residues
    union = candidate_residues | reference_residues
    exact_contacts = [[kind, residue] for kind, residue in sorted(candidate_pairs & reference_pairs, key=lambda x: (x[1], x[0]))]
    return {
        "shared": len(shared),
        "candidate_residue_count": len(candidate_residues),
        "reference_residue_count": len(reference_residues),
        "candidate_overlap": _ratio(len(shared), len(candidate_residues)),
        "reference_coverage": _ratio(len(shared), len(reference_residues)),
        "jaccard": _ratio(len(shared), len(union)),
        "exact": len(exact_contacts),
        "exact_contacts": exact_contacts,
    }


def parse_interactions(path: str | Path | None, options: ExportOptions) -> list[list[Any]]:
    if not path or not Path(path).is_file():
        return []
    contacts: list[list[Any]] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            code = INTERACTION_CODES.get((row.get("interaction_type") or "").strip().lower())
            residue = normalize_residue(row.get("reschain"), row.get("restype"), row.get("resnr"))
            if not code or not residue:
                continue
            contact: list[Any] = [code, residue]
            if code == "HB":
                distance = _first_number(row, "dist_d-a", "dist_h-a", "dist")
                angle = _first_number(row, "don_angle", "angle")
                if distance is not None:
                    _append(contact, distance, 2)
                    if angle is not None:
                        _append(contact, angle, 0)
                        contact.append(classify_hbond(distance, angle, options))
            elif code == "PI":
                distance = _first_number(row, "cent_dist", "dist")
                if distance is not None:
                    _append(contact, distance, 2)
                    _append(contact, _first_number(row, "angle"), 0)
            else:
                _append(contact, _first_number(row, "dist", "dist_a-w", "cent_dist"), 2)
                if code in {"XB", "WB"}:
                    _append(contact, _first_number(row, "angle", "don_angle", "water_angle"), 0)
            if options.include_coordinates:
                coords = _coordinates(row)
                if coords:
                    contact.append(coords)
            contacts.append(contact)
    return sorted(contacts, key=lambda x: (str(x[1]), str(x[0]), json.dumps(x, separators=(",", ":"))))


def classify_hbond(distance: float, angle: float, options: ExportOptions) -> str:
    if distance <= 2.7 and angle >= 150:
        return "S"
    if distance <= 3.0 and angle >= 130:
        return "G"
    if distance <= options.hb_preferred_max_da_A and angle >= options.hb_preferred_min_angle_deg:
        return "W"
    return "B"


def ligand_properties(sdf_path: str | Path | None) -> tuple[dict[str, Any], dict[str, str]]:
    if not sdf_path or not Path(sdf_path).is_file():
        return {}, {}
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
        mol = next((m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m is not None), None)
        if mol is None:
            return {}, {}
        props = {"MW": _r(Descriptors.MolWt(mol), 1), "HA": int(mol.GetNumHeavyAtoms()),
                 "HBD": int(Lipinski.NumHDonors(mol)), "HBA": int(Lipinski.NumHAcceptors(mol)),
                 "RB": int(Lipinski.NumRotatableBonds(mol)),
                 "charge": int(sum(a.GetFormalCharge() for a in mol.GetAtoms())),
                 "TPSA": _r(rdMolDescriptors.CalcTPSA(mol), 1), "cLogP": _r(Crippen.MolLogP(mol), 2)}
        identity = {"smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)}
        if mol.HasProp("PUBCHEM_IUPAC_NAME"):
            identity["name"] = mol.GetProp("PUBCHEM_IUPAC_NAME")
        if mol.HasProp("PUBCHEM_IUPAC_INCHIKEY"):
            identity["inchikey"] = mol.GetProp("PUBCHEM_IUPAC_INCHIKEY")
        return props, identity
    except (ImportError, OSError, ValueError):
        return {}, {}


def build_export(context: Mapping[str, Any], options: ExportOptions | None = None) -> dict[str, Any]:
    options = options or ExportOptions()
    if options.detail not in {"full", "standard", "compact"}:
        raise ValueError("detail must be full, standard, or compact")
    reference_source = context.get("reference") or {}
    reference_contacts = parse_interactions(reference_source.get("interaction_csv"), options)
    reference_residues = sorted({c[1] for c in reference_contacts})
    ligands: list[dict[str, Any]] = []
    requested_modes = _integer((context.get("experiment") or {}).get("num_modes"))
    for source in context.get("ligands", []):
        scores = [_r(float(v), 2) for v in source.get("scores", []) if _finite(v)]
        if not scores:
            continue
        contacts = parse_interactions(source.get("interaction_csv"), options)
        counts = dict(sorted(Counter(c[0] for c in contacts).items()))
        props, identity = ligand_properties(source.get("source_file"))
        ligand: dict[str, Any] = {"id": str(source["id"]), **identity,
            "scores": scores, "best": min(scores), "score_stats": score_stats(scores)}
        if source.get("name") and source.get("name") != source["id"]:
            ligand["name"] = source["name"]
        if props:
            ligand["props"] = props
            if props.get("HA"):
                ligand["score_per_HA"] = _r(min(scores) / props["HA"], 3)
        ligand["counts"] = counts
        residues = sorted({c[1] for c in contacts})
        ligand["residues"] = residues
        ligand["contacts"] = contacts
        if options.detail != "compact":
            pass
        elif contacts:
            ligand["key_contacts"] = [c for c in contacts if c[0] in {"HB", "SB", "PI", "PC", "XB", "MT", "CL"}]
        ligand["ref"] = reference_metrics(contacts, reference_contacts)
        quality_counts = Counter(c[4] for c in contacts if c[0] == "HB" and len(c) > 4)
        quality: dict[str, Any] = {"strong_HB": quality_counts["S"], "good_HB": quality_counts["G"],
                                  "weak_HB": quality_counts["W"], "borderline_HB": quality_counts["B"]}
        clashes = [c for c in contacts if c[0] == "CL"]
        if clashes:
            quality["clashes"] = len(clashes)
        ligand["quality"] = quality
        ligand["flags"] = ligand_flags(ligand, requested_modes, source, options)
        if options.include_all_pose_interactions:
            pose_contacts = []
            for pose_number, interaction_path in enumerate(source.get("pose_interaction_csvs", []), 1):
                pose_contacts.append({"pose": pose_number, "contacts": parse_interactions(interaction_path, options)})
            if pose_contacts:
                ligand["pose_interactions"] = pose_contacts
        if options.detail == "compact":
            ligand.pop("contacts", None)
        ligands.append(ligand)
    ligands.sort(key=lambda x: (x["best"], x["id"]))
    for rank, ligand in enumerate(ligands, 1):
        ligand["rank"] = rank
    ranking = [[x["id"], x["best"], x["rank"]] for x in ligands]
    frequent = Counter(r for ligand in ligands for r in ligand.get("residues", []))
    global_warnings = [[x["id"], flag] for x in ligands for flag in x["flags"]]
    result: dict[str, Any] = {"schema": SCHEMA_VERSION,
        "experiment": _experiment(context.get("experiment") or {}, options),
        "receptor": _clean(context.get("receptor") or {}),
        "site": _site(context.get("site") or {}, reference_residues),
        "reference": _reference(reference_source, reference_contacts),
        "ranking": ranking,
        "global_analysis": {"n_screened": len(ligands),
            "best_candidate": ligands[0]["id"] if ligands else None,
            "reference_score": _optional_round(reference_source.get("best_score"), 2),
            "best_candidate_score": ligands[0]["best"] if ligands else None,
            "frequent_residues": [[r, n] for r, n in sorted(frequent.items(), key=lambda x: (-x[1], x[0]))],
            "reference_like": [x["id"] for x in sorted(ligands, key=lambda x: (-x["ref"]["jaccard"], x["id"])) if x["ref"]["shared"] > 0],
            "warnings": global_warnings},
        "ligands": ligands, "codebook": _codebook(),
        "provenance": _provenance(context.get("provenance") or {})}
    return result


def export_experiment(context: Mapping[str, Any], output: str | Path, options: ExportOptions | None = None) -> Path:
    options = options or ExportOptions()
    payload = build_export(context, options)
    validate_export(payload)
    separators = None if options.pretty else (",", ":")
    text = json.dumps(payload, ensure_ascii=False, indent=2 if options.pretty else None,
                      separators=separators, allow_nan=False) + "\n"
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)
    return path


def validate_export(data: Mapping[str, Any]) -> None:
    required = ("schema", "experiment", "receptor", "site", "reference", "ranking", "global_analysis", "ligands", "codebook", "provenance")
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"missing required fields: {', '.join(missing)}")
    if data["schema"] != SCHEMA_VERSION:
        raise ValueError("unsupported DockLLM schema")
    ids = [x.get("id") for x in data["ligands"]]
    if any(not x for x in ids) or len(ids) != len(set(ids)):
        raise ValueError("ligand IDs must be present and unique")
    if any(not isinstance(score, (int, float)) for x in data["ligands"] for score in x.get("scores", [])):
        raise ValueError("docking scores must be numeric")
    if any(not RESIDUE_RE.fullmatch(c[1]) for x in data["ligands"] for c in x.get("contacts", [])):
        raise ValueError("interaction residue is not normalized")


def ligand_flags(ligand: Mapping[str, Any], requested_modes: int | None, source: Mapping[str, Any], options: ExportOptions) -> list[str]:
    flags: list[str] = []
    if ligand["best"] > 0:
        flags.append("POSITIVE_SCORE")
    if requested_modes and len(ligand["scores"]) < requested_modes:
        flags.append("FEW_POSES")
    if not ligand.get("residues"):
        flags.append("NO_INTERACTIONS")
    if ligand.get("ref", {}).get("candidate_overlap", 0) < 0.25 and ligand.get("residues"):
        flags.append("LOW_REFERENCE_OVERLAP")
    hb = [c for c in ligand.get("contacts", []) if c[0] == "HB"]
    if hb and all(len(c) > 4 and c[4] in {"W", "B"} for c in hb):
        flags.append("HB_ALL_WEAK_OR_BORDERLINE")
    sources = {k: float(v) for k, v in (source.get("score_sources") or {}).items() if _finite(v)}
    if sources and max(sources.values()) - min(sources.values()) > options.score_mismatch_tolerance:
        flags.append("SCORE_SOURCE_MISMATCH")
        ligand["score_sources"] = {k: _r(v, 2) for k, v in sorted(sources.items())}
    return flags


def _experiment(source: Mapping[str, Any], options: ExportOptions) -> dict[str, Any]:
    result = {"id": str(source.get("id") or "docking_experiment"), "engine": source.get("engine") or "AutoDock Vina",
              "score_unit": "kcal/mol", "lower_score_is_better": True, "export_detail": options.detail,
              "all_pose_interactions": options.include_all_pose_interactions,
              "coordinates_included": options.include_coordinates,
              "HB_quality_rules": {
                  "preferred_max_DA_A": options.hb_preferred_max_da_A,
                  "preferred_min_angle_deg": options.hb_preferred_min_angle_deg,
                  "outside_preferred_geometry": "detector-reported hydrogen bonds are retained and classified borderline",
              }}
    for key in ("created", "engine_version", "exhaustiveness", "num_modes", "energy_range", "seed"):
        if source.get(key) is not None:
            result[key] = source[key]
    return result


def _site(source: Mapping[str, Any], reference_residues: list[str]) -> dict[str, Any]:
    result = _clean(source)
    if "center" in result:
        result["center"] = [_r(float(x), 2) for x in result["center"]]
    if "size" in result:
        result["size"] = [_r(float(x), 2) for x in result["size"]]
    if reference_residues:
        result["reference_residues"] = reference_residues
    return result


def _reference(source: Mapping[str, Any], contacts: list[list[Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {"id": str(source.get("id") or "reference"), "role": "validation"}
    for key in ("name", "best_score"):
        if source.get(key) is not None:
            result[key] = _optional_round(source[key], 2) if key == "best_score" else source[key]
    rmsd = source.get("redock_rmsd")
    result["redock_rmsd"] = _optional_round(rmsd, 2)
    result["validation_threshold_A"] = _r(float(source.get("validation_threshold_A", 2.0)), 2)
    result["validation_pass"] = None if rmsd is None else float(rmsd) <= result["validation_threshold_A"]
    result["contacts"] = contacts
    result["contact_residues"] = sorted({c[1] for c in contacts})
    if rmsd is None:
        result["flags"] = ["MISSING_VALIDATION_RMSD"]
    return result


def _codebook() -> dict[str, Any]:
    return {"interaction": {"HB": "hydrogen_bond", "HP": "hydrophobic", "SB": "salt_bridge",
                "PI": "pi_stacking", "PC": "pi_cation", "XB": "halogen_bond", "WB": "water_bridge",
                "MT": "metal_coordination", "CL": "steric_clash"},
            "contact_fields": CONTACT_FIELDS, "HB_quality": {"S": "strong", "G": "good", "W": "weak", "B": "borderline"},
            "ranking_fields": ["ligand", "best_score", "rank"],
            "reference_match_fields": {"exact_contacts": ["interaction_type", "residue"]},
            "abbreviations": {"score_per_HA": "best docking score divided by heavy atom count; more negative is better when lower_score_is_better=true"},
            "flags": {"POSITIVE_SCORE": "Best docking score is positive.", "FEW_POSES": "Docking produced fewer poses than requested.",
                "NO_INTERACTIONS": "No best-pose interactions were available.", "SEVERE_CLASH": "One or more severe steric clashes were detected.",
                "LOW_REFERENCE_OVERLAP": "Candidate residue overlap with reference is below 0.25.", "HB_ALL_WEAK_OR_BORDERLINE": "All detected hydrogen bonds are weak or borderline; other interaction types may also exist.",
                "MISSING_VALIDATION_RMSD": "Reference redocking RMSD is unavailable.", "SCORE_SOURCE_MISMATCH": "Score sources differ beyond configured tolerance."}}


def _provenance(source: Mapping[str, Any]) -> dict[str, Any]:
    result = _clean(source)
    hashes: dict[str, str] = {}
    for key, value in source.items():
        if key.endswith("_file") and value and Path(str(value)).is_file():
            hashes[key] = hashlib.sha256(Path(str(value)).read_bytes()).hexdigest()
    if hashes:
        result["sha256"] = dict(sorted(hashes.items()))
    return result


def _coordinates(row: Mapping[str, str]) -> list[list[float]]:
    values: list[list[float]] = []
    for prefix in ("ligcoo", "protcoo"):
        xyz = [_number(row.get(f"{prefix}.{axis}")) for axis in "xyz"]
        if all(v is not None for v in xyz):
            values.append([_r(float(v), 2) for v in xyz])
    return values


def _clean(source: Mapping[str, Any]) -> dict[str, Any]:
    return {str(k): v for k, v in source.items() if v is not None and v != ""}


def _first_number(row: Mapping[str, str], *keys: str) -> float | None:
    for key in keys:
        value = _number(row.get(key))
        if value is not None:
            return value
    return None


def _number(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _finite(value: Any) -> bool:
    return _number(value) is not None


def _append(target: list[Any], value: float | None, digits: int) -> None:
    if value is not None:
        target.append(_r(value, digits))


def _ratio(numerator: int, denominator: int) -> float:
    return _r(numerator / denominator, 3) if denominator else 0.0


def _r(value: float, digits: int) -> float:
    return round(float(value), digits)


def _optional_round(value: Any, digits: int) -> float | None:
    return None if not _finite(value) else _r(float(value), digits)


def _integer(value: Any) -> int | None:
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None
