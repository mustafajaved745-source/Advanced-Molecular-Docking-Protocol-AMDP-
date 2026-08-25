"""
Molecular Docking Pipeline
==========================
hetatm_parser.py — HETATM parsing for active-site selection

This module reads the raw receptor structure (PDB or mmCIF) and extracts
candidate ligands for the deterministic, user-controlled "Active Site
Selection" step in the GUI.  It performs no model inference of any kind.

Responsibilities:
  1. Parse all HETATM atoms from a raw receptor structure
  2. Filter out water, ions, and common buffer / cryoprotectant molecules
  3. Group the surviving candidate ligands per chain so the GUI can
     populate its "Active Site Selection" panel

Public API
----------
  iter_hetatm_atoms(structure_file)    -> Iterator[dict]  (per HETATM atom)
  parse_hetatm_records(structure_file) -> List[dict]      (per ligand group)
  build_chain_ligand_map(hetatms)      -> Dict[str, List[dict]]
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Tuple

from mmcif_parser import is_mmcif, iter_atom_site_rows


@dataclass(frozen=True)
class LigandInstanceKey:
    """Full structural identity of one selected ligand conformer."""

    model_num: str
    component_id: str
    auth_chain_id: str
    label_asym_id: str
    auth_seq_id: str
    label_seq_id: str
    insertion_code: str = ""
    altloc: str = ""
    occupancy: float = 1.0

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping) -> "LigandInstanceKey":
        return cls(
            model_num=str(value.get("model_num", "1")),
            component_id=str(value.get("component_id", "")).upper(),
            auth_chain_id=str(value.get("auth_chain_id", "")),
            label_asym_id=str(value.get("label_asym_id", "")),
            auth_seq_id=str(value.get("auth_seq_id", "")),
            label_seq_id=str(value.get("label_seq_id", "")),
            insertion_code=str(value.get("insertion_code", "")),
            altloc=str(value.get("altloc", "")),
            occupancy=float(value.get("occupancy", 1.0)),
        )

    @property
    def chain_id(self) -> str:
        return self.auth_chain_id or self.label_asym_id

    @property
    def sequence_id(self) -> str:
        return self.auth_seq_id or self.label_seq_id

    def display(self) -> str:
        alt = self.altloc or "none"
        ins = self.insertion_code or "none"
        return (
            f"model={self.model_num}; component={self.component_id}; "
            f"auth_chain={self.auth_chain_id or '-'}; label_asym={self.label_asym_id or '-'}; "
            f"auth_seq={self.auth_seq_id or '-'}; label_seq={self.label_seq_id or '-'}; "
            f"ins={ins}; altloc={alt}; occupancy={self.occupancy:.3f}"
        )


@dataclass(frozen=True)
class LigandAtomReference:
    atom_id: str
    element: str
    xyz: Tuple[float, float, float]
    altloc: str
    occupancy: float


@dataclass(frozen=True)
class DockingTargetContext:
    """Immutable selected instance shared by grid, chemistry and RMSD stages."""

    source_file: str
    instance_key: LigandInstanceKey
    crystal_atoms: Tuple[LigandAtomReference, ...]

    @property
    def component_id(self) -> str:
        return self.instance_key.component_id

    @property
    def crystal_coords(self) -> List[Tuple[float, float, float]]:
        return [atom.xyz for atom in self.crystal_atoms]

    @property
    def crystal_atom_ids(self) -> List[str]:
        return [atom.atom_id for atom in self.crystal_atoms]


# ════════════════════════════════════════════════════════════════════════════
#  Reference sets
# ════════════════════════════════════════════════════════════════════════════

# Everything in this set is treated as a definite non-ligand and filtered
# out before anything reaches the selection panel.
NON_LIGAND_RESNS: frozenset[str] = frozenset({
    # ── Water ────────────────────────────────────────────────────────────
    "HOH", "WAT", "DOD", "H2O",

    # ── Monovalent / divalent ions ────────────────────────────────────────
    "NA",  "K",   "LI",  "RB",  "CS",           # alkali metals
    "MG",  "CA",  "SR",  "BA",                   # alkaline earth
    "ZN",  "CU",  "FE",  "MN",  "CO",           # transition metals (common)
    "NI",  "CD",  "HG",  "PB",  "CR",  "MO",
    "CL",  "BR",  "IOD", "I",   "F",            # halogens / anions
    "XE",  "KR",                                 # noble gases (cryo)

    # ── Inorganic polyatomic ions ─────────────────────────────────────────
    "SO4", "SUL", "PO4", "PHO", "NO3", "CO3",
    "BO4", "CLO", "SCN", "NCO",

    # ── Cryoprotectants ───────────────────────────────────────────────────
    "GOL", "EDO", "EG",  "PEG", "PG4", "PE4",
    "PE5", "P6G", "MPD", "MRD", "DMS", "DMSO",

    # ── Common crystallisation additives / buffers ────────────────────────
    "ACT", "ACY", "ACE",                         # acetate / acetyl
    "FMT",                                        # formate
    "TRS",                                        # tris
    "MES", "HEP", "EPE", "PIP",                  # HEPES / PIPES / piperazine
    "IMD",                                        # imidazole
    "DTT", "DTE", "BME", "MSE",                  # reducing agents
    "TAR", "TLA", "SUC", "TBU",                  # tartrate / succinate
    "NH4", "NH2",                                 # ammonium
    "EOH", "MOH", "IPH",                         # small alcohols
    "PGE", "PG",                                  # propylene glycol variants
    "BU1", "BU2", "BTB",                         # butanols / BIS-TRIS
    "MLA", "MLI", "MOL",                         # malonate / malate
    "CIT", "ISO",                                 # citrate / isopropanol

    # ── Detergents & lipids (often artefacts) ────────────────────────────
    "OLA", "PLM", "SDS", "LMT", "C8E",

    # ── Structural lipids (membrane components) ───────────────────────────
    "3PE", "POPC", "DPPC", "DPG", "CDL",
    "OLB", "BOG", "OCT", "D10", "P4G",
    "LMU", "LPP", "MYR", "STE", "PAM",

    # ── Glycans & carbohydrates (non-drug) ────────────────────────────────
    "NAG", "BMA", "MAN", "FUL", "GAL", "BGC", "NDG",
    "GLC", "SIA", "FUC", "RHA", "XYL", "RIB",
    "ARA", "GLA", "GLCN", "GALN", "SOR",
    "G1P", "G6P", "F6P",
    "MAL",  # maltose
    "LAC",  # lactose

    # ── Amino-acid fragments / protecting groups ──────────────────────────
    "NME",
})


# ════════════════════════════════════════════════════════════════════════════
#  HETATM atom iterator (PDB or mmCIF)
# ════════════════════════════════════════════════════════════════════════════

def iter_hetatm_atoms(structure_file: str | Path) -> Iterator[Dict]:
    """
    Yield one normalized dict per HETATM atom in *structure_file*.

    Format is chosen by ``mmcif_parser.is_mmcif``: mmCIF files are read
    through the ``_atom_site`` loop parser, PDB files through fixed-column
    record parsing.  Both branches yield identical dict shapes:

        group, resn, chain, resi, icode, x, y, z,
        element, atom_name, altloc, occ, bfac
    """
    path = Path(structure_file)

    if is_mmcif(path):
        for atom in iter_atom_site_rows(path):
            if atom["group"] == "HETATM":
                yield atom
        return

    # ── PDB fixed-column parsing ────────────────────────────────────────
    model_num = "1"
    in_first_model = True
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                model_num = line[10:14].strip() or line[5:].strip() or "1"
                in_first_model = model_num == "1"
                continue
            if line.startswith("ENDMDL"):
                if in_first_model:
                    break
                continue
            if not in_first_model:
                continue
            if not line.startswith("HETATM"):
                continue

            try:
                resn  = line[17:20].strip().upper()
                chain = line[21].strip() or "A"
                resi  = line[22:26].strip()
                icode = line[26].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                element = ""
                if len(line) >= 78:
                    element = line[76:78].strip().upper()

                atom_name = line[12:16].strip()
                if not element:
                    # Fallback parsing for element from atom name column (12-16)
                    clean_name = re.sub(r"^\d+", "", atom_name)
                    element = clean_name[:1].upper() if clean_name else "X"

                altloc = line[16].strip() if len(line) > 16 else ""
                occ  = _float_field(line, 54, 60, 1.0)
                bfac = _float_field(line, 60, 66, 0.0)

            except (ValueError, IndexError):
                continue  # malformed line — skip

            yield {
                "group":     "HETATM",
                "model_num": model_num,
                "resn":      resn,
                "component_id": resn,
                "auth_comp_id": resn,
                "label_comp_id": resn,
                "chain":     chain,
                "auth_chain_id": chain,
                "label_asym_id": chain,
                "resi":      resi,
                "auth_seq_id": resi,
                "label_seq_id": resi,
                "icode":     icode,
                "x":         x,
                "y":         y,
                "z":         z,
                "element":   element,
                "atom_name": atom_name,
                "auth_atom_id": atom_name,
                "label_atom_id": atom_name,
                "altloc":    altloc,
                "occ":       occ,
                "bfac":      bfac,
            }


def _float_field(line: str, start: int, end: int, default: float) -> float:
    """Parse a PDB numeric column, returning *default* when empty/broken."""
    if len(line) < end:
        return default
    raw = line[start:end].strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


# ════════════════════════════════════════════════════════════════════════════
#  HETATM record parser
# ════════════════════════════════════════════════════════════════════════════

def parse_hetatm_records(structure_file: str | Path) -> List[Dict]:
    """
    Read a receptor structure file (PDB or mmCIF) and return one dict per
    unique HETATM residue group.

    Each dict contains:
        resn               : 3-letter residue name  (upper-case)
        chain              : chain ID
        resi               : residue sequence number + insertion code (string)
        atom_count         : number of HETATM atoms in this group
        heavy_atom_count   : atoms where element != H
        centroid           : (cx, cy, cz) in Ångströms
        is_known_non_ligand: True if in the NON_LIGAND_RESNS reference set
    """
    base_groups: Dict[tuple, List[Dict]] = {}
    structure_path = Path(structure_file)

    if not structure_path.is_file():
        raise FileNotFoundError(
            f"Receptor structure file not found: {structure_file}"
        )

    for atom in iter_hetatm_atoms(structure_path):
        base_key = (
            atom.get("model_num", "1"),
            atom.get("component_id") or atom["resn"].upper(),
            atom.get("auth_chain_id", ""),
            atom.get("label_asym_id", ""),
            atom.get("auth_seq_id", ""),
            atom.get("label_seq_id", ""),
            atom.get("icode", ""),
        )
        base_groups.setdefault(base_key, []).append(atom)

    groups: List[Dict] = []
    for base_key, all_atoms in base_groups.items():
        altlocs = sorted({a.get("altloc", "") for a in all_atoms if a.get("altloc", "")}) or [""]
        for selected_altloc in altlocs:
            atoms = [
                a for a in all_atoms
                if not a.get("altloc", "") or a.get("altloc", "") == selected_altloc
            ]
            occupancies = [float(a.get("occ", 1.0)) for a in atoms]
            occupancy = min(occupancies, default=1.0)
            instance_key = LigandInstanceKey(
                model_num=str(base_key[0]),
                component_id=str(base_key[1]).upper(),
                auth_chain_id=str(base_key[2]),
                label_asym_id=str(base_key[3]),
                auth_seq_id=str(base_key[4]),
                label_seq_id=str(base_key[5]),
                insertion_code=str(base_key[6]),
                altloc=selected_altloc,
                occupancy=occupancy,
            )
            resn = instance_key.component_id
            chain = instance_key.chain_id or "A"
            resi = instance_key.sequence_id
            resi_full = f"{resi}{instance_key.insertion_code}" if instance_key.insertion_code else resi
            group = {
                "resn":                resn,
                "chain":               chain,
                "resi":                resi_full,
                "instance_key":        instance_key.to_dict(),
                "instance_id":         instance_key.display(),
                "atom_count":          0,
                "heavy_atom_count":    0,
                "_coords":             [],
                "is_known_non_ligand": resn in NON_LIGAND_RESNS,
            }
            for atom in atoms:
                group["atom_count"] += 1
                group["_coords"].append((atom["x"], atom["y"], atom["z"]))
                if atom["element"] not in ("H", "D"):
                    group["heavy_atom_count"] += 1
            groups.append(group)

    # Compute centroids and sort largest -> smallest
    result: List[Dict] = []
    for g in groups:
        coords = g.pop("_coords")
        if coords:
            n = len(coords)
            cx = sum(c[0] for c in coords) / n
            cy = sum(c[1] for c in coords) / n
            cz = sum(c[2] for c in coords) / n
            g["centroid"] = (round(cx, 3), round(cy, 3), round(cz, 3))
        else:
            g["centroid"] = (0.0, 0.0, 0.0)
        result.append(g)

    result.sort(key=lambda g: g["heavy_atom_count"], reverse=True)
    return result


def build_docking_target_context(
    structure_file: str | Path,
    instance: LigandInstanceKey | Mapping,
) -> DockingTargetContext:
    """Freeze one exact ligand instance and its labeled crystal heavy atoms."""
    key = instance if isinstance(instance, LigandInstanceKey) else LigandInstanceKey.from_dict(instance)
    atoms: List[LigandAtomReference] = []
    seen_ids: set[str] = set()

    for atom in iter_hetatm_atoms(structure_file):
        if str(atom.get("model_num", "1")) != key.model_num:
            continue
        if (atom.get("component_id") or atom["resn"]).upper() != key.component_id:
            continue
        if str(atom.get("auth_chain_id", "")) != key.auth_chain_id:
            continue
        if str(atom.get("label_asym_id", "")) != key.label_asym_id:
            continue
        if str(atom.get("auth_seq_id", "")) != key.auth_seq_id:
            continue
        if str(atom.get("label_seq_id", "")) != key.label_seq_id:
            continue
        if str(atom.get("icode", "")) != key.insertion_code:
            continue
        atom_altloc = str(atom.get("altloc", ""))
        if atom_altloc and atom_altloc != key.altloc:
            continue
        if atom["element"] in ("H", "D"):
            continue
        atom_id = str(atom.get("label_atom_id") or atom.get("auth_atom_id") or atom["atom_name"])
        if not atom_id or atom_id in seen_ids:
            raise RuntimeError(
                "INVALID_REFERENCE_MAPPING: selected ligand instance has missing or "
                f"duplicate heavy-atom ID '{atom_id}'."
            )
        seen_ids.add(atom_id)
        atoms.append(LigandAtomReference(
            atom_id=atom_id,
            element=str(atom["element"]).upper(),
            xyz=(float(atom["x"]), float(atom["y"]), float(atom["z"])),
            altloc=atom_altloc,
            occupancy=float(atom.get("occ", 1.0)),
        ))

    if not atoms:
        raise RuntimeError(
            "INVALID_LIGAND_INSTANCE: the exact selected ligand instance no longer "
            "exists in the receptor structure. Re-load and re-select the target."
        )
    return DockingTargetContext(str(Path(structure_file).resolve()), key, tuple(atoms))


def build_chain_ligand_map(hetatms: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Group candidate ligands by chain.

    Returns {chain_id: [list of candidate hetatm dicts]} for all
    non-trivial ligands only — i.e. every group that is NOT in
    NON_LIGAND_RESNS is excluded.

    The GUI calls this once at receptor load time to populate its
    chain and ligand selectors.  Each chain's list is sorted by
    heavy-atom count, largest first.
    """
    chain_map: Dict[str, List[Dict]] = {}
    for h in hetatms:
        if h.get("is_known_non_ligand"):
            continue
        chain_map.setdefault(h["chain"], []).append(h)

    for chain in chain_map:
        chain_map[chain].sort(key=lambda g: g["heavy_atom_count"], reverse=True)

    return chain_map
