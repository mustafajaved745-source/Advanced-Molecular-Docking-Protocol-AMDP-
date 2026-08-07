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
from pathlib import Path
from typing import Dict, Iterator, List

from mmcif_parser import is_mmcif, iter_atom_site_rows


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
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
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
                "resn":      resn,
                "chain":     chain,
                "resi":      resi,
                "icode":     icode,
                "x":         x,
                "y":         y,
                "z":         z,
                "element":   element,
                "atom_name": atom_name,
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
    groups: Dict[str, Dict] = {}
    structure_path = Path(structure_file)

    if not structure_path.is_file():
        raise FileNotFoundError(
            f"Receptor structure file not found: {structure_file}"
        )

    for atom in iter_hetatm_atoms(structure_path):
        resn  = atom["resn"].upper()
        chain = atom["chain"] or "A"
        resi  = atom["resi"]
        icode = atom["icode"]
        resi_full = f"{resi}{icode}" if icode else resi
        x, y, z = atom["x"], atom["y"], atom["z"]
        element  = atom["element"]

        key = f"{resn}|{chain}|{resi_full}"
        if key not in groups:
            groups[key] = {
                "resn":                resn,
                "chain":               chain,
                "resi":                resi_full,
                "atom_count":          0,
                "heavy_atom_count":    0,
                "_coords":             [],
                "is_known_non_ligand": resn in NON_LIGAND_RESNS,
            }

        groups[key]["atom_count"] += 1
        groups[key]["_coords"].append((x, y, z))
        if element not in ("H", "D"):
            groups[key]["heavy_atom_count"] += 1

    # Compute centroids and sort largest -> smallest
    result: List[Dict] = []
    for g in groups.values():
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
