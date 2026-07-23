"""
AI-Guided Molecular Docking Pipeline
=====================================
ai_identify.py — Active-site ligand identification via any OpenAI-compatible API

Responsibilities:
  1. Parse all HETATM records from a raw receptor PDB
  2. Filter out water, ions, and common buffer / cryoprotectant molecules
  3. Send the remaining candidates to a configurable language model API
  4. Parse the response and return the 3-letter residue code of the
     biologically relevant active-site ligand

Public API
----------
  parse_hetatm_records(pdb_file)              -> List[dict]
  identify_active_site_ligand(hetatms, ...)   -> (code: str, reason: str)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple


# ════════════════════════════════════════════════════════════════════════════
#  Reference sets
# ════════════════════════════════════════════════════════════════════════════

# Everything in this set is treated as a definite non-ligand and filtered
# before the API call.  The AI only sees candidates that survive this filter.
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

    # ── Amino-acid fragments / protecting groups ──────────────────────────
    "NME",
})

# Residue names that are almost certainly ligands regardless of size
KNOWN_DRUG_HINTS: frozenset[str] = frozenset({
    # Common well-known drugs / cofactors found in PDB structures
    "ATP", "ADP", "AMP", "GTP", "GDP", "NAD", "FAD", "FMN", "COA",
    "HEM", "HEC", "HEA",                           # haem variants
    "STI", "IMA",                                   # imatinib
    "TYR", "PHE", "TRP",                            # free amino acids as ligands
})


# ════════════════════════════════════════════════════════════════════════════
#  HETATM record parser
# ════════════════════════════════════════════════════════════════════════════

def parse_hetatm_records(pdb_file: str | Path) -> List[Dict]:
    """
    Read a PDB file and return one dict per unique HETATM residue group.

    Each dict contains:
        resn               : 3-letter residue name  (upper-case)
        chain              : chain ID
        resi               : residue sequence number + insertion code (string)
        atom_count         : number of HETATM atoms in this group
        heavy_atom_count   : atoms where element != H
        centroid           : (cx, cy, cz) in Ångströms
        is_known_non_ligand: True if in the NON_LIGAND_RESNS reference set
        is_drug_hint       : True if in KNOWN_DRUG_HINTS
    """
    groups: Dict[str, Dict] = {}
    pdb_path = Path(pdb_file)

    if not pdb_path.is_file():
        raise FileNotFoundError(f"Receptor PDB file not found: {pdb_file}")

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue

            try:
                resn  = line[17:20].strip().upper()
                chain = line[21].strip() or "A"
                resi  = line[22:26].strip()
                icode = line[26].strip()
                resi_full = f"{resi}{icode}" if icode else resi

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                element = ""
                if len(line) >= 78:
                    element = line[76:78].strip().upper()

                if not element:
                    # Fallback parsing for element from atom name column (12-16)
                    atom_name = line[12:16].strip()
                    # Strip leading numbers (e.g., "1HE2" -> "HE2")
                    clean_name = re.sub(r"^\d+", "", atom_name)
                    element = clean_name[:1].upper() if clean_name else "X"

            except (ValueError, IndexError):
                continue  # malformed line — skip

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
                    "is_drug_hint":        resn in KNOWN_DRUG_HINTS,
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


# ════════════════════════════════════════════════════════════════════════════
#  API Call Logic
# ════════════════════════════════════════════════════════════════════════════

def identify_active_site_ligand(
    hetatms:  List[Dict],
    api_key:  str,
    model:    str = "meta-llama/llama-3.1-70b-instruct",
    base_url: str = "https://openrouter.ai/api/v1",
) -> Tuple[str, str]:
    """
    Ask any OpenAI-compatible API which HETATM group is the
    biologically relevant active-site ligand.

    Returns
    -------
    (ligand_code, reasoning)
        ligand_code : 3-letter (or up to 5-letter) PDB residue name
        reasoning   : one-sentence AI explanation
    """

    # ── Separate candidates from definite non-ligands ────────────────────
    candidates = [h for h in hetatms if not h["is_known_non_ligand"]]
    excluded   = [h for h in hetatms if h["is_known_non_ligand"]]

    # ── Trivial cases: skip the API ──────────────────────────────────────
    if not candidates:
        raise ValueError(
            "No candidate ligands remain after filtering water, ions, and "
            "common buffer molecules.\n"
            "The receptor PDB may not contain a co-crystallised drug-like ligand.\n"
            "Check the HETATM records with a molecular viewer."
        )

    if len(candidates) == 1:
        c = candidates[0]
        return (
            c["resn"],
            f"Only one candidate ligand found after filtering: "
            f"{c['resn']} ({c['heavy_atom_count']} heavy atoms, "
            f"chain {c['chain']}, residue {c['resi']}).",
        )

    # ── Check for an unambiguous drug hint ───────────────────────────────
    hint_candidates = [c for c in candidates if c["is_drug_hint"]]
    if len(hint_candidates) == 1:
        c = hint_candidates[0]
        return (
            c["resn"],
            f"{c['resn']} is a known cofactor / drug molecule "
            f"({c['heavy_atom_count']} heavy atoms).  "
            "Identified without API call.",
        )

    # ── Build HETATM summary for the prompt ──────────────────────────────
    summary_lines = [
        f"{'ResName':<8}  {'Chain':<6}  {'ResID':<6}  "
        f"{'HeavyAtoms':>10}  {'Centroid (Å)'}",
        "-" * 65,
    ]
    for h in candidates:
        cx, cy, cz = h["centroid"]
        note = "  ← known drug/cofactor" if h["is_drug_hint"] else ""
        summary_lines.append(
            f"{h['resn']:<8}  {h['chain']:<6}  {h['resi']:<6}  "
            f"{h['heavy_atom_count']:>10}  "
            f"({cx:.1f}, {cy:.1f}, {cz:.1f}){note}"
        )

    if excluded:
        excl_names = sorted({h["resn"] for h in excluded})
        summary_lines.append(
            f"\n(Pre-filtered as non-ligands: {', '.join(excl_names)})"
        )

    hetatm_summary = "\n".join(summary_lines)

    # ── Prompt ────────────────────────────────────────────────────────────
    prompt = f"""You are an expert structural biologist and medicinal chemist \
analysing a protein crystal structure for virtual screening.

The following HETATM groups remain after removing water, common ions, \
and known buffer / cryoprotectant molecules:

{hetatm_summary}

Your task
─────────
Identify the ONE residue that is the biologically relevant co-crystallised \
ligand in the protein's active site — the molecule a medicinal chemist \
would use as the reference for docking.

Selection criteria (priority order):
  1. Drug-like size: prefer molecules with 5–70 heavy atoms.
  2. Chemical nature: organic, drug-like, or known cofactor / substrate.
  3. Context: active-site / binding-pocket binders over crystal contacts.
  4. Exclude: coordination metals, small fragments, residual additives.

IMPORTANT: Respond in EXACTLY this format — no other text:
LIGAND_CODE: XXX
REASONING: One sentence explanation.

Where XXX is the exact 3-letter (or up to 5-letter) PDB residue code \
from the table above."""

    # ── API call ─────────────────────────────────────────────────────────
    try:
        from openai import OpenAI, APIError
    except ImportError:
        raise ImportError(
            "The 'openai' Python package is required for API calls.\n"
            "Install it with:  pip install openai"
        )

    client = OpenAI(
        base_url=base_url,
        api_key=api_key,
    )

    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=350,
            temperature=0.05,       # near-deterministic for reproducibility
        )
    except APIError as err:
        raise RuntimeError(
            f"LLM API call failed ({err.__class__.__name__}): {err}"
        ) from err

    raw_text = response.choices[0].message.content.strip() if response.choices else ""

    if not raw_text:
        # Emergency fallback if API returns empty content
        biggest = candidates[0]
        return (
            biggest["resn"],
            f"LLM returned an empty response; falling back to largest candidate: {biggest['resn']}."
        )

    # ── Parse response ────────────────────────────────────────────────────
    code, reason = _parse_api_response(raw_text, candidates)
    return code, reason


# ════════════════════════════════════════════════════════════════════════════
#  Response parser
# ════════════════════════════════════════════════════════════════════════════

def _parse_api_response(
    text: str,
    candidates: List[Dict],
) -> Tuple[str, str]:
    """
    Extract LIGAND_CODE and REASONING from the model's reply.
    Falls back gracefully if the model doesn't follow the exact format.
    """
    candidate_codes = {h["resn"].upper() for h in candidates}
    code   = ""
    reason = ""

    # ── Try structured regex extraction ──────────────────────────────────
    code_match = re.search(r"LIGAND_CODE:\s*([A-Za-z0-9]+)", text, re.IGNORECASE)
    if code_match:
        code = code_match.group(1).upper()

    reason_match = re.search(r"REASONING:\s*(.+)", text, re.IGNORECASE)
    if reason_match:
        reason = reason_match.group(1).strip()

    # If regex missed or line-by-line is preferred:
    if not code or not reason:
        for line in text.splitlines():
            stripped = line.strip()
            upper    = stripped.upper()

            if not code and upper.startswith("LIGAND_CODE:"):
                raw_code = stripped.split(":", 1)[1].strip().upper()
                raw_code = re.sub(r"[^A-Z0-9]", "", raw_code)
                if raw_code:
                    code = raw_code

            elif not reason and upper.startswith("REASONING:"):
                reason = stripped.split(":", 1)[1].strip()

    # ── Validate the extracted code ───────────────────────────────────────
    if code and code in candidate_codes:
        return code, reason or text[:300]

    # ── Fallback 1: look for any candidate code mentioned in the reply ────
    text_upper = text.upper()
    for c in sorted(candidates,
                    key=lambda h: h["heavy_atom_count"],
                    reverse=True):
        if re.search(rf"\b{re.escape(c['resn'].upper())}\b", text_upper):
            return (
                c["resn"],
                f"Extracted from AI response (fallback match). "
                f"AI said: {text[:250]}",
            )

    # ── Fallback 2: if code looks valid but wasn't in candidate list ──────
    if code and 2 <= len(code) <= 5 and code.isalnum():
        return (
            code,
            f"AI returned [{code}] which was not in the candidate list. "
            f"Proceeding — verify manually. AI said: {text[:200]}",
        )

    # ── Last resort: return the largest candidate ─────────────────────────
    biggest = candidates[0]
    return (
        biggest["resn"],
        f"Could not parse AI response reliably; defaulting to the largest "
        f"candidate: {biggest['resn']} ({biggest['heavy_atom_count']} heavy atoms). "
        f"AI said: {text[:200]}",
    )