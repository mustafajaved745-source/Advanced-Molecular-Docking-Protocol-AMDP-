"""
Molecular Docking Pipeline
==========================
mmcif_parser.py — Minimal, dependency-free mmCIF atom-site parser

mmCIF (PDBx/mmCIF) is the successor of the PDB text format and is the
native format for RCSB and AlphaFold structure downloads.  Unlike PDB,
it is *column-delimited*: the `_atom_site.` loop declares its column
names explicitly, and the data rows follow in that order.  A parser must
therefore map column name -> index rather than assume a fixed layout.

This module provides just enough mmCIF reading for the docking pipeline:

  is_mmcif(path)              -> bool
  iter_atom_site_rows(path)   -> Iterator[dict]   (one dict per atom row)

The atom dicts expose the fields the rest of the pipeline needs to build
the active-site panel, the docking grid box, and the RMSD reference —
without ever writing an intermediate PDB file.

Only the first model is kept (`_atom_site.pdbx_PDB_model_num`).  The model
number and both label/auth identifiers are nevertheless retained on every
row so callers can construct an exact, auditable ligand-instance key.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterator, List, Optional


# ════════════════════════════════════════════════════════════════════════════
#  Format detection
# ════════════════════════════════════════════════════════════════════════════

def is_mmcif(path: str | Path) -> bool:
    """
    Return True if *path* is an mmCIF structure file.

    Decision by extension (`.cif`, `.mmcif`) with a content sniff as a
    fallback for files with an unexpected or missing extension: mmCIF
    files always start with a `data_` block name.
    """
    suffix = Path(path).suffix.lower()
    if suffix in (".cif", ".mmcif"):
        return True

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                return line.startswith("data_")
    except OSError:
        return False
    return False


# ════════════════════════════════════════════════════════════════════════════
#  Tokenizer
# ════════════════════════════════════════════════════════════════════════════

def _tokenize_line(line: str) -> List[str]:
    """
    Split one mmCIF line into values.

    - Whitespace separates unquoted values.
    - A value beginning with `'` or `"` runs to the matching quote and
      may contain spaces (an unterminated quote swallows the rest of the
      line — a safety net, this never occurs inside `_atom_site` loops).
    - `#` starts a comment that runs to the end of the line.
    """
    tokens: List[str] = []
    i, n = 0, len(line)

    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c == "#":
            break  # comment — ignore the remainder of the line
        if c in ("'", '"'):
            quote = c
            buf: List[str] = []
            i += 1
            while i < n and line[i] != quote:
                buf.append(line[i])
                i += 1
            tokens.append("".join(buf))
            i += 1  # skip closing quote (or run past end of line)
            continue

        j = i
        while j < n and line[j] not in " \t":
            j += 1
        tokens.append(line[i:j])
        i = j

    return tokens


def _token_stream(path: Path) -> Iterator[str]:
    """Yield every mmCIF token in the file, comments removed."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            yield from _tokenize_line(line)


# ════════════════════════════════════════════════════════════════════════════
#  Loop parser
# ════════════════════════════════════════════════════════════════════════════

_ATOM_SITE_LOOP_KEY = "_atom_site.group_PDB"


def _normalize(value: str) -> str:
    """Map an mmCIF '.' / '?' missing-value marker to an empty string."""
    v = value.strip()
    return "" if v in (".", "?") else v


def iter_atom_site_rows(path: str | Path) -> Iterator[Dict]:
    """
    Yield one normalized dict per `_atom_site` atom row in *path*.

    Only the first model is emitted (`pdbx_PDB_model_num` not in
    {"", ".", "?", "1"} is skipped).  Keys are the fields the pipeline
    consumes downstream:

        group     : "ATOM" or "HETATM"   (_atom_site.group_PDB)
        resn      : 3-letter residue / ligand code (auth_comp_id or label_comp_id)
        chain     : chain ID (auth_asym_id or label_asym_id)
        resi      : residue sequence number (auth_seq_id or label_seq_id)
        icode     : PDB insertion code (pdbx_PDB_ins_code)
        x, y, z   : Cartesian coordinates in Ångströms
        element   : element symbol, upper-cased
        atom_name : atom label (label_atom_id)
        altloc    : alternate-location ID (label_alt_id)
        occ       : occupancy (default 1.0)
        bfac      : B-factor (default 0.0)
    """
    header: List[str] = []      # column names of the current loop
    values: List[str] = []      # buffered data values of the current loop
    phase = "top"               # "top" | "header" | "data"
    atom_loop = False           # is the current loop the _atom_site loop?

    for tok in _token_stream(Path(path)):
        if tok == "loop_":
            yield from _flush(atom_loop, header, values)
            header, values = [], []
            phase = "header"
            atom_loop = False
            continue

        if phase == "top":
            continue  # top-level scalar (data_ name, etc.) — not needed

        if phase == "header":
            if tok.startswith("_"):
                header.append(tok)
            else:
                # First data value of the loop → data phase begins.
                atom_loop = any(c == _ATOM_SITE_LOOP_KEY for c in header)
                values.append(tok)
                phase = "data"
            continue

        # phase == "data"
        if tok.startswith("_"):
            # A new non-loop category line: the previous loop has ended.
            yield from _flush(atom_loop, header, values)
            header, values = [tok], []
            atom_loop = False
            continue

        values.append(tok)
        if len(header) == 1 and not atom_loop:
            # Single-tag assignment (e.g. `_cell.length_a 12.3`): one value.
            yield from _flush(atom_loop, header, values)
            header, values = [], []
            phase = "top"

    yield from _flush(atom_loop, header, values)


def _flush(
    atom_loop: bool,
    header: List[str],
    values: List[str],
) -> Iterator[Dict]:
    """Emit normalized atom dicts for a completed _atom_site loop chunk."""
    if not atom_loop or not header or not values:
        return

    ncols = len(header)
    full = len(values) - (len(values) % ncols)  # drop a partial trailing row

    for i in range(0, full, ncols):
        row = dict(zip(header, values[i:i + ncols]))
        model = _normalize(row.get("_atom_site.pdbx_PDB_model_num", ""))
        if model not in ("", "1"):
            continue  # keep only the first model

        atom = _normalize_row(row)
        if atom is not None:
            yield atom


def _normalize_row(row: Dict[str, str]) -> Optional[Dict]:
    """Map a raw `_atom_site` row to the pipeline's normalized atom dict."""
    group = _normalize(row.get("_atom_site.group_PDB", "")).upper()
    if not group:
        return None

    auth_comp_id = _normalize(row.get("_atom_site.auth_comp_id", "")).upper()
    label_comp_id = _normalize(row.get("_atom_site.label_comp_id", "")).upper()
    resn = auth_comp_id or label_comp_id

    auth_chain_id = _normalize(row.get("_atom_site.auth_asym_id", ""))
    label_asym_id = _normalize(row.get("_atom_site.label_asym_id", ""))
    chain = auth_chain_id or label_asym_id

    auth_seq_id = _normalize(row.get("_atom_site.auth_seq_id", ""))
    label_seq_id = _normalize(row.get("_atom_site.label_seq_id", ""))
    resi = auth_seq_id or label_seq_id

    icode = _normalize(row.get("_atom_site.pdbx_PDB_ins_code", ""))

    element = _normalize(row.get("_atom_site.type_symbol", "")).upper()
    label_atom_id = _normalize(row.get("_atom_site.label_atom_id", ""))
    auth_atom_id = _normalize(row.get("_atom_site.auth_atom_id", ""))
    atom_name = label_atom_id or auth_atom_id
    if not element and atom_name:
        # Fallback: derive the element from the leading alphabetic run of
        # the atom name (e.g. "CA" -> C, "1HE2" -> H).
        element = "".join(c for c in atom_name if c.isalpha())[:1].upper()

    try:
        x = float(row.get("_atom_site.Cartn_x", ""))
        y = float(row.get("_atom_site.Cartn_y", ""))
        z = float(row.get("_atom_site.Cartn_z", ""))
    except (TypeError, ValueError):
        return None  # malformed coordinate row — skip

    def _float_or(key: str, default: float) -> float:
        v = _normalize(row.get(key, ""))
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    return {
        "group":     group,
        "model_num": _normalize(row.get("_atom_site.pdbx_PDB_model_num", "")) or "1",
        "resn":      resn,
        "component_id": label_comp_id or auth_comp_id,
        "auth_comp_id": auth_comp_id,
        "label_comp_id": label_comp_id,
        "chain":     chain,
        "auth_chain_id": auth_chain_id,
        "label_asym_id": label_asym_id,
        "resi":      resi,
        "auth_seq_id": auth_seq_id,
        "label_seq_id": label_seq_id,
        "icode":     icode,
        "x":         x,
        "y":         y,
        "z":         z,
        "element":   element,
        "atom_name": atom_name,
        "auth_atom_id": auth_atom_id,
        "label_atom_id": label_atom_id,
        "altloc":    _normalize(row.get("_atom_site.label_alt_id", "")),
        "occ":       _float_or("_atom_site.occupancy", 1.0),
        "bfac":      _float_or("_atom_site.B_iso_or_equiv", 0.0),
    }
