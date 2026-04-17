"""Protein-change string parser for variant-effect queries.

Accepts common shorthands and normalizes to (original_aa, position, new_aa)
one-letter form.

Accepted forms (case-insensitive, whitespace-tolerant):
  ``R175H``, ``p.R175H``, ``Arg175His``, ``p.Arg175His``.

Raises :class:`ValueError` for unparsable strings — callers surface this
as a user-facing error since garbage input cannot be silently coerced.
"""

from __future__ import annotations

import re

# Three-letter → one-letter amino-acid code table. Covers the 20 standard
# amino acids; non-standard codes (Sec, Pyl, Xaa, etc.) intentionally
# excluded — a variant against those won't be queryable in ClinVar/VEP.
_AA3_TO_1: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
_AA1_TO_3: dict[str, str] = {v: k.title() for k, v in _AA3_TO_1.items()}

_ONE_LETTER_PATTERN = re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z])$")
_THREE_LETTER_PATTERN = re.compile(r"^(?:p\.)?([A-Za-z]{3})(\d+)([A-Za-z]{3})$")


def parse_protein_change(s: str) -> tuple[str, int, str]:
    """Parse a protein-change string into canonical (orig1, pos, new1).

    Args:
        s: Protein change string. Examples:
            ``"R175H"``, ``"p.R175H"``, ``"Arg175His"``, ``"p.Arg175His"``.

    Returns:
        ``(original_aa_one_letter, 1_indexed_position, new_aa_one_letter)``.
        Both amino-acid letters are uppercase, guaranteed to be one of the
        20 standard codes.

    Raises:
        ValueError: if the string doesn't match any recognized form or
            uses a non-standard amino-acid code.
    """
    text = (s or "").strip().replace(" ", "")
    if not text:
        raise ValueError("empty mutation string")

    m1 = _ONE_LETTER_PATTERN.match(text)
    if m1:
        orig, pos, new = m1.group(1), int(m1.group(2)), m1.group(3)
        if orig not in _AA1_TO_3 or new not in _AA1_TO_3:
            raise ValueError(f"unknown amino-acid code in '{s}'")
        return orig, pos, new

    m3 = _THREE_LETTER_PATTERN.match(text)
    if m3:
        orig3, pos, new3 = m3.group(1).upper(), int(m3.group(2)), m3.group(3).upper()
        if orig3 not in _AA3_TO_1 or new3 not in _AA3_TO_1:
            raise ValueError(f"unknown amino-acid code in '{s}'")
        return _AA3_TO_1[orig3], pos, _AA3_TO_1[new3]

    raise ValueError(
        f"unparsable protein change '{s}' — expected forms: 'R175H', "
        "'p.R175H', 'Arg175His', 'p.Arg175His'"
    )


def canonical_one_letter(orig: str, pos: int, new: str) -> str:
    """Return the canonical one-letter form like ``'R175H'``."""
    return f"{orig}{pos}{new}"


def canonical_three_letter(orig: str, pos: int, new: str) -> str:
    """Return the canonical three-letter HGVS form like ``'p.Arg175His'``."""
    return f"p.{_AA1_TO_3[orig]}{pos}{_AA1_TO_3[new]}"
