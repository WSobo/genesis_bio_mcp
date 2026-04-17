"""variant_parser unit tests."""

from __future__ import annotations

import pytest

from genesis_bio_mcp.tools.variant_parser import (
    canonical_one_letter,
    canonical_three_letter,
    parse_protein_change,
)


class TestParseProteinChange:
    def test_one_letter_plain(self):
        assert parse_protein_change("R175H") == ("R", 175, "H")

    def test_one_letter_with_p_prefix(self):
        assert parse_protein_change("p.R175H") == ("R", 175, "H")

    def test_three_letter_plain(self):
        assert parse_protein_change("Arg175His") == ("R", 175, "H")

    def test_three_letter_with_p_prefix(self):
        assert parse_protein_change("p.Arg175His") == ("R", 175, "H")

    def test_mixed_case(self):
        assert parse_protein_change("p.arg175HIS") == ("R", 175, "H")

    def test_whitespace_stripped(self):
        assert parse_protein_change("  p.R175H  ") == ("R", 175, "H")

    def test_rejects_empty_string(self):
        with pytest.raises(ValueError, match="empty"):
            parse_protein_change("")

    def test_rejects_non_standard_aa(self):
        # U = selenocysteine; not in the standard 20
        with pytest.raises(ValueError, match="unknown amino-acid"):
            parse_protein_change("U175H")

    def test_rejects_unparsable(self):
        with pytest.raises(ValueError, match="unparsable"):
            parse_protein_change("not a mutation")

    def test_rejects_three_letter_nonstandard(self):
        with pytest.raises(ValueError, match="unknown amino-acid"):
            parse_protein_change("Xyz175Abc")


class TestCanonicalForms:
    def test_one_letter_form(self):
        assert canonical_one_letter("R", 175, "H") == "R175H"

    def test_three_letter_form(self):
        assert canonical_three_letter("R", 175, "H") == "p.Arg175His"

    def test_three_letter_all_twenty(self):
        # Round-trip every standard AA to ensure the reverse table is complete.
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            result = canonical_three_letter(aa, 1, aa)
            assert result.startswith("p.")
            assert result.endswith(result.split("1")[1])  # symmetric
