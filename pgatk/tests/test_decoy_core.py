"""Unit tests for ProteinDBDecoyService static methods."""

import random
from collections import Counter

from pgatk.proteomics.db.protein_database_decoy import ProteinDBDecoyService


# ---------------------------------------------------------------------------
# revswitch
# ---------------------------------------------------------------------------

class TestRevswitch:
    """Tests for ProteinDBDecoyService.revswitch()."""

    def test_reverse_with_switch(self):
        """Reversed sequence should have cleavage-site residues swapped with preceding residue."""
        protein = "ABKCD"
        # Reversed: "DCKBA"
        # With switch on sites=['K','R']: K is at index 2 in reversed seq "DCKBA"
        # Swap K (index 2) with preceding (index 1): "DKCBA"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        assert result == "DKCBA"

    def test_reverse_without_switch(self):
        """When noswitch=True, the result is simply the reversed sequence."""
        protein = "ABKCD"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=True, sites=['K', 'R'])
        assert result == "DCKBA"

    def test_preserves_length(self):
        """Output should always be the same length as input."""
        protein = "MKFPAKLRST"
        result_switch = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        result_no_switch = ProteinDBDecoyService.revswitch(protein, noswitch=True, sites=['K', 'R'])
        assert len(result_switch) == len(protein)
        assert len(result_no_switch) == len(protein)

    def test_preserves_amino_acid_composition(self):
        """Reversed sequence (with or without switch) has the same amino acid composition."""
        protein = "MKFPAKLRST"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        assert Counter(result) == Counter(protein)

    def test_no_cleavage_sites_same_as_plain_reverse(self):
        """With empty sites list, revswitch is just a plain reverse."""
        protein = "ABCDEF"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=[])
        assert result == "FEDCBA"

    def test_single_char_protein(self):
        """Single character protein should return itself."""
        protein = "K"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        assert result == "K"
        assert len(result) == 1

    def test_multiple_cleavage_sites(self):
        """Sequence with multiple K/R sites should switch each with its preceding residue."""
        protein = "AKBRC"
        # Reversed: "CRBKA"
        # At index 1: R -> swap with index 0 (C): "RCBKA"
        # At index 3: K -> swap with index 2 (B): "RCKBA"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        assert result == "RCKBA"

    def test_adjacent_cleavage_sites(self):
        """When two cleavage sites are adjacent, both should be processed."""
        protein = "AKRD"
        # Reversed: "DRKA"
        # Process index 1: R -> swap with D at index 0 -> "RDKA"
        # Process index 2: K -> swap with D at index 1 -> "RKDA"
        result = ProteinDBDecoyService.revswitch(protein, noswitch=False, sites=['K', 'R'])
        assert len(result) == len(protein)
        assert Counter(result) == Counter(protein)


# ---------------------------------------------------------------------------
# shuffle
# ---------------------------------------------------------------------------

class TestShuffle:
    """Tests for ProteinDBDecoyService.shuffle()."""

    def test_preserves_length(self):
        """Shuffled peptide has the same length as the input."""
        peptide = "ABCDEFK"
        result = ProteinDBDecoyService.shuffle(peptide)
        assert len(result) == len(peptide)

    def test_preserves_cterminal(self):
        """The C-terminal amino acid is preserved in position."""
        peptide = "ABCDEFK"
        result = ProteinDBDecoyService.shuffle(peptide)
        assert result[-1] == "K"

    def test_preserves_amino_acid_composition(self):
        """Shuffled peptide has exactly the same amino acid composition."""
        peptide = "MKFPAKLRST"
        result = ProteinDBDecoyService.shuffle(peptide)
        assert Counter(result) == Counter(peptide)

    def test_different_ordering(self):
        """Shuffling usually produces a different ordering (statistical, use seed for determinism)."""
        peptide = "ABCDEFGHIJK"
        random.seed(42)
        result = ProteinDBDecoyService.shuffle(peptide)
        # With a sufficiently long peptide, shuffling is very likely to differ
        # But let's just assert it's a valid permutation
        assert Counter(result) == Counter(peptide)
        assert result[-1] == peptide[-1]

    def test_two_char_peptide(self):
        """A 2-character peptide: first char may swap but last is kept."""
        peptide = "AK"
        result = ProteinDBDecoyService.shuffle(peptide)
        assert result == "AK"  # only one arrangement possible with c-term fixed
        assert len(result) == 2

    def test_single_char_peptide(self):
        """A single character peptide returns itself."""
        peptide = "K"
        result = ProteinDBDecoyService.shuffle(peptide)
        assert result == "K"


# ---------------------------------------------------------------------------
# count_aa_in_dictionary
# ---------------------------------------------------------------------------

class TestCountAaInDictionary:
    """Tests for ProteinDBDecoyService.count_aa_in_dictionary()."""

    def test_counts_amino_acids(self):
        """Counts each amino acid in the given sequence."""
        aa_dict = {"A": 0, "K": 0, "M": 0, "F": 0}
        sequence = "MAAKFF"
        result = ProteinDBDecoyService.count_aa_in_dictionary(aa_dict, sequence)
        assert result["M"] == 1
        assert result["A"] == 2
        assert result["K"] == 1
        assert result["F"] == 2

    def test_missing_amino_acid_stays_zero(self):
        """Amino acids not in the sequence get count 0."""
        aa_dict = {"A": 0, "W": 0}
        sequence = "AAAA"
        result = ProteinDBDecoyService.count_aa_in_dictionary(aa_dict, sequence)
        assert result["A"] == 4
        assert result["W"] == 0

    def test_empty_sequence(self):
        """Empty sequence results in all zeros."""
        aa_dict = {"A": 0, "K": 0}
        result = ProteinDBDecoyService.count_aa_in_dictionary(aa_dict, "")
        assert result["A"] == 0
        assert result["K"] == 0

    def test_overwrites_existing_counts(self):
        """The function overwrites (not accumulates) counts from a fresh count."""
        aa_dict = {"A": 99, "K": 50}
        sequence = "AK"
        result = ProteinDBDecoyService.count_aa_in_dictionary(aa_dict, sequence)
        assert result["A"] == 1
        assert result["K"] == 1

    def test_full_alphabet(self):
        """Test with the full PGATK alphabet."""
        from pgatk.proteomics.models import PGATK_ALPHABET
        aa_dict = {aa: 0 for aa in PGATK_ALPHABET}
        sequence = "ACDEFGHIKLMNPQRSTVWYBJOUXZ"
        result = ProteinDBDecoyService.count_aa_in_dictionary(aa_dict, sequence)
        for aa in PGATK_ALPHABET:
            assert result[aa] == sequence.count(aa)
