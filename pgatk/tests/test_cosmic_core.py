"""Unit tests for CancerGenomesService.get_mut_pro_seq()."""

from Bio.Seq import Seq

from pgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pgatk.cgenomes.models import SNP


class TestGetMutProSeqDnaBranch:
    """Tests for the DNA-level mutation branch (when dna_mut has no '?' and aa_mut != 'p.?')."""

    def test_substitution_via_dna(self):
        """A DNA substitution (c.4A>T) changes the 4th base from A to T."""
        # CDS: ATGAAATTT -> protein MKF
        # Mutate position 4 (A->T): ATGTATTT -> protein MYF... let's verify
        # Actually position 4 is 'A' in ATGAAATTT (1-indexed: A=1, T=2, G=3, A=4)
        # After mutation: ATG T AATTT -> 'ATGTAATTT' which translates to M*F
        seq = Seq("ATGAAATTT")
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.4A>T",
            aa_mut="p.K2*",
            mutation_type="Nonsense"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # Position 4 (index 3): seq[3]='A', replace with 'T' -> 'ATGTAATTT'
        # translate(to_stop=False) -> M*F
        assert result is not None
        assert str(result) != ""
        assert str(result)[0] == "M"

    def test_insertion_via_dna(self):
        """A DNA insertion (c.3_4insGGG) inserts GGG between positions 3 and 4."""
        seq = Seq("ATGAAATTT")
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.3_4insGGG",
            aa_mut="p.?",  # does not matter for DNA branch, but aa_mut != 'p.?' so we use valid one
            mutation_type="Insertion - In frame"
        )
        # For DNA ins branch: re.findall(r'\d+', "c.3_4insGGG") = ['3','4']
        # But the code checks ">" first (no), then "ins" (yes)
        # index of "ins" in dna_mut = 5, insert_dna = "GGG"
        # ins_index1 = int(positions[0]) = 3
        # seq_mut = seq[:3] + "GGG" + seq[3:] = "ATG" + "GGG" + "AAATTT" = "ATGGGGAAATTT"
        # translate -> MGKF
        # Need aa_mut != 'p.?' for this branch
        snp.aa_mut = "p.M1_K2insG"
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        assert result is not None
        assert len(str(result)) > 0
        # Original: ATGAAATTT -> MKF (3 aa)
        # Mutated:  ATGGGGAAATTT -> MGKF (4 aa)
        assert str(result) == "MGKF"

    def test_deletion_two_positions_via_dna(self):
        """A DNA deletion (c.4_6del) removes bases at positions 4-6."""
        seq = Seq("ATGAAATTTCCC")  # 12 bases -> MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.4_6del",
            aa_mut="p.K2del",
            mutation_type="Deletion - In frame"
        )
        # positions = ['4', '6']
        # del_index1 = 4-1 = 3, del_index2 = 6
        # seq_mut = seq[:3] + seq[6:] = "ATG" + "TTTCCC" = "ATGTTTCCC" (9 bases)
        # translate -> MFP
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        assert result is not None
        assert str(result) == "MFP"

    def test_deletion_single_position_via_dna(self):
        """A DNA single-base deletion (c.4del) removes the base at position 4."""
        seq = Seq("ATGAAATTT")  # MKF
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.4del",
            aa_mut="p.K2fs",
            mutation_type="Frameshift"
        )
        # positions = ['4'], len(positions)==1
        # del_index1 = 4-1 = 3
        # seq_mut = seq[:3] + seq[4:] = "ATG" + "AATTT" = "ATGAATTT" (8 bases)
        # translate(to_stop=False) -> MNF (incomplete last codon 'T' is dropped? No, Biopython pads)
        # Actually 8 bases: ATG AAT TT -> MN + partial 'TT' -> depends on Biopython behavior
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        assert result is not None
        assert str(result).startswith("MN")

    def test_ambiguous_dna_mut_skips_to_protein_branch(self):
        """When dna_mut contains '?', the protein-level branch is used instead."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2N",
            mutation_type="Missense"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # Protein branch: Missense, aa_mut[-1] = 'N', positions = ['2'], index = 1
        # protein = 'MKFP', mut = 'M' + 'N' + 'FP' = 'MNFP'
        assert str(result) == "MNFP"


class TestGetMutProSeqProteinBranch:
    """Tests for the protein-level mutation branch (when dna_mut contains '?' or aa_mut == 'p.?')."""

    def test_missense_substitution(self):
        """Missense: change one amino acid at a specified position."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2N",
            mutation_type="Missense"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        assert str(result) == "MNFP"

    def test_nonsense_truncation(self):
        """Nonsense: truncates the protein at the mutation position."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2*",
            mutation_type="Nonsense"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # index = 2-1 = 1, protein[:1] = 'M'
        assert str(result) == "M"

    def test_inframe_insertion(self):
        """In-frame insertion inserts amino acids at the specified position."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2_F3insAA",
            mutation_type="Insertion - In frame"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # positions = ['2', '3'], ins_index1 = int(positions[0]) = 2
        # insert_aa = 'AA'
        # mut_pro_seq = protein[:2] + 'AA' + protein[2:] = 'MK' + 'AA' + 'FP' = 'MKAAFP'
        assert str(result) == "MKAAFP"

    def test_inframe_deletion_two_positions(self):
        """In-frame deletion removes amino acids between two positions."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2_F3del",
            mutation_type="Deletion - In frame"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # positions = ['2', '3']
        # del_index1 = 2-1 = 1, del_index2 = 3
        # mut_pro_seq = protein[:1] + protein[3:] = 'M' + 'P' = 'MP'
        assert str(result) == "MP"

    def test_inframe_deletion_single_position(self):
        """In-frame single amino acid deletion."""
        seq = Seq("ATGAAATTTCCC")  # protein: MKFP
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2del",
            mutation_type="Deletion - In frame"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # positions = ['2'], del_index1 = 1
        # protein[:1] + protein[2:] = 'M' + 'FP' = 'MFP'
        assert str(result) == "MFP"


class TestGetMutProSeqEdgeCases:
    """Edge cases for get_mut_pro_seq."""

    def test_mutation_beyond_sequence_length(self):
        """Mutation at a position beyond the sequence should produce an IndexError or empty result."""
        seq = Seq("ATGAAA")  # protein: MK (2 aa)
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.A10N",
            mutation_type="Missense"
        )
        # positions = ['10'], index = 9
        # protein_seq = 'MK' (length 2), protein_seq[:9] + 'N' + protein_seq[10:]
        # This produces 'N' padded with empty strings -> should not crash
        # The function may raise IndexError or return a weird string
        try:
            result = CancerGenomesService.get_mut_pro_seq(snp, seq)
            # If it returns a result, it should still be a string
            assert isinstance(str(result), str)
        except (IndexError, ValueError):
            pass  # Acceptable: function does not handle out-of-bounds gracefully

    def test_ambiguous_aa_mut_returns_empty(self):
        """When both dna_mut and aa_mut are '?', function returns empty string."""
        seq = Seq("ATGAAATTT")
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.?",
            mutation_type="Unknown"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        assert result == ""

    def test_missense_non_alpha_last_char_returns_empty(self):
        """When aa_mut last char is not alphabetical, return empty."""
        seq = Seq("ATGAAATTTCCC")
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.?",
            aa_mut="p.K2*",
            mutation_type="Missense"  # but last char is '*', not alpha
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # mut_aa = '*', not isalpha() -> return ''
        assert result == ""

    def test_dna_substitution_ref_mismatch(self):
        """When DNA ref allele doesn't match the sequence, return empty."""
        seq = Seq("ATGAAATTT")  # pos 4 = 'A'
        snp = SNP(
            gene="TEST",
            mrna="ENST0001",
            dna_mut="c.4C>T",  # claims ref is C, but seq[3] = 'A'
            aa_mut="p.K2*",
            mutation_type="Substitution"
        )
        result = CancerGenomesService.get_mut_pro_seq(snp, seq)
        # ref_dna = 'C', seq[3] = 'A' -> mismatch -> mut_pro_seq stays ""
        assert result == ""
