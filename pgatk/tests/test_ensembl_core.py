"""Unit tests for EnsemblDataService core static methods."""

from Bio.Seq import Seq

from pgatk.ensembl.ensembl import EnsemblDataService


# ---------------------------------------------------------------------------
# get_multiple_options
# ---------------------------------------------------------------------------

class TestGetMultipleOptions:
    """Tests for EnsemblDataService.get_multiple_options()."""

    def test_comma_separated_values(self):
        result = EnsemblDataService.get_multiple_options("option1, option2, option3")
        assert result == ["option1", "option2", "option3"]

    def test_single_value(self):
        result = EnsemblDataService.get_multiple_options("only_one")
        assert result == ["only_one"]

    def test_strips_whitespace(self):
        result = EnsemblDataService.get_multiple_options("  a , b , c  ")
        assert result == ["a", "b", "c"]

    def test_empty_string(self):
        result = EnsemblDataService.get_multiple_options("")
        assert result == [""]


# ---------------------------------------------------------------------------
# get_key
# ---------------------------------------------------------------------------

class TestGetKey:
    """Tests for EnsemblDataService.get_key()."""

    def test_pipe_separated_header(self):
        assert EnsemblDataService.get_key("ENST00000123|extra") == "ENST00000123"

    def test_space_separated_header(self):
        assert EnsemblDataService.get_key("ENST00000123 some description") == "ENST00000123"

    def test_pipe_before_space(self):
        assert EnsemblDataService.get_key("ENST00000123|gene description") == "ENST00000123"

    def test_no_separator(self):
        assert EnsemblDataService.get_key("ENST00000123") == "ENST00000123"


# ---------------------------------------------------------------------------
# check_overlap
# ---------------------------------------------------------------------------

class TestCheckOverlap:
    """Tests for EnsemblDataService.check_overlap()."""

    def test_variant_within_feature(self):
        """Variant fully inside a feature region."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(120, 130, features) is True

    def test_variant_outside_feature(self):
        """Variant completely outside all features."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(50, 60, features) is False

    def test_variant_overlaps_feature_start(self):
        """Variant overlaps the start boundary of a feature."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(90, 110, features) is True

    def test_variant_overlaps_feature_end(self):
        """Variant overlaps the end boundary of a feature."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(190, 210, features) is True

    def test_variant_covers_entire_feature(self):
        """Variant fully covers the feature."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(50, 300, features) is True

    def test_start_minus_one_always_true(self):
        """When var_start is -1, overlap is always True."""
        features = [[100, 200, 'exon']]
        assert EnsemblDataService.check_overlap(-1, 10, features) is True

    def test_multiple_features_overlap_second(self):
        """Variant overlaps only the second feature."""
        features = [[100, 200, 'exon'], [300, 400, 'CDS']]
        assert EnsemblDataService.check_overlap(350, 360, features) is True

    def test_multiple_features_no_overlap(self):
        """Variant falls between two features."""
        features = [[100, 200, 'exon'], [300, 400, 'CDS']]
        assert EnsemblDataService.check_overlap(210, 290, features) is False

    def test_default_features_info(self):
        """When features_info is None, uses the default [[0,1,'type']]."""
        # var_start=0, var_end=1 should overlap default [0,1,'type']
        assert EnsemblDataService.check_overlap(0, 1) is True


# ---------------------------------------------------------------------------
# get_altseq — plus-strand variants
# ---------------------------------------------------------------------------

class TestGetAltseqPlusStrand:
    """Tests for EnsemblDataService.get_altseq() on the plus strand."""

    def test_snp_substitution(self):
        """A single-base substitution on the plus strand modifies one base."""
        # A single-exon gene: positions 10..19, 10 bases
        ref_seq = Seq("ATGCATGCAT")  # 10 bases
        features_info = [[10, 19, 'exon']]  # exon from 10 to 19 (1-based genomic)
        # SNP at position 12: change the 3rd base (index 2)
        # ref_allele = 'G' (the base at position 12 in the exon), var_allele = 'T'
        ref_allele = Seq("G")
        var_allele = Seq("T")
        var_pos = 12

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info
        )

        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "ATTCATGCAT"
        assert len(coding_alt) == len(coding_ref)

    def test_insertion_on_plus_strand(self):
        """Insertion lengthens the resulting sequence."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        # Insertion at position 12: ref = 'G', var = 'GAA' (insert AA after G)
        ref_allele = Seq("G")
        var_allele = Seq("GAA")
        var_pos = 12

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info
        )

        assert str(coding_ref) == "ATGCATGCAT"
        # ref_allele[0] == var_allele[0] ('G' == 'G'), so the insertion logic is used
        # var_index_in_cds = 0 + (12 - 10) = 2; c = len(ref_allele) = 1
        # alt_seq = ref_seq[0:2] + 'GAA' + ref_seq[3:] = 'AT' + 'GAA' + 'CATGCAT' = 'ATGAACATGCAT'
        assert str(coding_alt) == "ATGAACATGCAT"
        assert len(coding_alt) == len(coding_ref) + 2

    def test_deletion_on_plus_strand(self):
        """Deletion shortens the resulting sequence."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        # Deletion at position 12: ref = 'GC', var = 'G' (delete C)
        ref_allele = Seq("GC")
        var_allele = Seq("G")
        var_pos = 12

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info
        )

        assert str(coding_ref) == "ATGCATGCAT"
        # ref_allele[0] == var_allele[0] ('G' == 'G'), insertion/deletion logic
        # var_index_in_cds = 0 + (12 - 10) = 2; c = len(ref_allele) = 2
        # alt_seq = ref_seq[0:2] + 'G' + ref_seq[4:] = 'AT' + 'G' + 'ATGCAT' = 'ATGATGCAT'
        assert str(coding_alt) == "ATGATGCAT"
        assert len(coding_alt) == len(coding_ref) - 1

    def test_variant_outside_features_returns_empty_alt(self):
        """If variant position is not inside any feature, alt_seq is empty."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        ref_allele = Seq("A")
        var_allele = Seq("T")
        var_pos = 5  # outside the exon 10..19

        _coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info
        )

        assert str(coding_alt) == ""

    def test_snp_with_cds_info(self):
        """SNP with explicit CDS coordinates to trim the coding region."""
        # ref_seq is the full transcript, CDS is a subset
        ref_seq = Seq("NNNNATGCATGCATNNN")  # 17 bases
        features_info = [[10, 23, 'exon']]
        cds_info = [5, 14]  # CDS from position 5 to 14 (1-based) -> index 4..14

        ref_allele = Seq("G")
        var_allele = Seq("T")
        var_pos = 12  # within the exon

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info, cds_info
        )

        # cds_info: start_coding_index = 5-1 = 4; stop_coding_index = 14
        # ref_seq[4:14] = 'ATGCATGCAT'
        assert str(coding_ref) == "ATGCATGCAT"


# ---------------------------------------------------------------------------
# get_altseq — minus-strand variants
# ---------------------------------------------------------------------------

class TestGetAltseqMinusStrand:
    """Tests for EnsemblDataService.get_altseq() on the minus strand.

    For minus-strand genes, the transcript sequence is stored in genomic order
    (5'→3' on the plus strand). The algorithm:
      1. Reverses ref_seq → puts it in minus-strand reading order
      2. Complements (not reverse-complements) alleles — the reversal on return
         completes the reverse-complement operation
      3. Applies variant using genomic coordinates on the reversed sequence
      4. Reverses both ref and alt back to genomic order on return

    This means: for a minus-strand gene at positions [10,19], genomic position 10
    (first in genomic order) is the LAST base of the transcript, so the variant
    should modify the last character of the returned sequence.
    """

    def test_snp_at_exon_start(self):
        """SNP at the first genomic position of a minus-strand exon.

        Position 10 is the first base genomically but the LAST base of the
        minus-strand transcript. The variant should modify the last character.
        """
        ref_seq = Seq("ATGCATGCAT")  # 10 bases
        features_info = [[10, 19, 'exon']]
        ref_allele = Seq("A")   # forward-strand allele at position 10
        var_allele = Seq("C")
        var_pos = 10

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '-', features_info
        )

        # Reversed ref: TACGTACGTA. complement(A)=T, complement(C)=G
        # var_index=0 on reversed → base T replaced by G → GACGTACGTA
        # Returned (reversed back): ref=ATGCATGCAT, alt=ATGCATGCAG
        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "ATGCATGCAG"
        assert len(coding_alt) == len(coding_ref)
        # Only the last base changed (position 10 = last in transcript)
        assert str(coding_ref)[:-1] == str(coding_alt)[:-1]
        assert str(coding_alt)[-1] == "G"

    def test_snp_at_exon_end(self):
        """SNP at the last genomic position of a minus-strand exon.

        Position 19 is the last base genomically but the FIRST base of the
        minus-strand transcript. The variant should modify the first character.
        """
        ref_seq = Seq("ATGCATGCAT")  # 10 bases
        features_info = [[10, 19, 'exon']]
        ref_allele = Seq("T")   # forward-strand base at position 19
        var_allele = Seq("C")
        var_pos = 19

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '-', features_info
        )

        # Reversed ref: TACGTACGTA. complement(T)=A, complement(C)=G
        # var_index = 0 + (19-10) = 9. Base at index 9 of reversed = A.
        # complement(T)=A matches. Replace with G → TACGTACGTG
        # Returned (reversed back): ref=ATGCATGCAT, alt=GTGCATGCAT
        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "GTGCATGCAT"
        assert len(coding_alt) == len(coding_ref)
        # Only the first base changed (position 19 = first in transcript)
        assert str(coding_ref)[1:] == str(coding_alt)[1:]

    def test_snp_at_middle_position(self):
        """SNP at a middle genomic position on minus strand."""
        ref_seq = Seq("ATGCATGCAT")  # 10 bases
        features_info = [[10, 19, 'exon']]
        ref_allele = Seq("A")
        var_allele = Seq("G")
        var_pos = 14  # middle of exon

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '-', features_info
        )

        # Reversed ref: TACGTACGTA. complement(A)=T, complement(G)=C
        # var_index = 0 + (14-10) = 4. Base at index 4 = T. Matches complement(A)=T.
        # Replace with C → TACGCACGTA. Reversed back → ATGCACGCAT
        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "ATGCACGCAT"
        # Position 14 maps to transcript index (19-14)=5 from the start
        assert str(coding_ref)[:5] == str(coding_alt)[:5]
        assert str(coding_alt)[5] == "C"
        assert str(coding_ref)[6:] == str(coding_alt)[6:]

    def test_insertion_on_minus_strand(self):
        """Insertion on minus strand: alleles are complemented, sequence lengthens."""
        ref_seq = Seq("ATGCATGCAT")  # 10 bases
        features_info = [[10, 19, 'exon']]
        # Insertion at position 12: ref='G', var='GAA' (forward strand)
        ref_allele = Seq("G")
        var_allele = Seq("GAA")
        var_pos = 12

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '-', features_info
        )

        assert len(coding_alt) == len(coding_ref) + 2
        assert str(coding_ref) == "ATGCATGCAT"

    def test_minus_strand_with_cds_info(self):
        """SNP on minus strand with explicit CDS coordinates."""
        ref_seq = Seq("NNNNATGCATGCATNNN")  # 17 bases
        features_info = [[10, 26, 'exon']]
        cds_info = [5, 14]  # CDS from position 5 to 14

        ref_allele = Seq("A")
        var_allele = Seq("G")
        var_pos = 14  # within the exon

        coding_ref, coding_alt = EnsemblDataService.get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '-', features_info, cds_info
        )

        # CDS: start=4 (0-based), stop=14. n=17.
        # Minus strand: reversed ref, then ref[17-14:17-4] = ref[3:13]
        assert len(coding_ref) == 10
        assert len(coding_alt) == 10
        assert str(coding_alt) != str(coding_ref)


# ---------------------------------------------------------------------------
# get_orfs_vcf
# ---------------------------------------------------------------------------

class TestGetOrfsVcf:
    """Tests for EnsemblDataService.get_orfs_vcf()."""

    def test_single_orf_translation(self):
        """With num_orfs=1, should translate from frame 0."""
        ref_seq = Seq("ATGAAATTT")  # M K F
        alt_seq = Seq("ATGAAAGTT")  # M K V
        ref_orfs, alt_orfs = EnsemblDataService.get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=1)
        assert len(ref_orfs) == 1
        assert len(alt_orfs) == 1
        assert str(ref_orfs[0]) == "MKF"
        assert str(alt_orfs[0]) == "MKV"

    def test_three_orfs(self):
        """With num_orfs=3, should return translations starting at offsets 0, 1, 2."""
        ref_seq = Seq("ATGAAATTTCCC")
        alt_seq = Seq("ATGAAATTTCCC")
        ref_orfs, alt_orfs = EnsemblDataService.get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=3)
        assert len(ref_orfs) == 3
        assert len(alt_orfs) == 3
        # Frame 0: ATGAAATTTCCC -> MKF(P)
        assert str(ref_orfs[0]) == "MKFP"

    def test_ref_and_alt_differ(self):
        """Ref and alt orfs differ when sequences differ."""
        ref_seq = Seq("ATGTTT")  # MF
        alt_seq = Seq("ATGGTT")  # MV
        ref_orfs, alt_orfs = EnsemblDataService.get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=1)
        assert str(ref_orfs[0]) != str(alt_orfs[0])
        assert str(ref_orfs[0]) == "MF"
        assert str(alt_orfs[0]) == "MV"

    def test_stop_codon_in_translation(self):
        """Stop codons are included in the translation (no to_stop)."""
        ref_seq = Seq("ATGTAATTT")  # M * F  (TAA is stop)
        alt_seq = Seq("ATGTAATTT")
        ref_orfs, _alt_orfs = EnsemblDataService.get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=1)
        assert "*" in str(ref_orfs[0])
