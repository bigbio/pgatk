"""Unit tests for pgatk.toolbox.vcf_utils — shared VCF-to-protein helpers."""

import io

from Bio.Seq import Seq

from pgatk.toolbox.vcf_utils import check_overlap, get_altseq, get_orfs_vcf, write_output


# ---------------------------------------------------------------------------
# check_overlap
# ---------------------------------------------------------------------------

class TestCheckOverlap:
    """Tests for check_overlap()."""

    def test_variant_fully_within_feature(self):
        """Variant is fully contained inside a single feature."""
        features = [[100, 200, 'exon']]
        assert check_overlap(120, 130, features) is True

    def test_variant_outside_feature(self):
        """Variant does not overlap any feature."""
        features = [[100, 200, 'exon']]
        assert check_overlap(50, 60, features) is False

    def test_start_minus_one_always_true(self):
        """When var_start is -1, overlap is unconditionally True."""
        features = [[100, 200, 'exon']]
        assert check_overlap(-1, 10, features) is True

    def test_partial_overlap_start(self):
        """Variant overlaps the start boundary of a feature."""
        features = [[100, 200, 'exon']]
        assert check_overlap(90, 110, features) is True

    def test_partial_overlap_end(self):
        """Variant overlaps the end boundary of a feature."""
        features = [[100, 200, 'exon']]
        assert check_overlap(190, 210, features) is True

    def test_variant_covers_entire_feature(self):
        """Variant fully covers the feature."""
        features = [[100, 200, 'exon']]
        assert check_overlap(50, 300, features) is True

    def test_default_features_info(self):
        """Default features_info [[0,1,'type']] is used when None."""
        assert check_overlap(0, 1) is True

    def test_multiple_features_overlap_second(self):
        """Variant overlaps only the second feature."""
        features = [[100, 200, 'exon'], [300, 400, 'CDS']]
        assert check_overlap(350, 360, features) is True

    def test_multiple_features_no_overlap(self):
        """Variant falls between two features — no overlap."""
        features = [[100, 200, 'exon'], [300, 400, 'CDS']]
        assert check_overlap(210, 290, features) is False


# ---------------------------------------------------------------------------
# get_altseq
# ---------------------------------------------------------------------------

class TestGetAltseq:
    """Tests for get_altseq()."""

    def test_snp_plus_strand(self):
        """Single-base substitution on plus strand returns correct alt."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        ref_allele = Seq("G")
        var_allele = Seq("T")
        var_pos = 12

        coding_ref, coding_alt = get_altseq(
            ref_seq, ref_allele, var_allele, var_pos, '+', features_info
        )

        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "ATTCATGCAT"

    def test_returns_tuple_of_length_two(self):
        """get_altseq always returns a 2-tuple."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        result = get_altseq(
            ref_seq, Seq("G"), Seq("T"), 12, '+', features_info
        )
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_variant_outside_features_returns_empty_alt(self):
        """Variant position outside all features yields empty alt_seq."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        _ref, alt = get_altseq(
            ref_seq, Seq("A"), Seq("T"), 5, '+', features_info
        )
        assert str(alt) == ""

    def test_snp_minus_strand(self):
        """SNP on minus strand modifies the expected position."""
        ref_seq = Seq("ATGCATGCAT")
        features_info = [[10, 19, 'exon']]
        coding_ref, coding_alt = get_altseq(
            ref_seq, Seq("A"), Seq("C"), 10, '-', features_info
        )
        assert str(coding_ref) == "ATGCATGCAT"
        assert str(coding_alt) == "ATGCATGCAG"


# ---------------------------------------------------------------------------
# get_orfs_vcf
# ---------------------------------------------------------------------------

class TestGetOrfsVcf:
    """Tests for get_orfs_vcf()."""

    def test_single_orf(self):
        """With num_orfs=1, produces one ORF per sequence."""
        ref_seq = Seq("ATGAAATTT")  # M K F
        alt_seq = Seq("ATGAAAGTT")  # M K V
        ref_orfs, alt_orfs = get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=1)
        assert len(ref_orfs) == 1
        assert len(alt_orfs) == 1
        assert str(ref_orfs[0]) == "MKF"
        assert str(alt_orfs[0]) == "MKV"

    def test_three_orfs(self):
        """With num_orfs=3, returns translations for three reading frames."""
        ref_seq = Seq("ATGAAATTTCCC")
        alt_seq = Seq("ATGAAATTTCCC")
        ref_orfs, alt_orfs = get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=3)
        assert len(ref_orfs) == 3
        assert len(alt_orfs) == 3
        # Frame 0: ATGAAATTTCCC -> MKFP
        assert str(ref_orfs[0]) == "MKFP"

    def test_ref_and_alt_differ(self):
        """Different input sequences produce different ORFs."""
        ref_seq = Seq("ATGTTT")  # MF
        alt_seq = Seq("ATGGTT")  # MV
        ref_orfs, alt_orfs = get_orfs_vcf(ref_seq, alt_seq, 1, num_orfs=1)
        assert str(ref_orfs[0]) == "MF"
        assert str(alt_orfs[0]) == "MV"


# ---------------------------------------------------------------------------
# write_output
# ---------------------------------------------------------------------------

class TestWriteOutput:
    """Tests for write_output()."""

    def test_writes_fasta_format(self):
        """Single ORF is written in FASTA format with correct header."""
        buf = io.StringIO()
        write_output(seq_id="prot1", desc="some desc", seqs=["MKFV"], prots_fn=buf)
        content = buf.getvalue()
        assert content == ">prot1 some desc\nMKFV\n"

    def test_multiple_orfs_get_numbered(self):
        """Multiple ORFs receive _1, _2, ... suffixes."""
        buf = io.StringIO()
        write_output(seq_id="prot1", desc="", seqs=["MKF", "KFP", "FPS"], prots_fn=buf)
        lines = buf.getvalue().strip().split('\n')
        assert lines[0] == ">prot1_1"
        assert lines[2] == ">prot1_2"
        assert lines[4] == ">prot1_3"

    def test_seqs_filter_skips_matching_orfs(self):
        """ORFs present in seqs_filter are not written."""
        buf = io.StringIO()
        write_output(
            seq_id="prot1", desc="", seqs=["MKF", "AAA"],
            prots_fn=buf, seqs_filter=["MKF"]
        )
        content = buf.getvalue()
        assert "MKF" not in content
        assert "AAA" in content

    def test_empty_desc_no_trailing_space(self):
        """When desc is empty, no trailing space appears after the ID."""
        buf = io.StringIO()
        write_output(seq_id="prot1", desc="", seqs=["MKFV"], prots_fn=buf)
        content = buf.getvalue()
        assert content == ">prot1\nMKFV\n"
