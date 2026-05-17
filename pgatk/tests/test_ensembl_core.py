"""Unit tests for EnsemblDataService core static methods."""

import logging
from pathlib import Path

import pytest
from Bio.Seq import Seq

from pgatk.ensembl.ensembl import EnsemblDataService, _split_vcf_into_batches

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_VCF_HEADER = (
    "##fileformat=VCFv4.1\n"
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
    'Format: ALLELE|CONSEQUENCE|FEATURE_TYPE|FEATURE|AMINO_ACIDS|SIFT">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _write_vcf(path, n_variants, header=_VCF_HEADER, chrom="1"):
    """Write a minimal VCF with n_variants data lines; return path as str."""
    with open(path, "w") as fh:
        fh.write(header)
        for i in range(n_variants):
            fh.write(f"{chrom}\t{100 + i}\t.\tA\tT\t.\t.\t"
                     f"CSQ=A|missense_variant|Transcript|ENST{i:06d}||\n")
    return str(path)


def _data_line_count(vcf_path):
    with open(vcf_path) as fh:
        return sum(1 for ln in fh if ln.strip() and not ln.startswith("#"))


def _header_lines(vcf_path):
    with open(vcf_path) as fh:
        return [ln for ln in fh if ln.startswith("#")]


class _StubSvc:
    """Duck-typed stub for _parse_annotation_indices unit tests.

    Provides only the instance attributes and get_logger() that the method
    reads, without going through ParameterConfiguration.__init__.
    """
    _annotation_field_name = "CSQ"
    _transcript_str = "FEATURE"
    _consequence_str = "CONSEQUENCE"
    _biotype_str = "transcript_biotype"

    @staticmethod
    def get_logger():
        return logging.getLogger("test_ensembl")


# Pipeline args that mirror the vcf-to-proteindb testdata invocation.
_CHUNK_PIPELINE_ARGS = {
    "annotation_field_name": "CSQ",
    "transcript_str": "FEATURE",
    "consequence_str": "CONSEQUENCE",
    "biotype_str": "feature_type",
    "include_biotypes": "mRNA,ncRNA",
    "af_field": "MAF",
    "protein_prefix": "ensvar",
}


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


# ---------------------------------------------------------------------------
# _split_vcf_into_batches
# ---------------------------------------------------------------------------

class TestSplitVcfIntoBatches:

    def test_header_only_returns_no_batches(self, tmp_path):
        """A VCF with no data lines produces an empty batch list."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=0)
        assert _split_vcf_into_batches(vcf, str(tmp_path), batch_size=10) == []

    def test_fewer_variants_than_batch_size(self, tmp_path):
        """When n_variants < batch_size all variants land in a single file."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=3)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=10)
        assert len(batches) == 1
        assert _data_line_count(batches[0]) == 3

    def test_exactly_batch_size_produces_one_file(self, tmp_path):
        """Exactly batch_size variants → exactly one batch file."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=5)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=5)
        assert len(batches) == 1
        assert _data_line_count(batches[0]) == 5

    def test_one_over_batch_size_produces_two_files(self, tmp_path):
        """batch_size + 1 variants must split into exactly 2 files (full + remainder)."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=6)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=5)
        assert len(batches) == 2
        assert _data_line_count(batches[0]) == 5
        assert _data_line_count(batches[1]) == 1

    def test_header_copied_to_every_batch(self, tmp_path):
        """Every batch file starts with the same header lines as the source VCF."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=12)
        expected = _header_lines(vcf)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=5)
        assert len(batches) > 1
        for batch in batches:
            assert _header_lines(batch) == expected

    def test_batch_files_are_ordered(self, tmp_path):
        """Returned list is in creation order (batch_0000 before batch_0001, …)."""
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=10)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=4)
        names = [Path(b).name for b in batches]
        assert names == sorted(names)

    def test_blank_lines_not_counted_as_variants(self, tmp_path):
        """Blank lines within the data section are skipped and not split into batches."""
        vcf_path = tmp_path / "blank.vcf"
        with vcf_path.open("w") as fh:
            fh.write(_VCF_HEADER)
            fh.write("1\t100\t.\tA\tT\t.\t.\t.\n")
            fh.write("\n")
            fh.write("1\t200\t.\tG\tC\t.\t.\t.\n")
        batches = _split_vcf_into_batches(str(vcf_path), str(tmp_path), batch_size=10)
        assert len(batches) == 1
        assert _data_line_count(batches[0]) == 2

    def test_total_variant_count_preserved_across_batches(self, tmp_path):
        """Sum of data lines across all batches equals the original variant count."""
        n = 23
        vcf = _write_vcf(tmp_path / "in.vcf", n_variants=n)
        batches = _split_vcf_into_batches(vcf, str(tmp_path), batch_size=7)
        assert sum(_data_line_count(b) for b in batches) == n

    def test_batches_may_span_chromosomes(self, tmp_path):
        """Batches are not bounded by chromosome; a single batch can contain
        variants from multiple chromosomes."""
        vcf_path = tmp_path / "multi_chrom.vcf"
        with vcf_path.open("w") as fh:
            fh.write(_VCF_HEADER)
            for chrom, pos in [("1", 100), ("1", 200), ("2", 100), ("2", 200)]:
                fh.write(f"{chrom}\t{pos}\t.\tA\tT\t.\t.\t.\n")
        batches = _split_vcf_into_batches(str(vcf_path), str(tmp_path), batch_size=10)
        assert len(batches) == 1
        assert _data_line_count(batches[0]) == 4


# ---------------------------------------------------------------------------
# _parse_annotation_indices
# ---------------------------------------------------------------------------

class TestParseAnnotationIndices:
    """Tests for EnsemblDataService._parse_annotation_indices().

    Uses _StubSvc (duck-typed) so the method can be tested in isolation
    without the ParameterConfiguration initialisation overhead.
    The testdata CSQ FORMAT is:
        ALLELE|CONSEQUENCE|FEATURE_TYPE|FEATURE|AMINO_ACIDS|SIFT
    After upper-casing: indices 0..5.  FEATURE → 3, CONSEQUENCE → 1.
    """

    # Metadata matching testdata/test.vcf CSQ FORMAT
    _META = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
        'Format: ALLELE|CONSEQUENCE|FEATURE_TYPE|FEATURE|AMINO_ACIDS|SIFT">\n',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ]

    def test_finds_consequence_and_feature_indices(self):
        """CONSEQUENCE at column 1 and FEATURE at column 3 are found correctly."""
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(_StubSvc(), self._META)
        assert ci == 1
        assert ti == 3

    def test_missing_biotype_returns_none(self):
        """transcript_biotype is absent from this FORMAT → bi is None."""
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(_StubSvc(), self._META)
        assert bi is None

    def test_biotype_found_when_present(self):
        """When biotype_str IS in the FORMAT, bi is set correctly."""
        meta = [
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
            'ALLELE|CONSEQUENCE|FEATURE_TYPE|FEATURE|TRANSCRIPT_BIOTYPE">\n',
            "#CHROM\n",
        ]
        stub = _StubSvc()
        stub._biotype_str = "transcript_biotype"
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(stub, meta)
        assert bi == 4

    def test_trailing_vcf_syntax_stripped_from_last_field(self):
        """The closing `">` of the INFO Description must not corrupt the last field.

        If TRANSCRIPT_BIOTYPE is the last FORMAT field, the raw split produces
        'TRANSCRIPT_BIOTYPE">\\n'; the strip must clean that so the index is found.
        """
        meta = [
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '
            'ALLELE|CONSEQUENCE|FEATURE|TRANSCRIPT_BIOTYPE">\n',
            "#CHROM\n",
        ]
        stub = _StubSvc()
        stub._biotype_str = "transcript_biotype"
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(stub, meta)
        assert bi is not None, "TRANSCRIPT_BIOTYPE should be found as the last field after stripping"
        assert bi == 3

    def test_transcript_absent_falls_back_to_zero(self):
        """When _transcript_str is not in the FORMAT, ti falls back to 0."""
        stub = _StubSvc()
        stub._transcript_str = "NOT_PRESENT"
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(stub, self._META)
        assert ti == 0

    def test_no_matching_info_header_all_indices_default(self):
        """When no ##INFO=<ID=CSQ,...> line exists, annotation_cols is empty;
        ti falls back to 0, ci and bi remain None."""
        meta = ["##fileformat=VCFv4.1\n", "#CHROM\tPOS\n"]
        ti, ci, bi = EnsemblDataService._parse_annotation_indices(_StubSvc(), meta)
        assert ti == 0
        assert ci is None
        assert bi is None

    def test_real_testdata_vcf_header(self, testdata_dir):
        """Parsed against the actual testdata/test.vcf header."""
        svc = EnsemblDataService({}, {})
        meta, _ = svc.vcf_from_file(str(testdata_dir / "test.vcf"))
        ti, ci, bi = svc._parse_annotation_indices(meta, str(testdata_dir / "test.vcf"))
        # CSQ FORMAT: ALLELE|CONSEQUENCE|FEATURE_TYPE|FEATURE|AMINO_ACIDS|SIFT
        assert ti == 3   # FEATURE
        assert ci == 1   # CONSEQUENCE
        assert bi is None  # transcript_biotype not in this VCF's CSQ FORMAT


# ---------------------------------------------------------------------------
# _vcf_to_proteindb_chunk
# ---------------------------------------------------------------------------

_STATS_KEYS = frozenset({
    "# variants with invalid record",
    "# variants not passing Filter",
    "# variants not passing AF threshold",
    "# feature IDs from VCF that are not found in the given FASTA file",
    "# variants successfully translated",
})


@pytest.fixture
def chunk_svc():
    """EnsemblDataService configured to match the testdata VCF."""
    return EnsemblDataService({}, _CHUNK_PIPELINE_ARGS)


def _run_chunk(svc, tmp_path, suffix="out", **kwargs):
    """Call _vcf_to_proteindb_chunk on testdata files; returns (path, stats)."""
    return svc._vcf_to_proteindb_chunk(
        "testdata/test.vcf",
        "testdata/test.fa",
        "testdata/test.gtf",
        str(tmp_path / f"{suffix}.fa"),
        **kwargs,
    )


class TestVcfToProteindbChunk:

    def test_return_type_is_tuple_of_str_and_dict(self, chunk_svc, tmp_path):
        """_vcf_to_proteindb_chunk returns (str, dict)."""
        result = _run_chunk(chunk_svc, tmp_path)
        assert isinstance(result, tuple) and len(result) == 2
        out_path, stats = result
        assert isinstance(out_path, str)
        assert isinstance(stats, dict)

    def test_stats_dict_has_all_expected_keys(self, chunk_svc, tmp_path):
        """The returned stats dict contains exactly the standard counter keys."""
        _, stats = _run_chunk(chunk_svc, tmp_path)
        assert set(stats.keys()) == _STATS_KEYS

    def test_all_stats_counts_are_non_negative(self, chunk_svc, tmp_path):
        """Every counter in the stats dict is ≥ 0."""
        _, stats = _run_chunk(chunk_svc, tmp_path)
        for key, val in stats.items():
            assert val >= 0, f"Negative count for '{key}': {val}"

    def test_log_summary_false_suppresses_info_line(self, chunk_svc, tmp_path, caplog):
        """With log_summary=False no 'Translation summary' INFO line is emitted."""
        with caplog.at_level(logging.INFO):
            _run_chunk(chunk_svc, tmp_path, log_summary=False)
        hits = [r for r in caplog.records
                if r.levelno == logging.INFO and "Translation summary" in r.message]
        assert hits == []

    def test_log_summary_true_emits_info_line(self, chunk_svc, tmp_path, caplog):
        """With log_summary=True (the default) a 'Translation summary' INFO line appears."""
        with caplog.at_level(logging.INFO):
            _run_chunk(chunk_svc, tmp_path, log_summary=True)
        hits = [r for r in caplog.records
                if r.levelno == logging.INFO and "Translation summary" in r.message]
        assert len(hits) >= 1

    def test_annotation_indices_precomputed_matches_auto(self, chunk_svc, tmp_path):
        """Pre-supplying annotation_indices produces identical stats to letting
        the chunk compute them from the VCF header itself."""
        meta, _ = chunk_svc.vcf_from_file("testdata/test.vcf")
        indices = chunk_svc._parse_annotation_indices(meta, "testdata/test.vcf")

        _, stats_auto = _run_chunk(chunk_svc, tmp_path, suffix="auto",
                                   log_summary=False)
        _, stats_pre = _run_chunk(chunk_svc, tmp_path, suffix="pre",
                                  annotation_indices=indices, log_summary=False)
        assert stats_auto == stats_pre

    def test_output_fasta_written_to_given_path(self, chunk_svc, tmp_path):
        """The returned path matches the argument and the file exists on disk."""
        out_path, _ = _run_chunk(chunk_svc, tmp_path)
        assert out_path == str(tmp_path / "out.fa")
        assert Path(out_path).exists()


# ---------------------------------------------------------------------------
# vcf_from_file — bgz/gzip-compressed input
# Regression test: previously raised UnicodeDecodeError on .gz/.bgz files
# because the file was opened with plain open() instead of gzip.open().
# ---------------------------------------------------------------------------

class TestBgzVcfReading:
    """vcf_from_file must transparently decompress .gz/.bgz VCF files."""

    _BGZ_VCF = "testdata/test_gnomad_bgz.vcf.gz"

    def test_vcf_from_file_reads_bgz_without_error(self):
        """Reading a gzip-compressed VCF must not raise UnicodeDecodeError."""
        meta, df = EnsemblDataService.vcf_from_file(self._BGZ_VCF)
        assert isinstance(meta, list)
        assert not df.empty

    def test_vcf_from_file_bgz_metadata_contains_vep_header(self):
        """The vep ##INFO header line must be present in the parsed metadata."""
        meta, _ = EnsemblDataService.vcf_from_file(self._BGZ_VCF)
        assert any("ID=vep" in line for line in meta)

    def test_vcf_from_file_bgz_data_has_expected_variant_count(self):
        """All four data lines in the compressed VCF are loaded into the DataFrame."""
        _, df = EnsemblDataService.vcf_from_file(self._BGZ_VCF)
        assert len(df) == 4

    def test_vcf_from_file_bgz_chrom_column_parsed(self):
        """The CHROM column is correctly parsed from the compressed file."""
        _, df = EnsemblDataService.vcf_from_file(self._BGZ_VCF)
        assert list(df["CHROM"].unique()) == ["22"]


# ---------------------------------------------------------------------------
# _vcf_to_proteindb_chunk — gnomAD bgz-style VCF with GENCODE transcript FASTA
#
# testdata/test_gnomad_gencode.fa uses GENCODE v44-style pipe-delimited headers
# with versioned transcript IDs (e.g. ENST00000643195.1).  The gnomAD VCF
# carries bare IDs (ENST00000643195); the pipeline version-strips FASTA keys
# on the first cache miss so both forms resolve to the same sequence.
# ---------------------------------------------------------------------------

_GNOMAD_CHUNK_ARGS = {
    "annotation_field_name": "vep",
    "transcript_str": "Feature",        # column name in VEP Format header
    "consequence_str": "Consequence",
    "biotype_str": "biotype",           # matches BIOTYPE column in VEP header
    "include_biotypes": "protein_coding",
    "include_consequences": "missense_variant,inframe_deletion",
    "af_field": "AF_afr",
    "af_threshold": 0.01,
    "protein_prefix": "gnomadvar",
}


@pytest.fixture
def gnomad_svc():
    """EnsemblDataService configured for gnomAD vep + AF_afr pipeline."""
    return EnsemblDataService({}, _GNOMAD_CHUNK_ARGS)


def _run_gnomad_chunk(svc, tmp_path, suffix="out", **kwargs):
    """Call _vcf_to_proteindb_chunk with the gnomAD bgz VCF and GENCODE v44 references."""
    return svc._vcf_to_proteindb_chunk(
        "testdata/test_gnomad_bgz.vcf.gz",
        "testdata/test_gnomad_gencode.fa",
        "testdata/test_gnomad_gencode.gtf",
        str(tmp_path / f"{suffix}.fa"),
        **kwargs,
    )


class TestVcfToProteindbChunkGnomad:
    """End-to-end pipeline tests using a gnomAD-style gzip-compressed VCF.

    Reference files are real GENCODE v44 records extracted from
    gencode.v44.transcripts.fa and gencode.v44.annotation.gtf:
      - test_gnomad_gencode.fa  -- pipe-delimited headers, versioned IDs
                                   (ENST00000643195.1, ENST00000550946.5);
                                   no CDS= token → pipeline uses exon features
                                   and 3-frame translation for all transcripts
      - test_gnomad_gencode.gtf -- GENCODE format with chr22 chromosome prefix
                                   and transcript_type attribute (vs Ensembl's
                                   transcript_biotype); chr prefix is stripped
                                   on comparison so no mismatch with VCF chrom

    Test VCF contains four variants for ENST00000643195 (OR11H1, chr22):
      - missense_variant  AF_afr=0.02  → passes AF + consequence filters → translated
      - inframe_deletion  AF_afr=0.015 → passes AF + consequence filters → translated
      - missense_variant  AF_afr=0.005 → below AF threshold (0.01)       → filtered
      - stop_gained       AF_afr=0.03  → fails consequence filter         → skipped silently
    """

    def test_gnomad_bgz_chunk_runs_without_error(self, gnomad_svc, tmp_path):
        """Pipeline completes without raising on a gzip-compressed gnomAD VCF."""
        out_path, stats = _run_gnomad_chunk(gnomad_svc, tmp_path)
        assert isinstance(out_path, str)
        assert isinstance(stats, dict)

    def test_gnomad_bgz_chunk_stats_has_all_expected_keys(self, gnomad_svc, tmp_path):
        """Stats dict contains the standard counter keys."""
        _, stats = _run_gnomad_chunk(gnomad_svc, tmp_path)
        assert set(stats.keys()) == _STATS_KEYS

    def test_gnomad_bgz_chunk_af_filter_applied(self, gnomad_svc, tmp_path):
        """The variant with AF_afr=0.005 (below threshold 0.01) is counted as filtered."""
        _, stats = _run_gnomad_chunk(gnomad_svc, tmp_path)
        assert stats["# variants not passing AF threshold"] == 1

    def test_gnomad_bgz_chunk_translates_passing_variants(self, gnomad_svc, tmp_path):
        """At least one variant passes all filters and is successfully translated."""
        _, stats = _run_gnomad_chunk(gnomad_svc, tmp_path)
        assert stats["# variants successfully translated"] >= 1

    def test_gnomad_bgz_chunk_output_fasta_exists(self, gnomad_svc, tmp_path):
        """The output protein FASTA file is created at the given path."""
        out_path, _ = _run_gnomad_chunk(gnomad_svc, tmp_path)
        assert Path(out_path).exists()

    def test_gnomad_bgz_chunk_annotation_indices_resolve_correctly(self, gnomad_svc):
        """vep FORMAT header is parsed: Feature→transcript_index, Consequence→consequence_index."""
        meta, _ = gnomad_svc.vcf_from_file("testdata/test_gnomad_bgz.vcf.gz")
        ti, ci, bi = gnomad_svc._parse_annotation_indices(meta, "testdata/test_gnomad_bgz.vcf.gz")
        assert ti == 6   # Feature is the 7th column (0-indexed)
        assert ci == 1   # Consequence is the 2nd column
        assert bi == 7   # BIOTYPE is the 8th column

    def test_gencode_versioned_id_key_extraction(self):
        """get_key strips at the first '|', returning the versioned ID from a GENCODE header."""
        header = "ENST00000643195.1|ENSG00000130538.6|OTTHUMG00000183615.4|OR11H1-201|OR11H1|948|protein_coding|"
        assert EnsemblDataService.get_key(header) == "ENST00000643195.1"

    def test_gencode_versioned_id_resolves_to_sequence(self):
        """Bare VCF transcript ID (no version) maps to the versioned GENCODE FASTA entry.

        The pipeline builds {stripped_id: versioned_id} on the first cache miss.
        This test verifies that ENST00000643195 (VCF) → ENST00000643195.1 (FASTA).
        """
        from Bio import SeqIO
        td = SeqIO.to_dict(
            SeqIO.parse("testdata/test_gnomad_gencode.fa", "fasta"),
            key_function=lambda r: EnsemblDataService.get_key(r.description),
        )
        # Keys in the dict are versioned
        assert "ENST00000643195.1" in td
        assert "ENST00000550946.5" in td
        # Version-strip mapping matches bare IDs from the VCF
        mapping = {k.split(".")[0]: k for k in td.keys()}
        assert mapping["ENST00000643195"] == "ENST00000643195.1"
        assert mapping["ENST00000550946"] == "ENST00000550946.5"
