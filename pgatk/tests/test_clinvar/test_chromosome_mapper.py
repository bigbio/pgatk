"""Tests for pgatk.clinvar.chromosome_mapper."""

from pathlib import Path

import pytest

from pgatk.clinvar.chromosome_mapper import ChromosomeMapper


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"
ASSEMBLY_REPORT = TESTDATA_DIR / "mini_assembly_report.txt"


class TestChromosomeMapper:
    """Tests for ChromosomeMapper."""

    @pytest.fixture
    def mapper(self):
        return ChromosomeMapper.from_assembly_report(str(ASSEMBLY_REPORT))

    def test_numeric_to_refseq(self, mapper):
        assert mapper.map_chrom("1", "refseq") == "NC_000001.11"

    def test_refseq_to_numeric(self, mapper):
        assert mapper.map_chrom("NC_000001.11", "numeric") == "1"

    def test_numeric_to_ucsc(self, mapper):
        assert mapper.map_chrom("1", "ucsc") == "chr1"

    def test_ucsc_to_numeric(self, mapper):
        assert mapper.map_chrom("chr1", "numeric") == "1"

    def test_refseq_to_ucsc(self, mapper):
        assert mapper.map_chrom("NC_000001.11", "ucsc") == "chr1"

    def test_x_chromosome(self, mapper):
        assert mapper.map_chrom("X", "refseq") == "NC_000023.11"
        assert mapper.map_chrom("NC_000023.11", "numeric") == "X"

    def test_mt_chromosome(self, mapper):
        assert mapper.map_chrom("MT", "refseq") == "NC_012920.1"
        assert mapper.map_chrom("NC_012920.1", "ucsc") == "chrM"

    def test_unknown_chromosome_returns_input(self, mapper):
        assert mapper.map_chrom("UNKNOWN_CONTIG", "refseq") == "UNKNOWN_CONTIG"

    def test_invalid_convention_raises(self, mapper):
        with pytest.raises(ValueError, match="convention"):
            mapper.map_chrom("1", "invalid_convention")

    def test_chr_prefixed_input(self, mapper):
        """Input with chr prefix should be recognized as UCSC convention."""
        assert mapper.map_chrom("chrX", "numeric") == "X"
        assert mapper.map_chrom("chrX", "refseq") == "NC_000023.11"
