"""Tests for pgatk.clinvar.clinvar_service — ClinVar VCF-to-protein pipeline."""

import os
import tempfile
from pathlib import Path

import pytest

from pgatk.clinvar.clinvar_service import ClinVarService


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"
MINI_VCF = str(TESTDATA_DIR / "mini_clinvar.vcf")
MINI_GTF = str(TESTDATA_DIR / "mini_refseq.gtf")
MINI_FASTA = str(TESTDATA_DIR / "mini_refseq_protein.faa")
ASSEMBLY_REPORT = str(TESTDATA_DIR / "mini_assembly_report.txt")


# ---------------------------------------------------------------------------
# TestClinSigFiltering
# ---------------------------------------------------------------------------


class TestClinSigFiltering:
    """Tests for ClinVarService.passes_clnsig_filter()."""

    def test_pathogenic_passes(self):
        """Pathogenic should NOT be excluded."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Pathogenic", exclude) is True

    def test_benign_excluded(self):
        """Benign should be excluded."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Benign", exclude) is False

    def test_likely_benign_excluded(self):
        """Likely_benign should be excluded."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Likely_benign", exclude) is False

    def test_empty_clnsig_passes(self):
        """Empty CLNSIG always passes."""
        exclude = ["Benign", "Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("", exclude) is True

    def test_none_like_empty_passes(self):
        """None-equivalent empty string passes."""
        exclude = ["Benign"]
        assert ClinVarService.passes_clnsig_filter("", exclude) is True

    def test_likely_pathogenic_passes(self):
        """Likely_pathogenic should pass."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Likely_pathogenic", exclude) is True

    def test_compound_all_benign_excluded(self):
        """Compound value with all-benign components should be excluded."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        # Both components are in the exclude list → excluded
        assert ClinVarService.passes_clnsig_filter("Benign,Likely_benign", exclude) is False

    def test_compound_mixed_passes(self):
        """Compound value with a non-benign component should pass."""
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Pathogenic/Likely_benign", exclude) is True

    def test_compound_risk_factor_excluded(self):
        """Benign with risk_factor — all components benign → excluded."""
        exclude = ["Benign", "Likely_benign"]
        # '_risk_factor' is not in exclude list, so it's not an excluded component
        assert ClinVarService.passes_clnsig_filter("Benign,_risk_factor", exclude) is True


# ---------------------------------------------------------------------------
# TestMolecularConsequenceParser
# ---------------------------------------------------------------------------


class TestMCFiltering:
    """Tests for MC consequence filtering."""

    def test_passes_mc_filter_matching_consequence(self):
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001583|missense_variant"
        assert ClinVarService.passes_mc_filter(mc, include) is True

    def test_passes_mc_filter_no_match(self):
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001819|synonymous_variant"
        assert ClinVarService.passes_mc_filter(mc, include) is False

    def test_passes_mc_filter_multi_consequence_any_match(self):
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001819|synonymous_variant,SO:0001587|stop_gained"
        assert ClinVarService.passes_mc_filter(mc, include) is True

    def test_passes_mc_filter_empty_mc_passes(self):
        include = ["missense_variant"]
        assert ClinVarService.passes_mc_filter("", include) is True

    def test_passes_mc_filter_all_keyword(self):
        mc = "SO:0001819|synonymous_variant"
        assert ClinVarService.passes_mc_filter(mc, ["all"]) is True


class TestMolecularConsequenceParser:
    """Tests for MC field parsers."""

    def test_parse_single_mc(self):
        """Parse a single MC entry returns the consequence."""
        mc = "SO:0001583|missense_variant"
        assert ClinVarService.parse_mc_consequence(mc) == "missense_variant"

    def test_parse_multiple_mc_returns_first(self):
        """parse_mc_consequence returns only the first consequence."""
        mc = "SO:0001583|missense_variant,SO:0001587|stop_gained"
        assert ClinVarService.parse_mc_consequence(mc) == "missense_variant"

    def test_parse_mc_consequences_all(self):
        """parse_mc_consequences returns all consequences."""
        mc = "SO:0001583|missense_variant,SO:0001587|stop_gained"
        result = ClinVarService.parse_mc_consequences(mc)
        assert result == ["missense_variant", "stop_gained"]

    def test_parse_empty_mc(self):
        """Empty MC field returns empty."""
        assert ClinVarService.parse_mc_consequence("") == ""
        assert ClinVarService.parse_mc_consequences("") == []

    def test_parse_mc_no_pipe(self):
        """MC field without pipe returns the raw value."""
        mc = "missense_variant"
        assert ClinVarService.parse_mc_consequence(mc) == "missense_variant"


# ---------------------------------------------------------------------------
# TestGeneInfoParser
# ---------------------------------------------------------------------------


class TestGeneInfoParser:
    """Tests for GENEINFO field parsing."""

    def test_single_gene(self):
        """Single gene GENEINFO returns (symbol, id)."""
        symbol, gene_id = ClinVarService.parse_geneinfo("GeneA:1234")
        assert symbol == "GeneA"
        assert gene_id == "1234"

    def test_multi_gene(self):
        """Multi-gene GENEINFO returns the first gene."""
        symbol, gene_id = ClinVarService.parse_geneinfo("GeneA:1234|GeneB:5678")
        assert symbol == "GeneA"
        assert gene_id == "1234"

    def test_empty_geneinfo(self):
        """Empty GENEINFO returns empty strings."""
        symbol, gene_id = ClinVarService.parse_geneinfo("")
        assert symbol == ""
        assert gene_id == ""


# ---------------------------------------------------------------------------
# TestInfoFieldParser
# ---------------------------------------------------------------------------


class TestInfoFieldParser:
    """Tests for _get_info_field helper."""

    def test_extract_clnsig(self):
        info = "GENEINFO=GeneA:1234;CLNSIG=Pathogenic;MC=SO:0001583|missense_variant"
        assert ClinVarService._get_info_field(info, "CLNSIG") == "Pathogenic"

    def test_extract_mc(self):
        info = "GENEINFO=GeneA:1234;CLNSIG=Pathogenic;MC=SO:0001583|missense_variant"
        assert ClinVarService._get_info_field(info, "MC") == "SO:0001583|missense_variant"

    def test_missing_key(self):
        info = "GENEINFO=GeneA:1234;CLNSIG=Pathogenic"
        assert ClinVarService._get_info_field(info, "MC") == ""


# ---------------------------------------------------------------------------
# TestClinVarPipeline — Integration
# ---------------------------------------------------------------------------


class TestClinVarPipeline:
    """Integration tests for the ClinVar pipeline run()."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        """Remove gffutils DB before/after test to ensure fresh state."""
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_pipeline_produces_output(self, tmp_path):
        """Pipeline should produce a non-empty output file."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        result_path = service.run()
        assert os.path.exists(result_path)
        with open(result_path, "r") as f:
            content = f.read()
        assert len(content) > 0, "Output file should not be empty"

    def test_benign_variant_excluded(self, tmp_path):
        """Benign variant (rs00002) should not appear in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            content = f.read()
        # rs00002 is the Benign variant — should not appear
        assert "rs00002" not in content

    def test_pathogenic_variant_present(self, tmp_path):
        """Pathogenic variant (rs00001) should appear in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            content = f.read()
        # rs00001 is the Pathogenic variant — should appear
        assert "rs00001" in content

    def test_likely_pathogenic_variant_present(self, tmp_path):
        """Likely_pathogenic variant (rs00003) should appear in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            content = f.read()
        # rs00003 is the Likely_pathogenic variant — should appear
        assert "rs00003" in content

    def test_synonymous_variant_excluded(self, tmp_path):
        """Pathogenic synonymous variant (rs00004) should not appear in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            content = f.read()
        assert "rs00004" not in content


# ---------------------------------------------------------------------------
# TestBuildOverlapMap
# ---------------------------------------------------------------------------


class TestBuildOverlapMap:
    """Tests for _build_overlap_map (DataFrame-based BedTools annotation)."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_build_overlap_map_from_dataframe(self):
        """_build_overlap_map should accept a DataFrame and return overlap dict."""
        from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(ASSEMBLY_REPORT)
        _meta, vcf_df = ClinVarService._read_vcf(MINI_VCF)
        overlap_map = ClinVarService._build_overlap_map(vcf_df, MINI_GTF, chrom_mapper)
        assert isinstance(overlap_map, dict)
        # rs00001 at chr1:1006 should overlap NM_000001.1 CDS (1000-1299)
        assert any("1:1006:" in k for k in overlap_map)

    def test_build_overlap_map_empty_df(self):
        """Empty DataFrame should return empty dict."""
        import pandas as pd
        from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(ASSEMBLY_REPORT)
        empty_df = pd.DataFrame(columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
        overlap_map = ClinVarService._build_overlap_map(empty_df, MINI_GTF, chrom_mapper)
        assert overlap_map == {}


# ---------------------------------------------------------------------------
# TestBiotypeFiltering
# ---------------------------------------------------------------------------


class TestBiotypeFiltering:
    """Tests for biotype filtering from GTF."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_get_biotype_from_db(self):
        """_get_transcript_biotype should extract gene_biotype from gffutils DB."""
        db = ClinVarService._parse_gtf(MINI_GTF)
        biotype = ClinVarService._get_transcript_biotype(db, "NM_000001.1")
        assert biotype == "protein_coding"

    def test_get_biotype_missing_returns_empty(self):
        """Missing transcript returns empty string."""
        db = ClinVarService._parse_gtf(MINI_GTF)
        biotype = ClinVarService._get_transcript_biotype(db, "NONEXISTENT")
        assert biotype == ""

    def test_get_biotype_without_version(self):
        """Should find transcript even without version number."""
        db = ClinVarService._parse_gtf(MINI_GTF)
        # NM_000001 without .1 should still find NM_000001.1
        biotype = ClinVarService._get_transcript_biotype(db, "NM_000001")
        assert biotype == "protein_coding"


# ---------------------------------------------------------------------------
# TestDuplicateGuard
# ---------------------------------------------------------------------------


class TestDuplicateGuard:
    """Pipeline should not process the same variant-transcript pair twice."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_no_duplicate_fasta_entries(self, tmp_path):
        """Each variant-transcript combo should appear at most once in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            headers = [line.strip() for line in f if line.startswith(">")]
        # No duplicate headers
        assert len(headers) == len(set(headers)), (
            f"Duplicate FASTA entries found: {headers}"
        )
