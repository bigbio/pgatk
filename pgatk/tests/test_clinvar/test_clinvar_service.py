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


# ---------------------------------------------------------------------------
# TestMolecularConsequenceParser
# ---------------------------------------------------------------------------


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
