"""End-to-end integration test for the ClinVar pipeline."""

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from pgatk.cli import cli


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"


@pytest.fixture()
def clinvar_testdata(tmp_path):
    """Copy ClinVar test data to tmp_path so gffutils .db files are isolated.

    gffutils creates a ``<name>.db`` file next to the GTF.  By copying the GTF
    into *tmp_path* we avoid polluting the source tree and ensure each test run
    starts from a clean state.
    """
    for name in (
        "mini_clinvar.vcf",
        "mini_refseq.gtf",
        "mini_refseq_protein.faa",
        "mini_assembly_report.txt",
    ):
        shutil.copy2(TESTDATA_DIR / name, tmp_path / name)
    return tmp_path


class TestClinVarCLIIntegration:
    """Test the clinvar-to-proteindb CLI command end-to-end."""

    def test_cli_produces_output(self, tmp_path, clinvar_testdata):
        output_file = tmp_path / "clinvar_output.fa"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "clinvar-to-proteindb",
            "--vcf", str(clinvar_testdata / "mini_clinvar.vcf"),
            "--gtf", str(clinvar_testdata / "mini_refseq.gtf"),
            "--fasta", str(clinvar_testdata / "mini_refseq_protein.faa"),
            "--assembly-report", str(clinvar_testdata / "mini_assembly_report.txt"),
            "--output", str(output_file),
        ])
        assert result.exit_code == 0, (
            f"CLI failed with exit code {result.exit_code}:\n"
            f"{result.output}\n{result.exception}"
        )
        assert output_file.exists()

    def test_cli_help_shows_options(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["clinvar-to-proteindb", "--help"])
        assert result.exit_code == 0
        assert "clinvar" in result.output.lower()
        assert "--vcf" in result.output
        assert "--assembly-report" in result.output

    def test_ncbi_downloader_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["ncbi-downloader", "--help"])
        assert result.exit_code == 0
        assert "--output-dir" in result.output
        assert "--force" in result.output

    def test_cli_filters_synonymous_variants(self, tmp_path, clinvar_testdata):
        """Synonymous variants should be filtered by MC even if Pathogenic."""
        output_file = tmp_path / "clinvar_output.fa"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "clinvar-to-proteindb",
            "--vcf", str(clinvar_testdata / "mini_clinvar.vcf"),
            "--gtf", str(clinvar_testdata / "mini_refseq.gtf"),
            "--fasta", str(clinvar_testdata / "mini_refseq_protein.faa"),
            "--assembly-report", str(clinvar_testdata / "mini_assembly_report.txt"),
            "--output", str(output_file),
        ])
        assert result.exit_code == 0
        content = output_file.read_text()
        # rs00004 is Pathogenic but synonymous — should NOT appear
        assert "rs00004" not in content
        # rs00001 is Pathogenic missense — SHOULD appear
        assert "rs00001" in content

    def test_pipeline_output_contains_clnsig_in_header(self, tmp_path, clinvar_testdata):
        """Output FASTA headers should contain CLNSIG and gene symbol."""
        output_file = tmp_path / "clinvar_output.fa"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "clinvar-to-proteindb",
            "--vcf", str(clinvar_testdata / "mini_clinvar.vcf"),
            "--gtf", str(clinvar_testdata / "mini_refseq.gtf"),
            "--fasta", str(clinvar_testdata / "mini_refseq_protein.faa"),
            "--assembly-report", str(clinvar_testdata / "mini_assembly_report.txt"),
            "--output", str(output_file),
        ])
        assert result.exit_code == 0
        content = output_file.read_text()
        assert "Pathogenic" in content
