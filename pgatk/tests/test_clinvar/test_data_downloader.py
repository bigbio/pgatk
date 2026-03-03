"""Tests for NCBI data downloader (no actual downloads)."""

from pathlib import Path

import pytest

from pgatk.clinvar.data_downloader import NcbiDataDownloader


class TestNcbiDataDownloader:
    """Tests for NcbiDataDownloader."""

    def test_build_refseq_urls(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        urls = downloader.get_refseq_urls()
        assert len(urls) == 3
        assert any("GRCh38_latest_genomic.gtf.gz" in u for u in urls)
        assert any("GRCh38_latest_rna.fna.gz" in u for u in urls)
        assert any("assembly_report.txt" in u for u in urls)

    def test_build_clinvar_urls(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        urls = downloader.get_clinvar_urls()
        assert len(urls) == 2
        assert any("clinvar.vcf.gz" in u for u in urls)

    def test_expected_files_list(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        files = downloader.expected_files()
        assert len(files) == 5  # 3 refseq + 2 clinvar

    def test_output_dir_created(self, tmp_path):
        out_dir = tmp_path / "ncbi_data"
        downloader = NcbiDataDownloader(output_dir=str(out_dir))
        downloader.ensure_output_dir()
        assert out_dir.exists()

    def test_custom_base_urls(self):
        downloader = NcbiDataDownloader(
            output_dir="/tmp/test",
            refseq_base_url="https://custom.example.com/refseq/",
            clinvar_base_url="https://custom.example.com/clinvar/",
        )
        urls = downloader.get_refseq_urls()
        assert all(u.startswith("https://custom.example.com/refseq/") for u in urls)
        clinvar_urls = downloader.get_clinvar_urls()
        assert all(u.startswith("https://custom.example.com/clinvar/") for u in clinvar_urls)

    def test_expected_files_include_output_dir(self):
        downloader = NcbiDataDownloader(output_dir="/my/output")
        files = downloader.expected_files()
        assert all(f.startswith("/my/output/") for f in files)
