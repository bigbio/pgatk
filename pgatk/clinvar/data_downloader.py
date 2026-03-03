"""Download NCBI RefSeq and ClinVar reference files.

Uses urllib for FTP/HTTPS downloads. Skips existing files unless force=True.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path

from pgatk.toolbox.general import download_file, check_create_folders

logger = logging.getLogger(__name__)

_DEFAULT_REFSEQ_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/"
    "GRCh38_latest/refseq_identifiers/"
)
_DEFAULT_CLINVAR_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
)

_REFSEQ_FILES = [
    "GRCh38_latest_genomic.gtf.gz",
    "GRCh38_latest_protein.faa.gz",
    "GRCh38_latest_assembly_report.txt",
]
_CLINVAR_FILES = [
    "clinvar.vcf.gz",
    "clinvar.vcf.gz.tbi",
]


class NcbiDataDownloader:
    """Download NCBI RefSeq and ClinVar reference data."""

    def __init__(
        self,
        output_dir: str,
        refseq_base_url: str = _DEFAULT_REFSEQ_BASE,
        clinvar_base_url: str = _DEFAULT_CLINVAR_BASE,
    ) -> None:
        self._output_dir = output_dir
        self._refseq_base = refseq_base_url
        self._clinvar_base = clinvar_base_url

    def ensure_output_dir(self) -> None:
        """Create output directory if it doesn't exist."""
        check_create_folders([self._output_dir])

    def get_refseq_urls(self) -> list[str]:
        """Return list of RefSeq file URLs to download."""
        return [self._refseq_base + f for f in _REFSEQ_FILES]

    def get_clinvar_urls(self) -> list[str]:
        """Return list of ClinVar file URLs to download."""
        return [self._clinvar_base + f for f in _CLINVAR_FILES]

    def expected_files(self) -> list[str]:
        """Return list of expected local file paths after download."""
        all_files = _REFSEQ_FILES + _CLINVAR_FILES
        return [os.path.join(self._output_dir, f) for f in all_files]

    def download_all(self, force: bool = False) -> list[str]:
        """Download all required files. Returns list of downloaded file paths.

        :param force: If True, re-download even if files exist.
        """
        self.ensure_output_dir()
        downloaded = []

        all_urls = self.get_refseq_urls() + self.get_clinvar_urls()
        all_names = _REFSEQ_FILES + _CLINVAR_FILES

        for url, name in zip(all_urls, all_names):
            local_path = os.path.join(self._output_dir, name)

            if os.path.exists(local_path) and not force:
                logger.info("File already exists, skipping: %s", local_path)
                downloaded.append(local_path)
                continue

            logger.info("Downloading %s -> %s", url, local_path)
            result = download_file(url, local_path, logger)
            if result:
                downloaded.append(result)
            else:
                logger.error("Failed to download: %s", url)

        return downloaded
