"""Download NCBI RefSeq and ClinVar reference files.

Uses urllib for FTP/HTTPS downloads. Skips existing files unless force=True.
"""
from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess  # nosec B404 - used only with hardcoded gffread argv lists, never shell=True

from pgatk.toolbox.general import download_file, check_create_folders


def _resolve_genome_fna(genome_fna: str) -> tuple[str, str | None]:
    """Return (path_for_gffread, temp_path_to_clean_up).

    gffread requires random-access into the genome FASTA.  Regular .gz does not
    support seeks, so if the caller passes a .gz path we reuse an already-
    decompressed sibling file or decompress to one on the fly.
    Decompression is written to a sibling ``.tmp`` file and atomically renamed
    so an interrupted run cannot leave behind a partial / corrupt FASTA that
    later invocations would silently reuse.
    Returns the second element only when *we* created the file (caller must remove it).
    """
    if not genome_fna.endswith(".gz"):
        return genome_fna, None
    plain = genome_fna[:-3]
    if os.path.exists(plain):
        return plain, None
    logger.info("Decompressing genome FASTA for gffread (random-access required): %s", genome_fna)
    tmp = plain + ".tmp"
    try:
        with gzip.open(genome_fna, "rb") as fi, open(tmp, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.replace(tmp, plain)
    except BaseException:
        # On any failure (including KeyboardInterrupt) remove the partial temp
        # file so the next run re-decompresses from scratch.
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        raise
    return plain, plain

logger = logging.getLogger(__name__)

_REFSEQ_BASE = {
    "GRCh38": (
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/"
        "GRCh38_latest/refseq_identifiers/"
    ),
    "GRCh37": (
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/"
        "GRCh37_latest/refseq_identifiers/"
    ),
}
_CLINVAR_BASE = {
    "GRCh38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/",
    "GRCh37": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/",
}

# Keep old names for backwards compatibility
_DEFAULT_REFSEQ_BASE = _REFSEQ_BASE["GRCh38"]
_DEFAULT_CLINVAR_BASE = _CLINVAR_BASE["GRCh38"]

_REFSEQ_FILES = {
    "GRCh38": [
        "GRCh38_latest_genomic.fna.gz",
        "GRCh38_latest_genomic.gff.gz",
        "GRCh38_latest_assembly_report.txt",
    ],
    "GRCh37": [
        "GRCh37_latest_genomic.fna.gz",
        "GRCh37_latest_genomic.gff.gz",
        "GRCh37_latest_assembly_report.txt",
    ],
}
_CLINVAR_FILES = [
    "clinvar.vcf.gz",
    "clinvar.vcf.gz.tbi",
]


def _require_gffread() -> None:
    if not shutil.which("gffread"):
        raise FileNotFoundError(
            "gffread not found in PATH. "
            "Install it with: conda install -c bioconda gffread"
        )


class NcbiDataDownloader:
    """Download NCBI RefSeq and ClinVar reference data."""

    def __init__(
        self,
        output_dir: str,
        refseq_base_url: str = _DEFAULT_REFSEQ_BASE,
        clinvar_base_url: str = _DEFAULT_CLINVAR_BASE,
        grch37: bool = False,
    ) -> None:
        self._output_dir = output_dir
        self._build = "GRCh37" if grch37 else "GRCh38"
        self._refseq_base = _REFSEQ_BASE[self._build] if grch37 else refseq_base_url
        self._clinvar_base = _CLINVAR_BASE[self._build] if grch37 else clinvar_base_url

    def ensure_output_dir(self) -> None:
        check_create_folders([self._output_dir])

    def get_refseq_urls(self) -> list[str]:
        return [self._refseq_base + f for f in _REFSEQ_FILES[self._build]]

    def get_clinvar_urls(self) -> list[str]:
        return [self._clinvar_base + f for f in _CLINVAR_FILES]

    def expected_files(self) -> list[str]:
        all_files = _REFSEQ_FILES[self._build] + _CLINVAR_FILES
        return [os.path.join(self._output_dir, f) for f in all_files]

    def download_all(self, force: bool = False) -> list[str]:
        """Download all required files. Returns list of downloaded file paths."""
        self.ensure_output_dir()
        downloaded = []

        all_urls = self.get_refseq_urls() + self.get_clinvar_urls()
        all_names = _REFSEQ_FILES[self._build] + _CLINVAR_FILES
        total = len(all_names)

        for i, (url, name) in enumerate(zip(all_urls, all_names), 1):
            local_path = os.path.join(self._output_dir, name)

            if not force and os.path.exists(local_path):
                logger.info("[%d/%d] Already exists, skipping: %s", i, total, os.path.basename(local_path))
                downloaded.append(local_path)
                continue

            logger.info("[%d/%d] Downloading: %s", i, total, name)
            result = download_file(url, local_path, logger, progress_prefix=f"  [{name}]")
            if result:
                logger.info("[%d/%d] Saved to: %s", i, total, result)
                downloaded.append(result)
            else:
                logger.error("[%d/%d] Failed to download: %s", i, total, url)

        return downloaded

    @staticmethod
    def generate_transcripts(genome_fna: str, gff_file: str, output_fasta: str) -> str:
        """Run gffread to produce a transcript FASTA with CDS= coordinate headers.

        Required for ``clinvar-to-proteindb``. The GFF3 (not the GTF) must be
        used because only GFF3 carries the explicit ID=/Parent= linkage that
        gffread needs to resolve transcript → gene hierarchy.
        """
        _require_gffread()
        genome_path, _tmp = _resolve_genome_fna(genome_fna)
        try:
            cmd = ["gffread", "-F", "-w", output_fasta, "-g", genome_path, gff_file]
            logger.info("Running: %s", " ".join(cmd))
            subprocess.run(cmd, check=True)  # nosec B603 - argv list with hardcoded "gffread"
        finally:
            if _tmp and os.path.exists(_tmp):
                os.remove(_tmp)
        logger.info("Transcripts written to %s", output_fasta)
        return output_fasta

    @staticmethod
    def generate_cds_fasta(genome_fna: str, gff_file: str, output_fasta: str) -> str:
        """Extract spliced CDS sequences into a plain FASTA for cbioportal-to-proteindb.

        Uses ``gffread -x`` to write only the coding sequence (ATG → stop codon)
        for every annotated mRNA. Record IDs in the output are RefSeq accessions
        (e.g. ``NM_001762.3``); the translation code strips the version suffix with
        ``split(".")[0]`` → ``NM_001762``, matching the Transcript_ID values
        returned by the cBioPortal REST API.

        Parameters
        ----------
        genome_fna : str
            Path to the decompressed genomic FASTA
            (``GRCh38_latest_genomic.fna`` or ``GRCh37_latest_genomic.fna``).
        gff_file : str
            Path to the decompressed GFF3 annotation file.
        output_fasta : str
            Output path for the CDS FASTA (e.g. ``refseq_cds.fa``).
        """
        _require_gffread()
        genome_path, _tmp = _resolve_genome_fna(genome_fna)
        try:
            cmd = ["gffread", "-x", output_fasta, "-g", genome_path, gff_file]
            logger.info("Running: %s", " ".join(cmd))
            subprocess.run(cmd, check=True)  # nosec B603 - argv list with hardcoded "gffread"
        finally:
            if _tmp and os.path.exists(_tmp):
                os.remove(_tmp)
        logger.info("CDS FASTA written to %s", output_fasta)
        return output_fasta
