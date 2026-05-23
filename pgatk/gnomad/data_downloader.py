"""Download GENCODE annotation files and gnomAD per-chromosome VCF files.

GencodeDownloader  — GENCODE GTF + primary-assembly genome FASTA for a given
                     release, plus an optional gffread step to generate a
                     transcript FASTA with CDS= headers.
GnomadVcfDownloader — gnomAD exomes VCF + TBI files, all chromosomes in parallel.
"""
from __future__ import annotations

import gzip
import logging
import os
import re
import shutil
import subprocess  # nosec B404 - used only with hardcoded gffread argv lists, never shell=True
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

from pgatk.toolbox.general import download_file, check_create_folders


def _resolve_genome_fna(genome_fna: str) -> tuple[str, str | None]:
    """Return (path_for_gffread, temp_path_to_clean_up).

    gffread requires random-access into the genome FASTA.  Regular .gz (deflate)
    does not support seeks, so if the caller passes a .gz path we either reuse
    an already-decompressed sibling file or decompress to one on the fly.
    Decompression is staged through a sibling ``.tmp`` file and atomically
    renamed so an interrupted run cannot leave behind a partial / corrupt
    FASTA that later invocations would silently reuse.
    Returns the second element as a path only when *we* created the file and
    the caller must delete it afterward.
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
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        raise
    return plain, plain

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GENCODE version compatibility check
# ---------------------------------------------------------------------------

def _extract_gencode_version_from_vcf(vcf_path: str) -> Optional[str]:
    """Scan VCF header lines for ``##gencode_version=`` and return the version number."""
    open_fn = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
    try:
        with open_fn(vcf_path, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                # e.g. ##gencode_version=GENCODE 44
                m = re.search(r"##gencode_version=.*?(\d+)", line, re.IGNORECASE)
                if m:
                    return m.group(1)
    except Exception as exc:
        logger.debug("Could not read VCF header for GENCODE version check: %s", exc)
    return None


def _extract_gencode_version_from_gtf(gtf_path: str) -> Optional[str]:
    """Scan GTF comment lines for the GENCODE release version number."""
    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    try:
        with open_fn(gtf_path, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                # e.g. ##description: ... version 44 (Ensembl 110)
                m = re.search(r"version\s+(\d+)", line, re.IGNORECASE)
                if m:
                    return m.group(1)
    except Exception as exc:
        logger.debug("Could not read GTF header for GENCODE version check: %s", exc)
    return None


def check_gencode_version_compatibility(vcf_path: str, gtf_path: str) -> None:
    """Warn if the GENCODE version recorded in the VCF header differs from the GTF.

    gnomAD VCF headers carry ``##gencode_version=GENCODE N`` reflecting the
    GENCODE release used for VEP annotation.  GENCODE GTF headers carry
    ``##description: ... version N ...``.  A mismatch means transcript
    coordinates in the GTF may not correspond to the VEP consequence calls in
    the VCF, which can produce incorrect translations.

    No warning is emitted when the version cannot be determined from either
    file (e.g. for non-gnomAD VCFs or hand-crafted test files).
    """
    vcf_ver = _extract_gencode_version_from_vcf(vcf_path)
    gtf_ver = _extract_gencode_version_from_gtf(gtf_path)
    if vcf_ver is None or gtf_ver is None:
        return
    if vcf_ver != gtf_ver:
        logger.warning(
            "GENCODE version mismatch: VCF was annotated with GENCODE v%s but "
            "the supplied GTF is GENCODE v%s. Transcript coordinates may not "
            "match VEP annotations, leading to incorrect translations. "
            "Re-run 'pgatk gencode-downloader --release %s' to download the "
            "matching GTF and genome FASTA.",
            vcf_ver, gtf_ver, vcf_ver,
        )
    else:
        logger.debug("GENCODE version check passed: VCF and GTF are both v%s", vcf_ver)

_GENCODE_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    "Gencode_human/release_{release}/"
)
_GNOMAD_BASE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/{version}/vcf/{dataset}/"
)

# Genome FASTA filename is stable across releases (no patch number in name)
_GENOME_FASTA_GZ = "GRCh38.primary_assembly.genome.fa.gz"
_GENOME_FASTA = "GRCh38.primary_assembly.genome.fa"

# Standard chromosome set used by gnomAD (UCSC-style names with "chr" prefix)
ALL_CHROMOSOMES: list[str] = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ---------------------------------------------------------------------------
# GENCODE
# ---------------------------------------------------------------------------

class GencodeDownloader:
    """Download GENCODE human GTF and primary-assembly genome FASTA.

    The pre-built ``gencode.v{release}.transcripts.fa.gz`` is intentionally
    NOT downloaded: its pipe-delimited headers have no ``CDS=`` token, so the
    pipeline would fall back to 3-frame exon translation.  Instead, the
    primary-assembly genomic FASTA is downloaded and ``generate_transcripts``
    runs ``gffread -F`` to produce a transcript FASTA with embedded
    ``CDS=start-end`` coordinates — enabling 1-frame CDS translation.

    Typical usage::

        downloader = GencodeDownloader(output_dir="gencode_data", release=44)
        downloader.download_all()
        downloader.generate_transcripts(
            genome_fna="gencode_data/GRCh38.primary_assembly.genome.fa",
            gtf_file="gencode_data/gencode.v44.annotation.gtf",
            output_fasta="gencode_data/transcripts.fa",
        )
    """

    def __init__(
        self,
        output_dir: str,
        release: int = 44,
        base_url: Optional[str] = None,
    ) -> None:
        self._output_dir = output_dir
        self._release = release
        self._base = (base_url or _GENCODE_BASE).format(release=release)

    def ensure_output_dir(self) -> None:
        check_create_folders([self._output_dir])

    def get_urls(self) -> list[str]:
        """Return download URLs for the GTF annotation and primary-assembly genome FASTA."""
        v = self._release
        return [
            self._base + f"gencode.v{v}.annotation.gtf.gz",
            self._base + _GENOME_FASTA_GZ,
        ]

    def expected_files(self) -> list[str]:
        """Return expected local paths after download (files remain compressed)."""
        v = self._release
        return [
            os.path.join(self._output_dir, f"gencode.v{v}.annotation.gtf.gz"),
            os.path.join(self._output_dir, _GENOME_FASTA_GZ),
        ]

    def download_all(self, force: bool = False) -> list[str]:
        """Download GTF and primary-assembly genome FASTA. Returns local file paths."""
        self.ensure_output_dir()
        v = self._release
        filenames = [
            f"gencode.v{v}.annotation.gtf.gz",
            _GENOME_FASTA_GZ,
        ]
        downloaded = []
        for url, name in zip(self.get_urls(), filenames):
            local_gz = os.path.join(self._output_dir, name)
            if not force and os.path.exists(local_gz):
                logger.info("Already exists, skipping: %s", local_gz)
                downloaded.append(local_gz)
                continue
            logger.info("Downloading %s -> %s", url, local_gz)
            result = download_file(url, local_gz, logger)
            if result:
                downloaded.append(result)
            else:
                logger.error("Failed to download: %s", url)
        return downloaded

    @staticmethod
    def generate_transcripts(genome_fna: str, gtf_file: str, output_fasta: str) -> str:
        """Run gffread -F to produce a transcript FASTA with CDS= coordinate headers.

        Parameters
        ----------
        genome_fna:
            Path to the decompressed primary-assembly FASTA
            (``GRCh38.primary_assembly.genome.fa``).
        gtf_file:
            Path to the GENCODE annotation GTF
            (``gencode.v{release}.annotation.gtf``).
        output_fasta:
            Output path for the transcript FASTA (e.g. ``transcripts.fa``).

        Returns
        -------
        str
            Path to the generated FASTA file.

        Raises
        ------
        FileNotFoundError
            If ``gffread`` is not found in PATH.
        subprocess.CalledProcessError
            If gffread exits with a non-zero status.
        """
        if not shutil.which("gffread"):
            raise FileNotFoundError(
                "gffread not found in PATH. "
                "Install it with: conda install -c bioconda gffread"
            )
        genome_for_gffread, _tmp = _resolve_genome_fna(genome_fna)
        try:
            cmd = ["gffread", "-F", "-w", output_fasta, "-g", genome_for_gffread, gtf_file]
            logger.info("Running: %s", " ".join(cmd))
            subprocess.run(cmd, check=True)  # nosec B603 - argv list with hardcoded "gffread"
        finally:
            if _tmp and os.path.exists(_tmp):
                os.remove(_tmp)
        logger.info("Transcripts written to %s", output_fasta)
        return output_fasta


# ---------------------------------------------------------------------------
# gnomAD VCF
# ---------------------------------------------------------------------------

class GnomadVcfDownloader:
    """Download gnomAD per-chromosome VCF (BGZF) and tabix index files in parallel.

    Supports both gnomAD genomes and exomes VCFs via the ``dataset`` parameter.
    Each chromosome yields two files: ``*.vcf.bgz`` and ``*.vcf.bgz.tbi``.
    Files are stored as-is — tabix requires them to stay BGZF-compressed.

    Both v3.x and v4.x ship genome VCFs; v4.x additionally ships exome VCFs.
    Check the gnomAD downloads page to confirm which dataset and GENCODE release
    your VCF files were annotated with before choosing ``--release`` for
    ``gencode-downloader``.
    """

    def __init__(
        self,
        output_dir: str,
        version: str = "3.1.2",
        dataset: str = "genomes",
        chromosomes: Optional[list[str]] = None,
        workers: int = 4,
        base_url: Optional[str] = None,
    ) -> None:
        self._output_dir = output_dir
        self._version = version
        self._dataset = dataset
        self._chromosomes = chromosomes if chromosomes is not None else ALL_CHROMOSOMES
        self._workers = workers
        self._base = (base_url or _GNOMAD_BASE).format(version=version, dataset=dataset)

    def ensure_output_dir(self) -> None:
        check_create_folders([self._output_dir])

    def get_urls(self) -> list[tuple[str, str]]:
        """Return ``(url, local_filename)`` pairs for every VCF and TBI file."""
        pairs = []
        v = self._version
        d = self._dataset
        for chrom in self._chromosomes:
            vcf_name = f"gnomad.{d}.v{v}.sites.{chrom}.vcf.bgz"
            tbi_name = vcf_name + ".tbi"
            pairs.append((self._base + vcf_name, vcf_name))
            pairs.append((self._base + tbi_name, tbi_name))
        return pairs

    def _download_one(self, url: str, name: str, force: bool) -> Optional[str]:
        local_path = os.path.join(self._output_dir, name)
        if not force and os.path.exists(local_path):
            logger.info("Already exists, skipping: %s", local_path)
            return local_path
        logger.info("Downloading %s -> %s", url, local_path)
        result = download_file(url, local_path, logger)
        if not result:
            logger.error("Failed to download: %s", url)
        return result

    def download_all(
        self,
        force: bool = False,
        progress_callback=None,
    ) -> list[str]:
        """Download all VCF + TBI files in parallel.

        Parameters
        ----------
        force:
            Re-download files that already exist.
        progress_callback:
            Optional callable invoked after each file completes:
            ``callback(completed: int, total: int, name: str)``.

        Returns
        -------
        list[str]
            Paths of successfully downloaded files.
        """
        self.ensure_output_dir()
        pairs = self.get_urls()
        total = len(pairs)
        results: list[str] = []
        completed = 0

        with ThreadPoolExecutor(max_workers=self._workers) as pool:
            futures = {
                pool.submit(self._download_one, url, name, force): name
                for url, name in pairs
            }
            for future in as_completed(futures):
                name = futures[future]
                completed += 1
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                except Exception as exc:
                    logger.error("Error downloading %s: %s", name, exc)
                if progress_callback is not None:
                    progress_callback(completed, total, name)

        return results
