import errno
import glob
import logging
import os
import subprocess  # nosec B404 - referenced only to catch CalledProcessError from gffread
import sys

import click

from pgatk.ensembl.data_downloader import EnsemblDataDownloadService
from pgatk.config.registry import load_config

log = logging.getLogger(__name__)

@click.command('ensembl-downloader', short_help='Command to download the ensembl information')
@click.option('-c', '--config_file', help='Configuration file for the ensembl data downloader pipeline')
@click.option('-o', '--output_directory',
              help='Output directory for the peptide databases')
@click.option('-t', '--taxonomy',
              help='Taxonomy identifiers (comma separated list can be given) that will be use to download the data from Ensembl')
@click.option('-fp', '--folder_prefix_release', help='Output folder prefix to download the data')
@click.option('-sg', '--skip_gtf', help="Skip the GTF file during the download", is_flag=True)
@click.option('-sp', '--skip_protein', help="Skip the protein fasta file during download", is_flag=True)
@click.option('-sc', '--skip_cds', help='Skip the CDS file download', is_flag=True)
@click.option('-sdn', '--skip_cdna', help='Skip the cDNA file download', is_flag=True)
@click.option('-sn', '--skip_ncrna', help='Skip the ncRNA file download', is_flag=True)
@click.option('-sd', '--skip_dna', help='Skip the DNA (reference genome assembly) file download', is_flag=True)
@click.option('-sv', '--skip_vcf', help='Skip the VCF variant file', is_flag=True)
@click.option('-en', '--ensembl_name',
              help='Ensembl name code to download, it can be use instead of taxonomy (e.g. homo_sapiens)')
@click.option('--grch37', is_flag=True, help='Download a previous version GRCh37 of ensembl genomes')
@click.option('--url_file', help='Add the url to a downloaded file')
@click.option(
    '--generate-transcripts', 'generate_transcripts',
    is_flag=True, default=False,
    help=(
        'Run gffread after download to generate transcripts.fa with CDS= headers. '
        'Requires the genome assembly (do not use --skip_dna) and gffread in PATH '
        '(conda install -c bioconda gffread). Has no effect when --skip_gtf is set.'
    ),
)
@click.pass_context
def ensembl_downloader(ctx, config_file, output_directory, taxonomy, folder_prefix_release, skip_gtf, skip_protein,
                       skip_cds, skip_cdna, skip_ncrna, skip_dna, skip_vcf,
                       ensembl_name, grch37, url_file, generate_transcripts):
    """ This tool enables to download from enseml ftp the FASTA and GTF files"""

    config_data = load_config("ensembl_downloader", config_file)

    if taxonomy is None and ensembl_name is None:
        raise click.UsageError("Either --taxonomy or --ensembl_name is required.")

    pipeline_arguments = {}

    if output_directory is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory

    if taxonomy is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_TAXONOMY] = taxonomy

    if ensembl_name is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_ENSEMBL_NAME] = ensembl_name

    if skip_protein is not None and skip_protein:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_PROTEIN] = skip_protein

    if skip_gtf is not None and skip_gtf:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_GTF] = skip_gtf

    if skip_cds is not None and skip_cds:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_CDS] = skip_cds

    if skip_ncrna is not None and skip_ncrna:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_NCRNA] = skip_ncrna

    if skip_dna is not None and skip_dna:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_DNA] = skip_dna

    if skip_cdna is not None and skip_cdna:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_CDNA] = skip_cdna

    if skip_vcf is not None and skip_vcf:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_SKIP_VCF] = skip_vcf

    if grch37 is not None and grch37:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_KEY_GRCh37] = grch37

    if folder_prefix_release is not None:
        pipeline_arguments[EnsemblDataDownloadService.CONFIG_PREFIX_RELEASE_FOLDER] = folder_prefix_release

    ensembl_download_service = EnsemblDataDownloadService(config_data, pipeline_arguments)

    if taxonomy is not None:
        if not ensembl_download_service.validate_taxonomies():
            sys.exit(errno.EFAULT)

    ensembl_download_service.download_database_by_species(url_file)

    if generate_transcripts:
        if skip_gtf:
            click.echo(
                "Warning: --generate-transcripts has no effect when --skip_gtf is set.",
                err=True,
            )
        else:
            out_dir = ensembl_download_service.get_local_path_root_ensembl_repo()
            gtf_files = sorted(glob.glob(os.path.join(out_dir, "*.gtf.gz")))
            genome_files = sorted(glob.glob(os.path.join(out_dir, "*.dna_sm.toplevel.fa.gz")))

            if not genome_files:
                click.echo(
                    f"Error: no genome FASTA (*.dna_sm.toplevel.fa.gz) found in {out_dir}. "
                    "Re-run without --skip_dna to download the genome assembly.",
                    err=True,
                )
                raise SystemExit(1)

            if not gtf_files:
                click.echo(f"Error: no GTF file found in {out_dir}.", err=True)
                raise SystemExit(1)

            # Map species+assembly prefix (e.g. "Homo_sapiens.GRCh38") to genome path.
            genome_map = {}
            for gf in genome_files:
                prefix = ".".join(os.path.basename(gf).split(".")[:2])
                genome_map[prefix] = gf

            for gtf_path in gtf_files:
                prefix = ".".join(os.path.basename(gtf_path).split(".")[:2])
                genome_fna = genome_map.get(prefix) or next(iter(genome_map.values()))
                output_fasta = os.path.join(
                    out_dir,
                    "transcripts.fa" if len(gtf_files) == 1 else f"{prefix}_transcripts.fa",
                )
                try:
                    click.echo(f"Running gffread to generate {os.path.basename(output_fasta)} ...")
                    EnsemblDataDownloadService.generate_transcripts(genome_fna, gtf_path, output_fasta)
                    click.echo(f"Transcripts written to {output_fasta}")
                except (FileNotFoundError, subprocess.CalledProcessError) as exc:
                    click.echo(f"Error: {exc}", err=True)
                    raise SystemExit(1)
