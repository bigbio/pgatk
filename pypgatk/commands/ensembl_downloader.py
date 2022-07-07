import errno
import logging
import sys

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.data_downloader import EnsemblDataDownloadService
from pypgatk.toolbox.general import read_yaml_from_file

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
@click.pass_context
def ensembl_downloader(ctx, config_file, output_directory, taxonomy, folder_prefix_release, skip_gtf, skip_protein,
                       skip_cds, skip_cdna, skip_ncrna, skip_dna, skip_vcf,
                       ensembl_name, grch37, url_file):
    """ This tool enables to download from enseml ftp the FASTA and GTF files"""

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if taxonomy is None and ensembl_name is None:
        print_help()

    # Parse pipelines parameters.
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
