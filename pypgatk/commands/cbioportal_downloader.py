import logging

import click

from pypgatk.cgenomes.cbioportal_downloader import CbioPortalDownloadService
from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


@click.command('cbioportal-downloader', short_help='Command to download the the cbioportal studies')
@click.option('-c', '--config_file', help='Configuration file for the ensembl data downloader pipeline')
@click.option('-o', '--output_directory', help='Output directory for the peptide databases')
@click.option('-l', '--list_studies',
              help='Print the list of all the studies in cBioPortal (https://www.cbioportal.org)', is_flag=True)
@click.option('-d', '--download_study',
              help="Download a specific Study from cBioPortal -- (all to download all studies)")
@click.option('-th', '--multithreading', help=' Enable multithreading to download multiple files ad the same time',
              is_flag=True)
@click.option('--url_file', help='Add the url to a downloaded file')
@click.pass_context
def cbioportal_downloader(ctx, config_file, output_directory, list_studies, download_study, multithreading, url_file):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    pipeline_arguments = {}

    if output_directory is not None:
        pipeline_arguments[CbioPortalDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
    if list_studies:
        pipeline_arguments[CbioPortalDownloadService.CONFIG_LIST_STUDIES] = list_studies
    if multithreading is not None:
        pipeline_arguments[CbioPortalDownloadService.CONFIG_MULTITHREADING] = multithreading

    cbioportal_downloader_service = CbioPortalDownloadService(config_data, pipeline_arguments)

    if list_studies:
        list_of_studies = cbioportal_downloader_service.get_cancer_studies()
        print(list_of_studies)

    if download_study is not None:
        cbioportal_downloader_service.download_study(download_study, url_file_name=url_file)
