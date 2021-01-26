import logging
import os

import click

from pypgatk.cgenomes.cbioportal_downloader import CbioPortalDownloadService
from pypgatk.toolbox.exceptions import AppConfigException

this_dir, this_filename = os.path.split(__file__)


@click.command('cbioportal-downloader', short_help=' Command to download the the cbioportal studies')
@click.option('-c', '--config_file', help='Configuration file for the ensembl data downloader pipeline',
              default=this_dir + '/../config/cbioportal_config.yaml')
@click.option('-o', '--output_directory', help='Output directory for the peptide databases',
              default="./database_cbioportal/")
@click.option('-l', '--list_studies',
              help='Print the list of all the studies in cBioPortal (https://www.cbioportal.org)', is_flag=True)
@click.option('-d', '--download_study',
              help="Download a specific Study from cBioPortal -- (all to download all studies)")
@click.option('-th', '--multithreading', help=' Enable multithreading to download multiple files ad the same time',
              is_flag=True)
@click.pass_context
def cbioportal_downloader(ctx, config_file, output_directory, list_studies, download_study, multithreading):
  if config_file is None:
    msg = "The config file for the pipeline is missing, please provide one "
    logging.error(msg)
    raise AppConfigException(msg)

  pipeline_arguments = {}
  if output_directory is not None:
    pipeline_arguments[CbioPortalDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
  if list_studies:
    pipeline_arguments[CbioPortalDownloadService.CONFIG_LIST_STUDIES] = list_studies
  if multithreading is not None:
    pipeline_arguments[CbioPortalDownloadService.CONFIG_MULTITHREADING] = multithreading

  cbioportal_downloader_service = CbioPortalDownloadService(config_file, pipeline_arguments)

  if list_studies:
    list_of_studies = cbioportal_downloader_service.get_cancer_studies()
    print(list_of_studies)

  if download_study is not None:
    cbioportal_downloader_service.download_study(download_study)
