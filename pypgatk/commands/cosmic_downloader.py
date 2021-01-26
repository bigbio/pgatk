import logging
import os

import click

from pypgatk.cgenomes.cosmic_downloader import CosmicDownloadService
from pypgatk.toolbox.exceptions import AppConfigException

this_dir, this_filename = os.path.split(__file__)


@click.command('cosmic-downloader', short_help='Command to download the cosmic mutation database')
@click.option('-c', '--config_file',
              help='Configuration file for the ensembl data downloader pipeline',
              default=this_dir + '/../config/cosmic_config.yaml')
@click.option('-o', '--output_directory',
              help='Output directory for the peptide databases',
              default="./database_cosmic/")
@click.option('-u', '--username',
              help="Username for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.option('-p', '--password',
              help="Password for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.pass_context
def cosmic_downloader(ctx, config_file, output_directory, username, password):
  if config_file is None:
    msg = "The config file for the pipeline is missing, please provide one "
    logging.error(msg)
    raise AppConfigException(msg)

  pipeline_arguments = {}
  if output_directory is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
  if username is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_USER] = username

  if password is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD] = password

  cosmic_downloader_service = CosmicDownloadService(config_file, pipeline_arguments)

  cosmic_downloader_service.download_mutation_file()
