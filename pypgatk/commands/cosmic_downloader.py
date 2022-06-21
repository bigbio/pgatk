import logging

import click

from pypgatk.cgenomes.cosmic_downloader import CosmicDownloadService
import pkgutil

from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file

log = logging.getLogger(__name__)
try:
    default_config_text = pkgutil.get_data(__name__, "../config/cosmic_config.yaml").decode()
except Exception:
    try:
        default_config_text = pkgutil.get_data(__name__, "config/cosmic_config.yaml").decode()
    except ValueError:
        log.info("The config file is now available !! ")

@click.command('cosmic-downloader', short_help='Command to download the cosmic mutation database')
@click.option('-c', '--config_file', help='Configuration file for the ensembl data downloader pipeline')
@click.option('-o', '--output_directory', help='Output directory for the peptide databases')
@click.option('-u', '--username',
              help="Username for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.option('-p', '--password',
              help="Password for cosmic database -- please if you don't have one register here (https://cancer.sanger.ac.uk/cosmic/register)")
@click.option("--url_file", help='Add the url to a downloaded file')
@click.pass_context
def cosmic_downloader(ctx, config_file, output_directory, username, password, url_file):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("cosmic_config.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  pipeline_arguments = {}
  if output_directory is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
  if username is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_USER] = username

  if password is not None:
    pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD] = password

  cosmic_downloader_service = CosmicDownloadService(config_data, pipeline_arguments)
  cosmic_downloader_service.download_mutation_file(url_file_name = url_file)
