import logging

import click

from pypgatk.cgenomes.cosmic_downloader import CosmicDownloadService

from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


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

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    pipeline_arguments = {}
    if output_directory is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY] = output_directory
    if username is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_USER] = username

    if password is not None:
        pipeline_arguments[CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD] = password

    cosmic_downloader_service = CosmicDownloadService(config_data, pipeline_arguments)
    cosmic_downloader_service.download_mutation_file(url_file_name=url_file)
