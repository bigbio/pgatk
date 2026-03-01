import logging

import click

from pgatk.commands.utils import print_help
from pgatk.ensembl.ensembl import EnsemblDataService
from pgatk.config.registry import load_config

log = logging.getLogger(__name__)


@click.command('threeframe-translation', short_help="Command to perform 3'frame translation")
@click.option('-c', '--config_file',
              help='Configuration to perform conversion between ENSEMBL Files')
@click.option('-in', '--input_fasta', help='input_fasta file to perform the translation')
@click.option('-t', '--translation_table', help='Translation table default value 1')
@click.option('-out', '--output', help='Output File')
@click.pass_context
def threeframe_translation(ctx, config_file, input_fasta, translation_table, output):

    config_data = load_config("ensembl_config", config_file)

    if input_fasta is None:
        print_help()

    pipeline_arguments = {}

    if translation_table is not None:
        pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table

    if output is not None:
        pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.three_frame_translation(input_fasta)
