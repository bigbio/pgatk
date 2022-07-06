import logging

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService
from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


@click.command('threeframe-translation', short_help="Command to perform 3'frame translation")
@click.option('-c', '--config_file',
              help='Configuration to perform conversion between ENSEMBL Files')
@click.option('-in', '--input_fasta', help='input_fasta file to perform the translation')
@click.option('-t', '--translation_table', help='Translation table default value 1')
@click.option('-out', '--output', help='Output File')
@click.pass_context
def threeframe_translation(ctx, config_file, input_fasta, translation_table, output):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_fasta is None:
        print_help()

    pipeline_arguments = {}

    if translation_table is not None:
        pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table

    if output is not None:
        pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.three_frame_translation(input_fasta)
