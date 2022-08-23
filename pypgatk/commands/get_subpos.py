import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.subpos.get_subpos import GetSubPos
from pypgatk.commands.utils import print_help


log = logging.getLogger(__name__)


@click.command('get_subpos', short_help='Get position in peptide indicate which amino acid is substituted')
@click.option('-c', '--config_file', help='Configuration file for the get subpos')
@click.option('-in', '--input_psm_table', help='Input variant peptide PSMs table')
@click.option('-i', '--input_fasta', help='Protein sequence used')
@click.option('-o', '--output_psm_table', help='Output variant peptide PSMs table')
@click.pass_context
def get_subpos(ctx, config_file, input_psm_table, input_fasta, output_psm_table):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_psm_table is None or input_fasta is None or output_psm_table is None:
        print_help()

    pipeline_arguments = {}  

    if input_psm_table is not None:
        pipeline_arguments[GetSubPos.CONFIG_INPUT_PSM_TABLE] = input_psm_table
    if input_fasta is not None:
        pipeline_arguments[GetSubPos.CONFIG_INPUT_FASTA] = input_fasta
    if output_psm_table is not None:
        pipeline_arguments[GetSubPos.CONFIG_OUTPUT_PSM_TABLE] = output_psm_table

    get_subpos_service = GetSubPos(config_data, pipeline_arguments)
    get_subpos_service.get_subpos(input_psm_table, input_fasta, output_psm_table)


