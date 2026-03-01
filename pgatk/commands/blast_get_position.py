import logging

import click

from pgatk.toolbox.general import read_yaml_from_file
from pgatk.proteogenomics.blast_get_position import BlastGetPositionService

log = logging.getLogger(__name__)

@click.command('blast_get_position', short_help='Blast peptide and refence protein database to find variation sites.')
@click.option('-c', '--config_file', help='Configuration file for the fdr peptides pipeline.')
@click.option('-i', '--input_psm_to_blast', help='The file name of the input PSM table to blast.', required=True)
@click.option('-o', '--output_psm', help='The file name of the output PSM table.', required=True)
@click.option('-r', '--input_reference_database', help='The file name of the refence protein database to blast. The reference database includes Uniprot Proteomes with isoforms, ENSEMBL, RefSeq, etc.', required=True)
@click.option('-n', '--number_of_processes', help='Used to specify the number of processes. Default is 40.')

@click.pass_context
def blast_get_position(ctx, config_file, input_psm_to_blast, output_psm, input_reference_database, number_of_processes):
    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    pipeline_arguments = {}
    if input_reference_database is not None:
        pipeline_arguments[BlastGetPositionService.CONFIG_INPUT_REFERENCE_DATABASE] = input_reference_database
    if number_of_processes is not None:
        pipeline_arguments[BlastGetPositionService.CONFIG_NUMBER_OF_PROCESSES] = number_of_processes

    blast_get_position_service = BlastGetPositionService(config_data, pipeline_arguments)
    blast_get_position_service.blast(input_psm_to_blast, output_psm)