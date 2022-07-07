import logging

import click

from pypgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pypgatk.commands.utils import print_help
from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


@click.command('cosmic-to-proteindb', short_help='Command to translate Cosmic mutation data into proteindb')
@click.option('-c', '--config_file',
              help='Configuration file for the cosmic data pipelines')
@click.option('-in', '--input_mutation', help='Cosmic Mutation data file')
@click.option('-fa', '--input_genes', help='All Cosmic genes')
@click.option('-out', '--output_db',
              help='Protein database including all the mutations')
@click.option('-f', '--filter_column',
              help='Column in the VCF file to be used for filtering or splitting mutations')
@click.option('-a', '--accepted_values',
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered (e.g. "")')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.pass_context
def cosmic_to_proteindb(ctx, config_file, input_mutation, input_genes, output_db,
                        filter_column, accepted_values, split_by_filter_column):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_mutation is None or input_genes is None or output_db is None:
        print_help()

    pipeline_arguments = {}

    if input_mutation is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE] = input_mutation

    if input_genes is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_COMPLETE_GENES_FILE] = input_genes

    if output_db is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_OUTPUT_FILE] = output_db

    if filter_column is not None:
        pipeline_arguments[CancerGenomesService.FILTER_COLUMN] = filter_column
    elif config_data is None:
        pipeline_arguments[CancerGenomesService.FILTER_COLUMN] = 'Primary site'

    if accepted_values is not None:
        pipeline_arguments[CancerGenomesService.ACCEPTED_VALUES] = accepted_values

    if split_by_filter_column is not None:
        pipeline_arguments[CancerGenomesService.SPLIT_BY_FILTER_COLUMN] = split_by_filter_column

    cosmic_to_proteindb_service = CancerGenomesService(config_data, pipeline_arguments)
    cosmic_to_proteindb_service.cosmic_to_proteindb()
