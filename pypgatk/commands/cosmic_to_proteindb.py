import os

import click

from pypgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pypgatk.commands.utils import print_help

this_dir, this_filename = os.path.split(__file__)


@click.command('cosmic-to-proteindb', short_help='Command to translate Cosmic mutation data into proteindb')
@click.option('-c', '--config_file',
              help='Configuration file for the cosmic data pipelines',
              default=this_dir + '/../config/cosmic_config.yaml')
@click.option('-in', '--input_mutation', help='Cosmic Mutation data file')
@click.option('-fa', '--input_genes', help='All Cosmic genes')
@click.option('-out', '--output_db',
              help='Protein database including all the mutations')
@click.option('-f', '--filter_column', default='Primary site',
              help='Column in the VCF file to be used for filtering or splitting mutations')
@click.option('-a', '--accepted_values', default='all',
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.pass_context
def cosmic_to_proteindb(ctx, config_file, input_mutation, input_genes, output_db,
                        filter_column, accepted_values, split_by_filter_column):
  if input_mutation is None or input_genes is None or output_db is None:
    print_help()

  pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                        CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_genes,
                        CancerGenomesService.CONFIG_OUTPUT_FILE: output_db,
                        CancerGenomesService.FILTER_COLUMN: filter_column,
                        CancerGenomesService.ACCEPTED_VALUES: accepted_values,
                        CancerGenomesService.SPLIT_BY_FILTER_COLUMN: split_by_filter_column}

  cosmic_to_proteindb_service = CancerGenomesService(config_file, pipeline_arguments)
  cosmic_to_proteindb_service.cosmic_to_proteindb()
