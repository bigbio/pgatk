import logging
import click

from pypgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pypgatk.commands.utils import print_help
import pkgutil

from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file

log = logging.getLogger(__name__)
try:
    default_config_text = pkgutil.get_data(__name__, "../config/cosmic_config.yaml").decode()
except ValueError:
    try:
        default_config_text = pkgutil.get_data(__name__, "config/cosmic_config.yaml").decode()
    except ValueError:
        log.info("Configuration file not available !!!")

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
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.pass_context
def cosmic_to_proteindb(ctx, config_file, input_mutation, input_genes, output_db,
                        filter_column, accepted_values, split_by_filter_column):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("cosmic_config.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_mutation is None or input_genes is None or output_db is None:
    print_help()

  pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                        CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_genes,
                        CancerGenomesService.CONFIG_OUTPUT_FILE: output_db,
                        CancerGenomesService.FILTER_COLUMN: filter_column,
                        CancerGenomesService.ACCEPTED_VALUES: accepted_values,
                        CancerGenomesService.SPLIT_BY_FILTER_COLUMN: split_by_filter_column}

  cosmic_to_proteindb_service = CancerGenomesService(config_data, pipeline_arguments)
  cosmic_to_proteindb_service.cosmic_to_proteindb()
