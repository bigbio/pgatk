import click
import logging

from pypgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pypgatk.commands.utils import print_help
import pkgutil

from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file

log = logging.getLogger(__name__)

try:
    default_config_text = pkgutil.get_data(__name__, "../config/cbioportal_config.yaml").decode()
except Exception:
    try:
        default_config_text = pkgutil.get_data(__name__, "config/cbioportal_config.yaml").decode()
    except ValueError:
        log.info("The config file is now available !! ")


@click.command('cbioportal-to-proteindb', short_help='Command to translate cbioportal mutation data into proteindb')
@click.option('-c', '--config_file', help='Configuration for cbioportal to proteindb tool')
@click.option('-in', '--input_mutation', help='Cbioportal mutation file')
@click.option('-fa', '--input_cds', help='CDS genes from ENSEMBL database')
@click.option('-out', '--output_db', help='Protein database including all the mutations')
@click.option('-f', '--filter_column', help='Column in the VCF file to be used for filtering or splitting mutations')
@click.option('-a', '--accepted_values',
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.option('-cl', '--clinical_sample_file',
              help='File to get tissue type from for the samples in input_mutation file')
@click.pass_context
def cbioportal_to_proteindb(ctx, config_file, input_mutation, input_cds, output_db,
                            clinical_sample_file, filter_column, accepted_values, split_by_filter_column):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file for cbioportal to proteindb is used: {}".format("cbioportal_config.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_mutation is None or input_cds is None or output_db is None:
    print_help()

  pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                        CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_cds,
                        CancerGenomesService.CONFIG_OUTPUT_FILE: output_db,
                        CancerGenomesService.CLINICAL_SAMPLE_FILE: clinical_sample_file,
                        CancerGenomesService.FILTER_COLUMN: filter_column,
                        CancerGenomesService.ACCEPTED_VALUES: accepted_values,
                        CancerGenomesService.SPLIT_BY_FILTER_COLUMN: split_by_filter_column}

  cosmic_to_proteindb_service = CancerGenomesService(config_data, pipeline_arguments)
  cosmic_to_proteindb_service.cbioportal_to_proteindb()
