import os

import click

from pypgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pypgatk.commands.utils import print_help

this_dir, this_filename = os.path.split(__file__)


@click.command('cbioportal-to-proteindb', short_help='Command to translate cbioportal mutation data into proteindb')
@click.option('-c', '--config_file',
              help='Configuration for cbioportal to proteindb tool',
              default=this_dir + '/../config/cbioportal_config.yaml')
@click.option('-in', '--input_mutation', help='Cbioportal mutation file')
@click.option('-fa', '--input_cds', help='CDS genes from ENSEMBL database')
@click.option('-out', '--output_db', help='Protein database including all the mutations')
@click.option('-f', '--filter_column', default='Tumor_Sample_Barcode',
              help='Column in the VCF file to be used for filtering or splitting mutations')
@click.option('-a', '--accepted_values', default='all',
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.option('-cl', '--clinical_sample_file',
              help='File to get tissue type from for the samples in input_mutation file')
@click.pass_context
def cbioportal_to_proteindb(ctx, config_file, input_mutation, input_cds, output_db,
                            clinical_sample_file, filter_column, accepted_values, split_by_filter_column):
  if input_mutation is None or input_cds is None or output_db is None:
    print_help()

  pipeline_arguments = {CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: input_mutation,
                        CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: input_cds,
                        CancerGenomesService.CONFIG_OUTPUT_FILE: output_db,
                        CancerGenomesService.CLINICAL_SAMPLE_FILE: clinical_sample_file,
                        CancerGenomesService.FILTER_COLUMN: filter_column,
                        CancerGenomesService.ACCEPTED_VALUES: accepted_values,
                        CancerGenomesService.SPLIT_BY_FILTER_COLUMN: split_by_filter_column}

  cosmic_to_proteindb_service = CancerGenomesService(config_file, pipeline_arguments)
  cosmic_to_proteindb_service.cbioportal_to_proteindb()
