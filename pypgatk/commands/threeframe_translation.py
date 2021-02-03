import os

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService

this_dir, this_filename = os.path.split(__file__)


@click.command('threeframe-translation', short_help="Command to perform 3'frame translation")
@click.option('-c', '--config_file',
              help='Configuration to perform conversion between ENSEMBL Files',
              default= this_dir + '/../config/ensembl_config.yaml')
@click.option('-in', '--input_fasta', help='input_fasta file to perform the translation')
@click.option('-t', '--translation_table', help='Translation table default value 1', default='1')
@click.option('-out', '--output', help='Output File', default="peptide-database.fa")
@click.pass_context
def threeframe_translation(ctx, config_file, input_fasta, translation_table, output):
    if input_fasta is None:
        print_help()
    pipeline_arguments = {EnsemblDataService.TRANSLATION_TABLE: translation_table,
                          EnsemblDataService.PROTEIN_DB_OUTPUT: output}

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.three_frame_translation(input_fasta)
