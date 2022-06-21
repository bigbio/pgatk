import click
import logging
from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService

import pkgutil

from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file

log = logging.getLogger(__name__)
try:
    default_config_text = pkgutil.get_data(__name__, "../config/ensembl_config.yaml").decode()
except ValueError:
    try:
        default_config_text = pkgutil.get_data(__name__, "config/ensembl_config.yaml").decode()
    except ValueError:
        log.info("Configuration file not available !! ")


@click.command('threeframe-translation', short_help="Command to perform 3'frame translation")
@click.option('-c', '--config_file',
              help='Configuration to perform conversion between ENSEMBL Files')
@click.option('-in', '--input_fasta', help='input_fasta file to perform the translation')
@click.option('-t', '--translation_table', help='Translation table default value 1', default='1')
@click.option('-out', '--output', help='Output File', default="peptide-database.fa")
@click.pass_context
def threeframe_translation(ctx, config_file, input_fasta, translation_table, output):
    if config_file is None:
        config_data = read_yaml_from_text(default_config_text)
        msg = "The default configuration file is used: {}".format("ensembl_config.yaml")
        log.info(msg)
    else:
        config_data = read_yaml_from_file(config_file)

    if input_fasta is None:
        print_help()

    pipeline_arguments = {EnsemblDataService.TRANSLATION_TABLE: translation_table,
                          EnsemblDataService.PROTEIN_DB_OUTPUT: output}

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.three_frame_translation(input_fasta)
