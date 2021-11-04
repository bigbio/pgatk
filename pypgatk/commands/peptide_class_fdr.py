import click
import logging
from pypgatk.commands.utils import print_help

import pkgutil

from pypgatk.proteomics.openms import OpenmsDataService
from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file
default_config_text = pkgutil.get_data(__name__, "../config/openms_analysis.yaml").decode()


@click.command('peptide-class-fdr', short_help="Command to compute the Peptide class FDR")
@click.option('-c', '--config_file', help='Configuration to perform Peptide Class FDR')
@click.option('-in', '--input_idxml', help='idxml from openms with peptides and proteins')
@click.option('-out', '--output_idxml', help='idxml from openms with filtered peptides and proteins')
@click.pass_context
def peptide_class_fdr(ctx, config_file, input_idxml, output_idxml):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("ensembl_config.yaml")
    logging.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_idxml is None or output_idxml is None:
    print_help()

  pipeline_arguments = {}


  openms_analyzer = OpenmsDataService(config_data, pipeline_arguments)
  openms_analyzer.filter_peptide_class_fdr(input_idxml, output_idxml)
