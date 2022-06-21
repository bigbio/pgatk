import click
import logging
from pypgatk.commands.utils import print_help

import pkgutil

from pypgatk.proteomics.openms import OpenmsDataService
from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file

log = logging.getLogger(__name__)
try:
  default_config_text = pkgutil.get_data(__name__, "../config/openms_analysis.yaml").decode()
except ValueError:
  try:
    default_config_text = pkgutil.get_data(__name__, "config/openms_analysis.yaml").decode()
  except ValueError:
    log.info("Configuration file not available !!")

@click.command('msrescore-configuration', short_help="Command to generate the msrescore configuration file from idXML")
@click.option('-c', '--config_file', help='Configuration to perform msrescore configuration file')
@click.option('-in', '--input-file', help='input idxml containing the peptide identifications')
@click.option('-out', '--output-file', help='Output json configuration file for msrescore')
@click.option('--quant_method', help = 'Quantification method TMT or LFQ', type=click.Choice(['LFQ', 'TMT']))
@click.option('-d', '--decoy_prefix', help='Set accession prefix for decoy proteins in output. Default=DECOY_', default = 'DECOY_')
@click.pass_context
def msrescore_configuration(ctx, config_file, input_file, output_file, quant_method, decoy_prefix: str):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("openms_analysis.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_file is None or output_file is None:
    print_help()

  pipeline_arguments = {}

  openms_analyzer = OpenmsDataService(config_data, pipeline_arguments)
  openms_analyzer._generate_msrescore_file(input_xml = input_file, quant_method= quant_method, output_json = output_file, decoy_pattern = decoy_prefix)
