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
        log.info("No configuration file available -- ")


@click.command('generate-deeplc', short_help="Generate input for deepLC tool from idXML,mzTab or consensusXML")
@click.option('-c', '--config_file', help='Configuration to perform deepLC configuration file')
@click.option('-in', '--input-file', help='input idxml/consensusxml/mztab containing the peptide identifications')
@click.option('-out', '--output-file', help='Output for DeepLC tool')
@click.option('-d', '--decoy_prefix', help='Set accession prefix for decoy proteins in output. Default=DECOY_',
              default='DECOY_')
@click.option('--peptide-classes-prefix',
              help='Peptides classes e.g. \"altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs\"')
@click.option('--novel-peptides', help='This parameter will allow to remove from the peptide list the novel peptides',
              is_flag=True)
@click.pass_context
def generate_deeplc(ctx, config_file, input_file, output_file, decoy_prefix: str, peptide_classes_prefix: str,
                    novel_peptides: bool):
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
    openms_analyzer._generate_deepLC_file(input_xml=input_file, output_deepLC=output_file, decoy_pattern=decoy_prefix,
                                          peptide_class_prefix=peptide_classes_prefix, novel_peptides=novel_peptides)
