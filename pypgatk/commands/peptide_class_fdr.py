import click
import logging
from pypgatk.commands.utils import print_help

import pkgutil

from pypgatk.proteomics.openms import OpenmsDataService
from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file, parse_peptide_classes, \
  parse_peptide_groups

default_config_text = pkgutil.get_data(__name__, "../config/openms_analysis.yaml").decode()

log = logging.getLogger(__name__)

@click.command('peptide-class-fdr', short_help="Command to compute the Peptide class FDR")
@click.option('-c', '--config_file', help='Configuration to perform Peptide Class FDR')
@click.option('-in', '--input_idxml', help='idxml from openms with peptides and proteins')
@click.option('-out', '--output_idxml', help='idxml from openms with filtered peptides and proteins')
@click.option('--min-peptide-length', help='minimum peptide length')
@click.option('--psm-pep-fdr-cutoff', help="PSM peptide FDR cutoff or threshold", default=0.01)
@click.option('--psm-pep-class-fdr-cutoff', help="PSM class peptide FDR cutoff or threshold", default=0.01)
@click.option('--peptide_groups_prefix', help="Peptide class "
              "groups e.g. \"{non_canonical:[altorf,pseudo,ncRNA];mutations:[COSMIC,cbiomut];variants:[var_mut,var_rs]}\"")
@click.option('--peptide_classes_prefix', help='Peptides classes e.g. \"altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs\"')
@click.option("--enable_class_fdr", help="Enable Class-FDR over Global PSM FDR (default true)", default = True)
@click.pass_context
def peptide_class_fdr(ctx, config_file, input_idxml, output_idxml, min_peptide_length, psm_pep_fdr_cutoff, psm_pep_class_fdr_cutoff,
                      peptide_groups_prefix, peptide_classes_prefix, enable_class_fdr):

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("openms_analysis.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_idxml is None or output_idxml is None:
    print_help()

  pipeline_arguments = {}

  if min_peptide_length is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_MIN_PEPTIDE_LENGTH] = min_peptide_length

  if psm_pep_fdr_cutoff is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_FDR_CUTOFF] = psm_pep_fdr_cutoff

  if peptide_classes_prefix is not None and peptide_groups_prefix is not None:
    raise ValueError("The tool can't be use to compute class groups and classes FDR at the same time")

  if peptide_classes_prefix is not None:
    data = parse_peptide_classes(peptide_classes_prefix)
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_GROUP_PREFIX] = data

  if psm_pep_class_fdr_cutoff is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_CLASS_FDR_CUTOFF] = psm_pep_class_fdr_cutoff

  if peptide_groups_prefix is not None:
    data = parse_peptide_groups(peptide_groups_prefix)
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_GROUP_PREFIX] = data

  if enable_class_fdr is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_APPLY_CLASS_FDR] = enable_class_fdr

  openms_analyzer = OpenmsDataService(config_data, pipeline_arguments)
  openms_analyzer.filter_peptide_class_fdr(input_idxml, output_idxml)
