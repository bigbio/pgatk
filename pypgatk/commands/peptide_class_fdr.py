import click
import logging
from pypgatk.commands.utils import print_help

import pkgutil

from pypgatk.proteomics.openms import OpenmsDataService
from pypgatk.toolbox.general import read_yaml_from_text, read_yaml_from_file, parse_peptide_classes, \
  parse_peptide_groups

log = logging.getLogger(__name__)

try:
  default_config_text = pkgutil.get_data(__name__, "../config/openms_analysis.yaml").decode()
except ValueError:
  try:
    default_config_text = pkgutil.get_data(__name__, "config/openms_analysis.yaml").decode()
  except ValueError:
    log.info("Configuration file not available !!!")


@click.command('peptide-class-fdr', short_help="Command to compute the Peptide class FDR")
@click.option('-c', '--config_file', help='Configuration to perform Peptide Class FDR')
@click.option('-in', '--input-file', help='input file with the peptides and proteins')
@click.option('-out', '--output-file', help='idxml from openms with filtered peptides and proteins')
@click.option("--file-type")
@click.option('--min-peptide-length', help='minimum peptide length')
@click.option('--psm-pep-fdr-cutoff', help="PSM peptide FDR cutoff or threshold")
@click.option('--psm-pep-class-fdr-cutoff', help="PSM class peptide FDR cutoff or threshold")
@click.option('--peptide-groups-prefix', help="Peptide class "
              "groups e.g. \"{non_canonical:[altorf,pseudo,ncRNA];mutations:[COSMIC,cbiomut];variants:[var_mut,var_rs]}\"")
@click.option('--peptide-classes-prefix', help='Peptides classes e.g. \"altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs\"')
@click.option("--disable-class-fdr", help="Disable Class-FDR, only compute Global FDR", is_flag = True)
@click.pass_context
def peptide_class_fdr(ctx, config_file, input_file, output_file, file_type, min_peptide_length, psm_pep_fdr_cutoff, psm_pep_class_fdr_cutoff,
                      peptide_groups_prefix, peptide_classes_prefix, disable_class_fdr):
  """
  The peptide_class_fdr allows to filter the peptide psm files (IdXML files) using two different FDR threshold types:
   - Global FDR
   - Global FDR + Peptide Class FDR
  The peptide classes can be defined in two ways as simple class:
   - "altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs"
  where each class represent only one kind of peptide source pseudo gene, ncRNA, etc. The second for of representing
  peptide classes is using groups of classes:
   - "{non_canonical:[altorf,pseudo,ncRNA];mutations:[COSMIC,cbiomut];variants:[var_mut,var_rs]}"
  in this case a class is a group of peptide sources for example: mutations with two difference sources as COSMIC
  and cbiomut (CBioportal mutation) .

  :param ctx:
  :param config_file: Configuration file
  :param input_file: Input idXML/Triqler containing peptide identifications
  :param output_file: Output idXML/Triqler containing peptide identifications after filtering
  :param min_peptide_length: Minimum peptide length
  :param psm_pep_fdr_cutoff:  Global FDR cutoff
  :param psm_pep_class_fdr_cutoff: Peptide class FDR cutoff
  :param peptide_groups_prefix: Peptide groups prefix for the Peptide classes FDR
  :param peptide_classes_prefix: Peptide classes
  :param file_type: File type to compute the FDR and class FDR.
  :param disable_class_fdr: Do not compute class FDR and not filtering the PSMs
  :return:
  """

  if config_file is None:
    config_data = read_yaml_from_text(default_config_text)
    msg = "The default configuration file is used: {}".format("openms_analysis.yaml")
    log.info(msg)
  else:
    config_data = read_yaml_from_file(config_file)

  if input_file is None or output_file is None:
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

  if disable_class_fdr is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_PEPTIDE_DISABLE_CLASS_FDR] = disable_class_fdr

  if file_type is not None:
    pipeline_arguments[OpenmsDataService.CONFIG_FILE_TYPE] = file_type

  openms_analyzer = OpenmsDataService(config_data, pipeline_arguments)
  openms_analyzer.filter_peptide_class_fdr(input_file, output_file)
