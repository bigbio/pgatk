import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.commands.utils import print_help
from pypgatk.proteogenomics.mztab_class_fdr import MzTabClassFdr

log = logging.getLogger(__name__)

@click.command('mztab_class_fdr', short_help='Extract psms from mzTab for global-fdr and class-fdr filtering')
@click.option('-c', '--config_file', help='Configuration file for the fdr peptides pipeline')
@click.option('-i', '--input_mztab', help='The file name of the input mzTab')
@click.option('-o', '--outfile_name', help='The file name of the psm table filtered by global-fdr and class-fdr')
@click.option('-d', '--decoy_prefix', help='Default is "decoy"')
@click.option('-gf', '--global_fdr_cutoff', help='PSM peptide global-fdr cutoff or threshold. Default is 0.01')
@click.option('-cf', '--class_fdr_cutoff', help='PSM peptide class-fdr cutoff or threshold. Default is 0.01')
@click.option('-g', '--peptide_groups_prefix', help="Peptide class "
                                              "groups e.g. \"{non_canonical:[altorf,pseudo,ncRNA];mutations:[COSMIC,cbiomut];variants:[var_mut,var_rs]}\"")
@click.pass_context
def mztab_class_fdr(ctx, config_file, input_mztab, outfile_name, decoy_prefix, global_fdr_cutoff, class_fdr_cutoff, peptide_groups_prefix):
    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_mztab is None or outfile_name is None:
        print_help()
    pipeline_arguments = {}
    if decoy_prefix is not None:
        pipeline_arguments[MzTabClassFdr.CONFIG_DECOY_PREFIX] = decoy_prefix
    if global_fdr_cutoff is not None:
        pipeline_arguments[MzTabClassFdr.CONFIG_GLOBAL_FDR_CUTOFF] = global_fdr_cutoff
    if class_fdr_cutoff is not None:
        pipeline_arguments[MzTabClassFdr.CONFIG_CLASS_FDR_CUTOFF] = class_fdr_cutoff
    if peptide_groups_prefix is not None:
        pipeline_arguments[MzTabClassFdr.CONFIG_PEPTIDE_GROUPS_PREFIX] = peptide_groups_prefix

    mzTab_class_fdr = MzTabClassFdr(config_data, pipeline_arguments)
    mzTab_class_fdr.form_mztab_class_fdr(input_mztab, outfile_name)
