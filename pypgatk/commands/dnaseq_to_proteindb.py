import logging

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService
from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


@click.command("dnaseq-to-proteindb", short_help="Generate peptides based on DNA sequences")
@click.option('-c', '--config_file', help='Configuration to perform conversion between ENSEMBL Files')
@click.option('--input_fasta', help='Path to sequences fasta')
@click.option('--translation_table', type=int, help='Translation Table (default 1)')
@click.option('--num_orfs', type=int, help='Number of ORFs (default 3)')
@click.option('--num_orfs_complement', type=int,
              help='Number of ORFs from the reverse side (default 0)')
@click.option('--output_proteindb', help="Output file name, exits if already exists")
@click.option('-p', '--var_prefix', help="String to add before the variant peptides")
@click.option('--skip_including_all_cds',
              help="By default any transcript that has a defined CDS will be translated, this option disables this features instead it only depends on the biotypes",
              is_flag=True)
@click.option('--include_biotypes', help="Include Biotypes")
@click.option('--exclude_biotypes', help="Exclude Biotypes")
@click.option('--biotype_str', type=str,
              help='String used to identify gene/transcript biotype in the gtf file.')
@click.option('--transcript_description_sep', type=str,
              help='Separator used to separate features in the fasta headers, usually either (space, / or semicolon).')
@click.option('--expression_str', type=str,
              help='String to be used for extracting expression value (TPM, FPKM, etc).')
@click.option('--expression_thresh', type=float,
              help='Threshold used to filter transcripts based on their expression values')
@click.pass_context
def dnaseq_to_proteindb(ctx, config_file, input_fasta, translation_table, num_orfs, num_orfs_complement,
                        output_proteindb, var_prefix,
                        skip_including_all_cds, include_biotypes, exclude_biotypes, biotype_str,
                        transcript_description_sep, expression_str, expression_thresh):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_fasta is None:
        print_help()

    pipeline_arguments = {}

    if translation_table is not None:
        pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table

    if output_proteindb is not None:
        pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output_proteindb

    if var_prefix is not None:
        pipeline_arguments[EnsemblDataService.HEADER_VAR_PREFIX] = var_prefix

    if exclude_biotypes is not None:
        pipeline_arguments[EnsemblDataService.EXCLUDE_BIOTYPES] = exclude_biotypes

    if skip_including_all_cds is not None:
        pipeline_arguments[EnsemblDataService.SKIP_INCLUDING_ALL_CDS] = skip_including_all_cds

    if include_biotypes is not None:
        pipeline_arguments[EnsemblDataService.INCLUDE_BIOTYPES] = include_biotypes

    if biotype_str is not None:
        pipeline_arguments[EnsemblDataService.BIOTYPE_STR] = biotype_str

    if transcript_description_sep is not None:
        pipeline_arguments[EnsemblDataService.TRANSCRIPT_DESCRIPTION_SEP] = transcript_description_sep

    if num_orfs is not None:
        pipeline_arguments[EnsemblDataService.NUM_ORFS] = num_orfs

    if num_orfs_complement is not None:
        pipeline_arguments[EnsemblDataService.NUM_ORFS_COMPLEMENT] = num_orfs_complement

    if expression_str is not None:
        pipeline_arguments[EnsemblDataService.EXPRESSION_STR] = expression_str

    if expression_thresh is not None:
        pipeline_arguments[EnsemblDataService.EXPRESSION_THRESH] = expression_thresh

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.dnaseq_to_proteindb(input_fasta)
