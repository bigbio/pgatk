import os

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService

this_dir, this_filename = os.path.split(__file__)


@click.command("dnaseq-to-proteindb", short_help="Generate peptides based on DNA sequences")
@click.option('-c', '--config_file', help='Configuration to perform conversion between ENSEMBL Files',
              default= this_dir + '/../config/ensembl_config.yaml')
@click.option('--input_fasta', help='Path to sequences fasta')
@click.option('--translation_table', default=1, type=int, help='Translation Table (default 1)')
@click.option('--num_orfs', default=3, type=int, help='Number of ORFs (default 0)')
@click.option('--num_orfs_complement', default=0, type=int,
              help='Number of ORFs from the reverse side (default 0)')
@click.option('--output_proteindb', default="peptide-database.fa", help="Output file name, exits if already exists")
@click.option('-p', '--var_prefix', default="", help="String to add before the variant peptides")

@click.option('--skip_including_all_cds',
              help="By default any transcript that has a defined CDS will be translated, this option disables this features instead it only depends on the biotypes",
              is_flag=True)
@click.option('--include_biotypes', default='', help="Include Biotypes")
@click.option('--exclude_biotypes', default='', help="Exclude Biotypes")
@click.option('--biotype_str', default='transcript_biotype', type=str,
              help='String used to identify gene/transcript biotype in the gtf file.')
@click.option('--expression_str', default="", type=str,
              help='String to be used for extracting expression value (TPM, FPKM, etc).')
@click.option('--expression_thresh', default=5.0, type=float,
              help='Threshold used to filter transcripts based on their expression values')
@click.pass_context
def dnaseq_to_proteindb(ctx, config_file, input_fasta, translation_table, num_orfs, num_orfs_complement,
                        output_proteindb, var_prefix,
                        skip_including_all_cds, include_biotypes, exclude_biotypes, biotype_str, expression_str,
                        expression_thresh):
    if input_fasta is None:
        print_help()

    pipeline_arguments = {EnsemblDataService.TRANSLATION_TABLE: translation_table,
                          EnsemblDataService.PROTEIN_DB_OUTPUT: output_proteindb,
                          EnsemblDataService.HEADER_VAR_PREFIX: var_prefix,
                          EnsemblDataService.EXCLUDE_BIOTYPES: exclude_biotypes,
                          EnsemblDataService.SKIP_INCLUDING_ALL_CDS: skip_including_all_cds,
                          EnsemblDataService.INCLUDE_BIOTYPES: include_biotypes,
                          EnsemblDataService.BIOTYPE_STR: biotype_str, EnsemblDataService.NUM_ORFS: num_orfs,
                          EnsemblDataService.NUM_ORFS_COMPLEMENT: num_orfs_complement,
                          EnsemblDataService.EXPRESSION_STR: expression_str,
                          EnsemblDataService.EXPRESSION_THRESH: expression_thresh}

    ensembl_data_service = EnsemblDataService(config_file, pipeline_arguments)
    ensembl_data_service.dnaseq_to_proteindb(input_fasta)

