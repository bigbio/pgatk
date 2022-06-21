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
        log.info("Configuration file not available !!!")


@click.command('vcf-to-proteindb', short_help="Generate peptides based on DNA variants VCF files")
@click.option('-c', '--config_file', help='Configuration to perform conversion between ENSEMBL Files')
@click.option('-f', '--input_fasta', help='Path to the transcript sequence')
@click.option('-v', '--vcf', help='Path to the VCF file')
@click.option('-g', '--gene_annotations_gtf', help='Path to the gene annotations file')
@click.option('-t', '--translation_table', default=1, type=int, help="Translation table (Default 1) ")
@click.option('-m', '--mito_translation_table', default=2, type=int, help='Mito_trans_table (default 2)')
@click.option('-p', '--var_prefix', default="var", help="String to add before the variant peptides")
@click.option('--report_ref_seq', help='In addition to var peps, also report all ref peps', is_flag=True)
@click.option('-o', '--output_proteindb', default="peptide-database.fa",
              help="Output file name, exits if already exists")
@click.option('--annotation_field_name', default="CSQ",
              help='''Annotation field name found in the INFO column,
              e.g CSQ or vep; if empty it will identify overlapping transcripts
              from the given GTF file and no aa consequence will be considered''')
@click.option('--af_field', default="",
              help="field name in the VCF INFO column to use for filtering on AF, (Default None)")
@click.option('--af_threshold', default=0.01, help='Minium AF threshold for considering common variants')
@click.option('--transcript_str', default='FEATURE', type=str,
              help='String that is used for transcript ID in the VCF header INFO field')
@click.option('--biotype_str', default='BIOTYPE', type=str,
              help='String that is used for biotype in the VCF header INFO field')
@click.option('--exclude_biotypes',
              default='',
              help="Excluded Biotypes", show_default=True)
@click.option('--include_biotypes',
              default='protein_coding,polymorphic_pseudogene,non_stop_decay,nonsense_mediated_decay,IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,TEC,mRNA',
              help="included_biotypes, default all")
@click.option('--consequence_str', default='CONSEQUENCE', type=str,
              help='String that is used for consequence in the VCF header INFO field')
@click.option('--exclude_consequences',
              default='downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant',
              help="Excluded Consequences", show_default=True)
@click.option('-s', '--skip_including_all_cds',
              help="by default any transcript that has a defined CDS will be used, this option disables this features instead",
              is_flag=True)
@click.option('--include_consequences', default='all', help="included_consequences, default all")
@click.option('--ignore_filters',
              help="enabling this option causes or variants to be parsed. By default only variants that have not failed any filters will be processed (FILTER column is PASS, None, .) or if the filters are subset of the accepted filters. (default is False)",
              is_flag=True)
@click.option('--accepted_filters', default='PASS', help="Accepted filters for variant parsing")
@click.pass_context
def vcf_to_proteindb(ctx, config_file, input_fasta, vcf, gene_annotations_gtf, translation_table,
                     mito_translation_table,
                     var_prefix, report_ref_seq, output_proteindb, annotation_field_name,
                     af_field, af_threshold, transcript_str, biotype_str, exclude_biotypes,
                     include_biotypes, consequence_str, exclude_consequences,
                     skip_including_all_cds, include_consequences,
                     ignore_filters, accepted_filters):
    if config_file is None:
        config_data = read_yaml_from_text(default_config_text)
        msg = "The default configuration file is used: {}".format("ensembl_config.yaml")
        log.info(msg)
    else:
        config_data = read_yaml_from_file(config_file)

    if input_fasta is None or vcf is None or gene_annotations_gtf is None:
        print_help()

    pipeline_arguments = {EnsemblDataService.MITO_TRANSLATION_TABLE: mito_translation_table,
                          EnsemblDataService.TRANSLATION_TABLE: translation_table,
                          EnsemblDataService.HEADER_VAR_PREFIX: var_prefix,
                          EnsemblDataService.REPORT_REFERENCE_SEQ: report_ref_seq,
                          EnsemblDataService.PROTEIN_DB_OUTPUT: output_proteindb,
                          EnsemblDataService.ANNOTATION_FIELD_NAME: annotation_field_name,
                          EnsemblDataService.AF_FIELD: af_field, EnsemblDataService.AF_THRESHOLD: af_threshold,
                          EnsemblDataService.TRANSCRIPT_STR: transcript_str,
                          EnsemblDataService.BIOTYPE_STR: biotype_str,
                          EnsemblDataService.EXCLUDE_BIOTYPES: exclude_biotypes,
                          EnsemblDataService.INCLUDE_BIOTYPES: include_biotypes,
                          EnsemblDataService.CONSEQUENCE_STR: consequence_str,
                          EnsemblDataService.EXCLUDE_CONSEQUENCES: exclude_consequences,
                          EnsemblDataService.SKIP_INCLUDING_ALL_CDS: skip_including_all_cds,
                          EnsemblDataService.INCLUDE_CONSEQUENCES: include_consequences,
                          EnsemblDataService.IGNORE_FILTERS: ignore_filters,
                          EnsemblDataService.ACCEPTED_FILTERS: accepted_filters}

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.vcf_to_proteindb(vcf, input_fasta, gene_annotations_gtf)
