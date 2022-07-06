import logging

import click

from pypgatk.commands.utils import print_help
from pypgatk.ensembl.ensembl import EnsemblDataService
from pypgatk.toolbox.general import read_yaml_from_file

log = logging.getLogger(__name__)


@click.command('vcf-to-proteindb', short_help="Generate peptides based on DNA variants VCF files")
@click.option('-c', '--config_file', help='Configuration to perform conversion between ENSEMBL Files')
@click.option('-f', '--input_fasta', help='Path to the transcript sequence')
@click.option('-v', '--vcf', help='Path to the VCF file')
@click.option('-g', '--gene_annotations_gtf', help='Path to the gene annotations file')
@click.option('-t', '--translation_table', type=int, help="Translation table (Default 1) ")
@click.option('-m', '--mito_translation_table', type=int, help='Mito_trans_table (default 2)')
@click.option('-p', '--var_prefix', default="var", help="String to add before the variant peptides")
@click.option('--report_ref_seq', help='In addition to var peps, also report all ref peps', is_flag=True)
@click.option('-o', '--output_proteindb', help="Output file name, exits if already exists")
@click.option('--annotation_field_name', help='''Annotation field name found in the INFO column,
              e.g CSQ or vep; if empty it will identify overlapping transcripts
              from the given GTF file and no aa consequence will be considered''')
@click.option('--af_field', help="field name in the VCF INFO column to use for filtering on AF, (Default None)")
@click.option('--af_threshold', help='Minium AF threshold for considering common variants')
@click.option('--transcript_str', type=str, help='String that is used for transcript ID in the VCF header INFO field')
@click.option('--biotype_str', type=str,
              help='String that is used for biotype in the VCF header INFO field')
@click.option('--exclude_biotypes', help="Excluded Biotypes")
@click.option('--include_biotypes', help="included_biotypes, default all")
@click.option('--consequence_str', type=str, help='String that is used for consequence in the VCF header INFO field')
@click.option('--exclude_consequences', help="Excluded Consequences")
@click.option('-s', '--skip_including_all_cds',
              help="by default any transcript that has a defined CDS will be used, this option disables this features instead",
              is_flag=True)
@click.option('--include_consequences', help="included_consequences, default all")
@click.option('--ignore_filters',
              help="enabling this option causes or variants to be parsed. By default only variants that have not failed any filters will be processed (FILTER column is PASS, None, .) or if the filters are subset of the accepted filters. (default is False)",
              is_flag=True)
@click.option('--accepted_filters', help="Accepted filters for variant parsing")
@click.pass_context
def vcf_to_proteindb(ctx, config_file, input_fasta, vcf, gene_annotations_gtf, translation_table,
                     mito_translation_table,
                     var_prefix, report_ref_seq, output_proteindb, annotation_field_name,
                     af_field, af_threshold, transcript_str, biotype_str, exclude_biotypes,
                     include_biotypes, consequence_str, exclude_consequences,
                     skip_including_all_cds, include_consequences,
                     ignore_filters, accepted_filters):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_fasta is None or vcf is None or gene_annotations_gtf is None:
        print_help()

    pipeline_arguments = {}

    if mito_translation_table is not None:
        pipeline_arguments[EnsemblDataService.MITO_TRANSLATION_TABLE] = mito_translation_table

    if translation_table is not None:
        pipeline_arguments[EnsemblDataService.TRANSLATION_TABLE] = translation_table

    if var_prefix is not None:
        pipeline_arguments[EnsemblDataService.HEADER_VAR_PREFIX] = var_prefix

    if report_ref_seq is not None:
        pipeline_arguments[EnsemblDataService.REPORT_REFERENCE_SEQ] = report_ref_seq

    if output_proteindb is not None:
        pipeline_arguments[EnsemblDataService.PROTEIN_DB_OUTPUT] = output_proteindb

    if annotation_field_name is not None:
        pipeline_arguments[EnsemblDataService.ANNOTATION_FIELD_NAME] = annotation_field_name

    if af_field is not None:
        pipeline_arguments[EnsemblDataService.AF_FIELD] = af_field

    if af_threshold is not None:
        pipeline_arguments[EnsemblDataService.AF_THRESHOLD] = af_threshold

    if transcript_str is not None:
        pipeline_arguments[EnsemblDataService.TRANSCRIPT_STR] = transcript_str

    if biotype_str is not None:
        pipeline_arguments[EnsemblDataService.BIOTYPE_STR] = biotype_str

    if exclude_biotypes is not None:
        pipeline_arguments[EnsemblDataService.EXCLUDE_BIOTYPES] = exclude_biotypes

    if include_biotypes is not None:
        pipeline_arguments[EnsemblDataService.INCLUDE_BIOTYPES] = include_biotypes

    if consequence_str is not None:
        pipeline_arguments[EnsemblDataService.CONSEQUENCE_STR] = consequence_str

    if exclude_consequences is not None:
        pipeline_arguments[EnsemblDataService.EXCLUDE_CONSEQUENCES] = exclude_consequences

    if skip_including_all_cds is not None:
        pipeline_arguments[EnsemblDataService.SKIP_INCLUDING_ALL_CDS] = skip_including_all_cds

    if include_consequences is not None:
        pipeline_arguments[EnsemblDataService.INCLUDE_CONSEQUENCES] = include_consequences

    if ignore_filters is not None:
        pipeline_arguments[EnsemblDataService.IGNORE_FILTERS] = ignore_filters

    if accepted_filters is not None:
        pipeline_arguments[EnsemblDataService.ACCEPTED_FILTERS] = accepted_filters

    ensembl_data_service = EnsemblDataService(config_data, pipeline_arguments)
    ensembl_data_service.vcf_to_proteindb(vcf, input_fasta, gene_annotations_gtf)
