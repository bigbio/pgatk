import logging
import sys

import click

from pgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pgatk.config.registry import load_config


@click.command('cbioportal-to-proteindb', short_help='Command to translate cbioportal mutation data into proteindb')
@click.option('-c', '--config_file', help='Configuration for cbioportal to proteindb tool')
@click.option('-in', '--input_mutation', help='Cbioportal mutation file', required=True)
@click.option('-fa', '--input_fasta',
              help='Transcript FASTA with CDS= and gene_biotype= headers (from ncbi-downloader --generate-transcripts)',
              required=True)
@click.option('-out', '--output_db', help='Protein database including all the mutations', required=True)
@click.option('-f', '--filter_column', help='Column in the VCF file to be used for filtering or splitting mutations')
@click.option('-a', '--accepted_values',
              help='Limit mutations to values (tissue type, sample name, etc) considered for generating proteinDBs, by default mutations from all records are considered (e.g. "")')
@click.option('-s', '--split_by_filter_column',
              help='Use this flag to generate a proteinDB per group as specified in the filter_column, default is False',
              is_flag=True)
@click.option('-cl', '--clinical_sample_file',
              help='File to get tissue type from for the samples in input_mutation file')
@click.option('--include_biotypes',
              help='Comma-separated biotypes to include (default: protein_coding). Use "all" to include all biotypes.')
@click.option('--exclude_biotypes',
              help='Comma-separated biotypes to exclude (default: none)')
@click.option('--skip_including_all_cds', is_flag=True, default=False,
              help='Apply biotype filter to CDS-defined transcripts as well (default: CDS transcripts are always included)')
@click.option('--include_variant_classifications',
              help='Comma-separated Variant_Classification values to include (default: all)')
@click.option('--exclude_variant_classifications',
              help='Comma-separated Variant_Classification values to exclude (default: Nonsense_Mutation)')
@click.option('--gff', 'gff_file',
              help='GFF3/GTF annotation file (from ncbi-downloader) for coordinate-based mutation application '
                   'when HGVSc is absent. Indexed as <gff>.db on first use.')
@click.option('-w', '--workers', type=int, default=1, show_default=True,
              help='Number of parallel worker processes for the translation phase (default: 1 = sequential).')
@click.pass_context
def cbioportal_to_proteindb(ctx, config_file, input_mutation, input_fasta, output_db,
                            clinical_sample_file, filter_column, accepted_values, split_by_filter_column,
                            include_biotypes, exclude_biotypes, skip_including_all_cds,
                            include_variant_classifications, exclude_variant_classifications, gff_file,
                            workers):

    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stderr,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    config_data = load_config("cbioportal", config_file)

    pipeline_arguments = {}

    if input_mutation is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE] = input_mutation

    if input_fasta is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_COMPLETE_GENES_FILE] = input_fasta

    if output_db is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_OUTPUT_FILE] = output_db

    if clinical_sample_file is not None:
        pipeline_arguments[CancerGenomesService.CLINICAL_SAMPLE_FILE] = clinical_sample_file

    if filter_column is not None:
        pipeline_arguments[CancerGenomesService.FILTER_COLUMN] = filter_column
    elif config_data is None:
        pipeline_arguments[CancerGenomesService.FILTER_COLUMN] = 'CANCER_TYPE'

    if accepted_values is not None:
        pipeline_arguments[CancerGenomesService.ACCEPTED_VALUES] = accepted_values

    if split_by_filter_column is not None:
        pipeline_arguments[CancerGenomesService.SPLIT_BY_FILTER_COLUMN] = split_by_filter_column

    if include_biotypes is not None:
        pipeline_arguments[CancerGenomesService.INCLUDE_BIOTYPES] = include_biotypes

    if exclude_biotypes is not None:
        pipeline_arguments[CancerGenomesService.EXCLUDE_BIOTYPES] = exclude_biotypes

    if skip_including_all_cds:
        pipeline_arguments[CancerGenomesService.SKIP_INCLUDING_ALL_CDS] = True

    if include_variant_classifications is not None:
        pipeline_arguments[CancerGenomesService.INCLUDE_VARIANT_CLASSIFICATIONS] = include_variant_classifications

    if exclude_variant_classifications is not None:
        pipeline_arguments[CancerGenomesService.EXCLUDE_VARIANT_CLASSIFICATIONS] = exclude_variant_classifications

    if gff_file is not None:
        pipeline_arguments[CancerGenomesService.CONFIG_GFF_FILE] = gff_file

    if workers is not None and workers > 1:
        pipeline_arguments[CancerGenomesService.WORKERS] = workers

    cosmic_to_proteindb_service = CancerGenomesService(config_data, pipeline_arguments)
    cosmic_to_proteindb_service.cbioportal_to_proteindb()
