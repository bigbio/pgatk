import os

import click

from pgatk.clinvar.data_downloader import NcbiDataDownloader
from pgatk.config.registry import load_config


@click.command("ncbi-downloader", short_help="Download NCBI RefSeq and ClinVar reference files")
@click.option("-c", "--config_file", help="Configuration YAML file (optional)")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.option(
    "--generate-transcripts",
    "generate_transcripts",
    is_flag=True,
    default=False,
    help=(
        "Run gffread after download to extract transcript sequences with CDS= "
        "coordinate headers (required for clinvar-to-proteindb). Requires gffread "
        "in PATH (conda install -c bioconda gffread)."
    ),
)
@click.pass_context
def ncbi_downloader(ctx, config_file, output_dir, force, generate_transcripts):
    """Download NCBI RefSeq GFF3 annotation, genomic FASTA, assembly report, and ClinVar VCF.

    Files are downloaded to the specified output directory. Existing files
    are skipped unless --force is used.

    After downloading, use --generate-transcripts to run gffread and produce
    transcripts.fa with CDS= headers needed by clinvar-to-proteindb.
    """
    config = load_config("clinvar", config_file)
    clinvar_cfg = config.get("clinvar_translation", {})

    refseq_base = clinvar_cfg.get(
        "ncbi_refseq_ftp",
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/",
    )
    clinvar_base = clinvar_cfg.get(
        "clinvar_ftp",
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/",
    )

    downloader = NcbiDataDownloader(
        output_dir=output_dir,
        refseq_base_url=refseq_base,
        clinvar_base_url=clinvar_base,
    )
    downloaded = downloader.download_all(force=force)
    click.echo(f"Downloaded {len(downloaded)} files to {output_dir}")

    if generate_transcripts:
        genome_fna = os.path.join(output_dir, "GRCh38_latest_genomic.fna")
        gff_file = os.path.join(output_dir, "GRCh38_latest_genomic.gff")
        output_transcripts = os.path.join(output_dir, "transcripts.fa")
        try:
            click.echo("Running gffread to generate transcripts.fa with CDS= headers...")
            NcbiDataDownloader.generate_transcripts(genome_fna, gff_file, output_transcripts)
            click.echo(f"Transcripts written to {output_transcripts}")
        except FileNotFoundError as exc:
            click.echo(f"Error: {exc}", err=True)
            raise SystemExit(1)
