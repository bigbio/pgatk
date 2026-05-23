import logging
import os
import sys

import click

from pgatk.clinvar.data_downloader import NcbiDataDownloader
from pgatk.config.registry import load_config


@click.command("ncbi-downloader", short_help="Download NCBI RefSeq and ClinVar reference files")
@click.option("-c", "--config_file", help="Configuration YAML file (optional)")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.option(
    "--grch37",
    is_flag=True,
    default=False,
    help="Download GRCh37 files instead of the default GRCh38.",
)
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
@click.option(
    "--generate-cds",
    "generate_cds",
    is_flag=True,
    default=False,
    help=(
        "Run gffread after download to extract CDS-only sequences into refseq_cds.fa "
        "(required for cbioportal-to-proteindb with RefSeq transcript IDs). "
        "Requires gffread in PATH (conda install -c bioconda gffread)."
    ),
)
@click.pass_context
def ncbi_downloader(ctx, config_file, output_dir, force, grch37, generate_transcripts, generate_cds):
    """Download NCBI RefSeq GFF3 annotation, genomic FASTA, assembly report, and ClinVar VCF.

    Supports both GRCh38 (default) and GRCh37 (--grch37). Existing files are
    skipped unless --force is used.

    \b
    For clinvar-to-proteindb:
      --generate-transcripts   produces transcripts.fa with CDS= headers

    \b
    For cbioportal-to-proteindb (RefSeq transcript IDs from the cBioPortal API):
      --generate-cds           produces refseq_cds.fa with spliced CDS sequences
    """
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stderr,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    config = load_config("clinvar", config_file)
    clinvar_cfg = config.get("clinvar_translation", {})

    refseq_base = clinvar_cfg.get("ncbi_refseq_ftp", None)
    clinvar_base = clinvar_cfg.get("clinvar_ftp", None)

    downloader = NcbiDataDownloader(
        output_dir=output_dir,
        **({"refseq_base_url": refseq_base} if refseq_base and not grch37 else {}),
        **({"clinvar_base_url": clinvar_base} if clinvar_base and not grch37 else {}),
        grch37=grch37,
    )
    downloaded = downloader.download_all(force=force)
    click.echo(f"Downloaded {len(downloaded)} files to {output_dir}")

    build = "GRCh37" if grch37 else "GRCh38"
    prefix = f"{build}_latest_genomic"
    genome_fna = os.path.join(output_dir, f"{prefix}.fna.gz")
    gff_file = os.path.join(output_dir, f"{prefix}.gff.gz")

    if generate_transcripts:
        output_transcripts = os.path.join(output_dir, "transcripts.fa")
        try:
            click.echo(f"Running gffread to generate transcripts.fa ({build}) ...")
            NcbiDataDownloader.generate_transcripts(genome_fna, gff_file, output_transcripts)
            click.echo(f"Transcripts written to {output_transcripts}")
        except (FileNotFoundError, Exception) as exc:
            click.echo(f"Error: {exc}", err=True)
            raise SystemExit(1)

    if generate_cds:
        output_cds = os.path.join(output_dir, "refseq_cds.fa")
        try:
            click.echo(f"Running gffread to generate refseq_cds.fa ({build}) ...")
            NcbiDataDownloader.generate_cds_fasta(genome_fna, gff_file, output_cds)
            click.echo(f"CDS FASTA written to {output_cds}")
        except (FileNotFoundError, Exception) as exc:
            click.echo(f"Error: {exc}", err=True)
            raise SystemExit(1)
