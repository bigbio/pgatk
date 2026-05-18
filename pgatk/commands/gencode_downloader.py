import os

import click

from pgatk.gnomad.data_downloader import GencodeDownloader


@click.command("gencode-downloader", short_help="Download GENCODE GTF and genome FASTA")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option(
    "--release",
    default=44,
    show_default=True,
    type=int,
    help="GENCODE release number (e.g. 44 → gencode.v44.*).",
)
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.option(
    "--generate-transcripts",
    "generate_transcripts",
    is_flag=True,
    default=False,
    help=(
        "Run gffread after download to generate transcripts.fa with CDS= headers. "
        "Required for vcf-to-proteindb (gnomAD). Requires gffread in PATH "
        "(conda install -c bioconda gffread)."
    ),
)
@click.pass_context
def gencode_downloader(ctx, output_dir, release, force, generate_transcripts):
    """Download GENCODE human GTF annotation and primary-assembly genome FASTA.

    Downloads two files from the EBI GENCODE FTP:

    \b
      gencode.v<RELEASE>.annotation.gtf          (annotation)
      GRCh38.primary_assembly.genome.fa          (genomic DNA)

    Use --generate-transcripts to run gffread and produce transcripts.fa with
    CDS= coordinate headers.  This file is required by vcf-to-proteindb
    (--input_fasta) when processing gnomAD VCFs.  Without it the pipeline
    falls back to 3-frame exon translation, which is incorrect for transcripts
    with a 5' UTR.

    \b
    Example — download only:
      pgatk gencode-downloader -o gencode_data --release 44

    \b
    Example — download + generate transcripts:
      pgatk gencode-downloader -o gencode_data --release 44 --generate-transcripts
    """
    downloader = GencodeDownloader(output_dir=output_dir, release=release)
    click.echo(f"Downloading GENCODE release {release} to {output_dir} ...")
    downloaded = downloader.download_all(force=force)
    click.echo(f"Downloaded {len(downloaded)} file(s):")
    for path in downloaded:
        click.echo(f"  {path}")

    if generate_transcripts:
        genome_fna = os.path.join(output_dir, "GRCh38.primary_assembly.genome.fa")
        gtf_file = os.path.join(output_dir, f"gencode.v{release}.annotation.gtf")
        output_fasta = os.path.join(output_dir, "transcripts.fa")
        try:
            click.echo("Running gffread to generate transcripts.fa with CDS= headers ...")
            GencodeDownloader.generate_transcripts(genome_fna, gtf_file, output_fasta)
            click.echo(f"Transcripts written to {output_fasta}")
        except FileNotFoundError as exc:
            click.echo(f"Error: {exc}", err=True)
            raise SystemExit(1)
