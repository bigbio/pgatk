import click

from pgatk.gnomad.data_downloader import GnomadVcfDownloader, ALL_CHROMOSOMES


@click.command("gnomad-vcf-downloader", short_help="Download gnomAD per-chromosome VCF files in parallel")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option(
    "--version",
    default="3.1.2",
    show_default=True,
    help="gnomAD release version (e.g. 3.1.2 for genomes, 4.1.1 for exomes).",
)
@click.option(
    "--dataset",
    default="genomes",
    show_default=True,
    type=click.Choice(["genomes", "exomes"], case_sensitive=False),
    help="gnomAD dataset type. Use 'genomes' for v3.x and 'exomes' for v4.x.",
)
@click.option(
    "--chromosomes",
    default=None,
    help=(
        "Comma-separated chromosome list (e.g. chr1,chr22,chrX). "
        "Defaults to all autosomes + chrX + chrY."
    ),
)
@click.option(
    "--workers",
    default=4,
    show_default=True,
    type=int,
    help="Number of parallel download threads.",
)
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.pass_context
def gnomad_vcf_downloader(ctx, output_dir, version, dataset, chromosomes, workers, force):
    """Download gnomAD per-chromosome VCF (BGZF) and tabix index files in parallel.

    Downloads one .vcf.bgz + .vcf.bgz.tbi pair per chromosome from the
    gnomAD Google Cloud Storage bucket.  By default all 24 chromosomes are
    fetched concurrently using --workers threads.

    \b
    Example — all chromosomes (genomes, v3.1.2):
      pgatk gnomad-vcf-downloader -o gnomad_data --version 3.1.2 --dataset genomes --workers 8

    \b
    Example — single chromosome:
      pgatk gnomad-vcf-downloader -o gnomad_data --chromosomes chr22
    """
    chrom_list = (
        [c.strip() for c in chromosomes.split(",") if c.strip()]
        if chromosomes
        else ALL_CHROMOSOMES
    )

    n_files = len(chrom_list) * 2  # VCF + TBI per chromosome
    click.echo(
        f"Downloading gnomAD v{version} — {len(chrom_list)} chromosome(s), "
        f"{n_files} files, {workers} parallel thread(s) ..."
    )

    completed_count = [0]

    def _progress(completed, total, name):
        completed_count[0] = completed
        click.echo(f"  [{completed}/{total}] {name}", err=False)

    downloader = GnomadVcfDownloader(
        output_dir=output_dir,
        version=version,
        dataset=dataset,
        chromosomes=chrom_list,
        workers=workers,
    )
    downloaded = downloader.download_all(force=force, progress_callback=_progress)
    click.echo(f"Done. {len(downloaded)}/{n_files} file(s) ready in {output_dir}")
