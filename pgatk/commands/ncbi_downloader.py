import click

from pgatk.clinvar.data_downloader import NcbiDataDownloader
from pgatk.config.registry import load_config


@click.command("ncbi-downloader", short_help="Download NCBI RefSeq and ClinVar reference files")
@click.option("-c", "--config_file", help="Configuration YAML file (optional)")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.pass_context
def ncbi_downloader(ctx, config_file, output_dir, force):
    """Download NCBI RefSeq GTF, protein FASTA, assembly report, and ClinVar VCF.

    Files are downloaded to the specified output directory. Existing files
    are skipped unless --force is used.
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
