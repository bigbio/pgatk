import click

from pgatk.clinvar.clinvar_service import ClinVarService


@click.command("clinvar-to-proteindb", short_help="Generate protein database from ClinVar VCF + RefSeq GFF3")
@click.option("-c", "--config_file", help="Configuration YAML file (optional, uses bundled defaults)")
@click.option("-v", "--vcf", required=True, help="ClinVar VCF file path")
@click.option("-g", "--gff", required=True, help="NCBI RefSeq GFF3 annotation file path")
@click.option("-f", "--fasta", required=True, help="RefSeq transcript nucleotide FASTA file path")
@click.option("-a", "--assembly-report", required=True, help="NCBI assembly report file path")
@click.option("-o", "--output", required=True, help="Output protein FASTA file path")
@click.pass_context
def clinvar_to_proteindb(ctx, config_file, vcf, gff, fasta, assembly_report, output):
    """Generate a variant protein database from ClinVar VCF and NCBI RefSeq GFF3.

    This command does NOT require VEP annotations. It uses BedTools interval
    overlap to find transcripts affected by each ClinVar variant, then applies
    the variant and translates to protein.

    Use ``ncbi-downloader --generate-transcripts`` to produce the GFF3 annotation
    and transcript FASTA required by this command.
    """
    service = ClinVarService(
        vcf_file=vcf,
        gff_file=gff,
        fasta_file=fasta,
        assembly_report=assembly_report,
        output_file=output,
        config_file=config_file,
    )
    service.run()
