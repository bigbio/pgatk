import logging

import click

from pgatk.db.map_peptide2genome import map_peptides_to_genome

log = logging.getLogger(__name__)


@click.command("map-peptide2genome", short_help="Map peptides to genomic coordinates (GFF3 output)")
@click.option("-i", "--input", "input_file", required=True, help="Input peptide identification TSV")
@click.option("-g", "--gtf", "gtf_file", required=True, help="GTF gene annotation file")
@click.option("-f", "--fasta", "fasta_file", required=True, help="Protein FASTA file")
@click.option("-m", "--idmap", "idmap_file", required=True, help="Protein-to-transcript ID mapping file")
@click.option("-o", "--output", "output_file", required=True, help="Output GFF3 file")
@click.option("--pep-col", default=0, type=int, help="Peptide column index (0-based, default: 0)")
@click.option("--prot-col", default=1, type=int, help="Protein column index (0-based, default: 1)")
@click.pass_context
def map_peptide2genome(ctx, input_file, gtf_file, fasta_file, idmap_file, output_file, pep_col, prot_col):
    """Map proteomics-identified peptides to genomic coordinates.

    Reads a peptide identification table (TSV), an Ensembl GTF annotation, a
    protein FASTA, and a protein-to-transcript ID mapping file to produce a
    GFF3 file with peptide genomic coordinates.
    """
    map_peptides_to_genome(
        input_file=input_file,
        gtf_file=gtf_file,
        fasta_file=fasta_file,
        idmap_file=idmap_file,
        output_file=output_file,
        pep_col=pep_col,
        prot_col=prot_col,
    )
