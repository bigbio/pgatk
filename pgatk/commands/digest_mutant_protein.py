import logging

import click

from pgatk.db.digest_mutant_protein import digest_mutant_proteins

log = logging.getLogger(__name__)


@click.command("digest-mutant-protein", short_help="Digest mutant proteins and filter against canonical proteome")
@click.option("-i", "--input", "input_file", required=True,
              help="Input mutant protein FASTA file(s), comma-separated")
@click.option("-f", "--fasta", "fasta_file", required=True,
              help="Reference canonical protein FASTA")
@click.option("-o", "--output", "output_file", required=True,
              help="Output file for unique variant peptides")
@click.option("--prefix", default="Mutation",
              help="Header prefix for output FASTA entries (default: Mutation)")
@click.option("--min-len", default=7, type=int,
              help="Minimum peptide length (default: 7)")
@click.option("--max-len", default=40, type=int,
              help="Maximum peptide length (default: 40)")
@click.option("--missed-cleavages", default=0, type=int,
              help="Number of missed cleavages (default: 0)")
@click.pass_context
def digest_mutant_protein(ctx, input_file, fasta_file, output_file, prefix, min_len, max_len, missed_cleavages):
    """Digest mutant proteins and filter peptides against a canonical proteome.

    Reads a reference FASTA to build a canonical peptidome, then digests mutant
    protein sequences and writes only variant-specific peptides to the output.
    """
    digest_mutant_proteins(
        input_file=input_file,
        fasta_file=fasta_file,
        output_file=output_file,
        header_prefix=prefix,
        min_length=min_len,
        max_length=max_len,
        missed_cleavages=missed_cleavages,
    )
