import logging
import re
from collections import OrderedDict

from Bio import SeqIO

log = logging.getLogger(__name__)


def trypsin_cleavage(proseq: str, miss_cleavage: int):
    peptides = []
    peptide = ''
    for c, aa in enumerate(proseq):
        peptide += aa
        next_aa = ''
        try:
            next_aa = proseq[c + 1]
        except IndexError:
            pass

        if aa in ['K', 'R'] and next_aa != 'P':  # for trypsin peptides
            if len(peptide) > 0:
                peptides.append(peptide)
            peptide = ''
            continue

    if len(peptide) > 0:
        peptides.append(peptide)

    peptides_with_miss_cleavage = []
    for i in range(1, miss_cleavage + 1):
        for j in range(len(peptides)):
            if j + i < len(peptides):
                peptide = ''.join([x for x in (peptides[j:j + i + 1])])
                peptides_with_miss_cleavage.append(peptide)

    peptides.extend(peptides_with_miss_cleavage)

    return peptides


def digest_mutant_proteins(
    input_file: str,
    fasta_file: str,
    output_file: str,
    header_prefix: str = "Mutation",
    min_length: int = 7,
    max_length: int = 40,
    missed_cleavages: int = 0,
) -> None:
    """
    Digest mutant proteins and filter peptides against a canonical proteome.

    Reads a reference FASTA to build a canonical peptidome, then digests mutant
    protein sequences from the input file(s) and writes only those peptides that
    are not found in the canonical peptidome.

    :param input_file: Comma-separated list of mutant protein FASTA file paths.
    :param fasta_file: Reference canonical protein FASTA file path.
    :param output_file: Output file path for unique variant peptides.
    :param header_prefix: Prefix for the output FASTA headers (default: "Mutation").
    :param min_length: Minimum peptide length to retain (default: 7).
    :param max_length: Maximum peptide length to retain (default: 40).
    :param missed_cleavages: Number of missed cleavages for trypsin digestion (default: 0).
    """
    # Build canonical peptidome from reference FASTA
    handle1 = SeqIO.parse(fasta_file, 'fasta')
    peptidome = {}

    for record in handle1:
        aa_seq = record.seq
        peptide_list = trypsin_cleavage(str(aa_seq), missed_cleavages)
        for peptide in peptide_list:
            if min_length <= len(peptide) <= max_length:
                if peptide not in peptidome:
                    peptidome[peptide.replace("I", "L")] = 1

    log.info("Known peptides number: %d", len(peptidome))
    handle1.close()

    # Parse mutant protein input files
    filelist = input_file.split(",")
    handle_list = []
    for f in filelist:
        handle_list.append(SeqIO.parse(f.strip(), 'fasta'))

    var_peptidome = OrderedDict()

    for h in handle_list:
        for record in h:
            proseq = record.seq
            descrip = record.description
            if len(proseq) >= min_length:
                peptide_list = trypsin_cleavage(str(proseq), missed_cleavages)
                for peptide in peptide_list:
                    if min_length <= len(peptide) <= max_length:
                        peptide1 = peptide.replace("I", "L")
                        if peptide1 not in peptidome:
                            des_list = descrip.split(":")
                            des_list[0] = header_prefix
                            mut_type = des_list[-1]
                            snp = des_list[-2]
                            if "Missense" in mut_type:
                                try:
                                    mut_pos = int(re.findall(r'\d+', snp)[0])
                                    index = str(proseq).index(peptide)
                                    pos = mut_pos - index
                                    des_list.append(str(pos))
                                except IndexError:
                                    continue

                            new_description = ":".join(des_list).replace("*", "-")
                            if peptide not in var_peptidome:
                                var_peptidome[peptide] = [new_description]
                            else:
                                var_peptidome[peptide].append(new_description)
        h.close()
        log.info("File processing done")

    log.info("Mutant peptide numbers: %d", len(var_peptidome))
    with open(output_file, 'w', encoding='utf-8') as output:
        for pep in var_peptidome.keys():
            acc = ";".join(set(var_peptidome[pep]))
            output.write(">%s\n%s\n" % (acc, pep))


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4:
        print("Usage: python digest_mutant_protein.py <input_file> <fasta_file> <output_file>")
        sys.exit(1)
    digest_mutant_proteins(sys.argv[1], sys.argv[2], sys.argv[3])
