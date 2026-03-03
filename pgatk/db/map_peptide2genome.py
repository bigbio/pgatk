"""
Map proteomics-identified peptides to genomic coordinates and output a GFF3 file.

The input peptide table should contain at minimum two columns: peptide sequence
and Ensembl protein ID.  Three additional files are required:

1. Ensembl protein database (same one used in the database search)
2. Ensembl GTF annotation
3. ID-mapping file with three columns: Ensembl Gene ID, Ensembl Transcript ID,
   and Ensembl Protein ID.
"""

import logging
import re

from Bio import SeqIO

log = logging.getLogger(__name__)


class EXON(object):
    def __init__(self, number=0, gene=None, variant=None, chromosome=None, strand=None, start=0, end=0, exon_length=0,
                 trans_start=0, trans_end=0):
        self.gene = gene
        self.variant = variant
        self.number = number
        self.start = start  # chromosome start coordinate
        self.end = end  # chromosome end coordinate
        self.strand = strand
        self.chr = chromosome
        self.trans_start = trans_start
        self.trans_end = trans_end
        self.exon_length = exon_length

    def calculate_length(self):
        self.exon_length = self.trans_end - self.trans_start + 1


def cal_trans_pos(exon_list):  # calculate transcript position of exon start & end, exon_list is a list of exon objects
    strand = exon_list[0].strand
    if strand == "+":
        new_exon_list = sorted(exon_list, key=lambda x: x.start)
    else:
        new_exon_list = sorted(exon_list, key=lambda x: x.start, reverse=True)

    sumExonlength = 0
    for exon in new_exon_list:
        exon_length = exon.end - exon.start + 1
        exon.trans_start = 1 + sumExonlength
        exon.trans_end = exon.trans_start + exon_length - 1
        sumExonlength += exon_length

    return new_exon_list


def get_pep_cor(exon_object_list, n1,
                n2):  # return peptide's chromosome start and end cor given peptide's trans_start (n1) and trans_end (n2)
    pep_chr = ""
    pep_strand = ""
    pep_chr_start = 0
    pep_chr_end = 0
    pep_start_exon = 0
    pep_end_exon = 0
    for i, exon in enumerate(exon_object_list):
        if exon.trans_end >= n1 >= exon.trans_start:
            pep_chr = exon.chr
            pep_strand = exon.strand
            pep_start_exon = i + 1
            if pep_strand == '+':
                pep_chr_start = exon.start + (n1 - exon.trans_start)
            else:
                pep_chr_end = exon.end - (n1 - exon.trans_start)

        if exon.trans_end >= n2 >= exon.trans_start:
            pep_chr = exon.chr
            pep_strand = exon.strand
            pep_end_exon = i + 1
            if pep_strand == '+':
                pep_chr_end = exon.start + (n2 - exon.trans_start)
            else:  # chr_cor of n2 is pep_chr_start
                pep_chr_start = exon.end - (n2 - exon.trans_start)

    return pep_chr, pep_strand, pep_chr_start, pep_chr_end, pep_start_exon, pep_end_exon


def parse_gtf(infile):
    dic = {}
    with open(infile, "r", encoding='utf-8') as infile_object:
        for line in infile_object:
            if line[0] != "#":  # skip lines commented out
                row = line.strip().split("\t")
                if row[2] == "CDS":
                    attri_list = row[8].split(";")
                    transID = ""
                    exon = EXON(start=int(row[3]), end=int(row[4]), chromosome=row[0], strand=row[6])
                    for attri in attri_list:
                        if "transcript_id" in attri:
                            transID = attri.strip().replace("transcript_id ", "").replace('\"', "")

                    if transID not in dic:
                        dic[transID] = [exon]
                    else:
                        dic[transID].append(exon)
    return dic


def map_peptides_to_genome(
    input_file: str,
    gtf_file: str,
    fasta_file: str,
    idmap_file: str,
    output_file: str,
    pep_col: int = 0,
    prot_col: int = 1,
) -> None:
    """Map identified peptides to genomic coordinates and write a GFF3 file.

    Parameters
    ----------
    input_file : str
        Path to the input peptide identification TSV (header row expected).
    gtf_file : str
        Path to the Ensembl GTF gene annotation file.
    fasta_file : str
        Path to the Ensembl protein FASTA file.
    idmap_file : str
        Path to the protein-to-transcript ID mapping file (Gene / Transcript / Protein columns).
    output_file : str
        Path to the output GFF3 file.
    pep_col : int
        Zero-based column index for the peptide sequence (default 0).
    prot_col : int
        Zero-based column index for the protein accession (default 1).
    """
    # Convert to 1-based internally (original script used 1-based columns)
    _pep_col = pep_col + 1
    _prot_col = prot_col + 1

    log.info("Reading GTF input file")
    feature_dic = parse_gtf(gtf_file)
    log.info("Number of unique transcripts in GTF file: %d", len(feature_dic))

    with open(idmap_file, "r", encoding='utf-8') as IDlist_input:
        id_dic = {}
        for line in IDlist_input:
            row = line.strip().split("\t")
            if len(row) == 3:
                enst = row[1]
                ensp = row[2]
                if ensp not in id_dic:
                    id_dic[ensp] = enst
    log.info("Number of unique ENSP IDs in ID table: %d", len(id_dic))

    pep_dic = {}
    with open(input_file, 'r', encoding='utf-8') as input_stream:
        # peptide table with two columns, peptide sequence in first column, protein ID in second column
        input_stream.readline()
        for line in input_stream:
            row = line.strip().split("\t")
            pep = re.sub(r"[\W\d]", "", row[_pep_col - 1].strip())
            acc = row[_prot_col - 1].split(";")[0]  # in case there are multiple IDs
            if pep not in pep_dic:
                pep_dic[pep] = acc

    seq_dic = SeqIO.index(fasta_file, 'fasta')
    log.info("Number of unique protein sequences in FASTA file: %d", len(seq_dic))

    non_mapped_pep = 0

    with open(output_file, 'w', encoding='utf-8') as output:
        for peptide, ensp in pep_dic.items():
            enst = id_dic[ensp]
            try:
                exons = feature_dic[enst]
            except KeyError:
                non_mapped_pep += 1
                continue

            aa_seq = str(seq_dic[ensp].seq)
            pep_index = aa_seq.index(peptide)

            pep_trans_start = 3 * pep_index + 1
            pep_trans_end = pep_trans_start + 3 * len(peptide) - 1

            exons = cal_trans_pos(exons)

            pep_chr, pep_strand, pep_chr_start, pep_chr_end, pep_start_exon, pep_end_exon = get_pep_cor(
                exons, pep_trans_start, pep_trans_end
            )

            # handle exceptions
            if pep_chr_start > pep_chr_end:
                non_mapped_pep += 1
                continue
            if pep_chr_start <= 0:
                non_mapped_pep += 1
                continue

            pep_chr = "chr" + pep_chr.replace("MT", "M")
            if pep_start_exon == pep_end_exon:  # if peptide maps to one exon
                gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                    "ID=" + peptide]
                gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, pep_chr_end, ".", pep_strand, "0",
                                    "Parent=" + peptide]
                output.write("\t".join(map(str, gff_format_line1)) + "\n")
                output.write("\t".join(map(str, gff_format_line2)) + "\n")
            elif abs(pep_start_exon - pep_end_exon) == 1:  # splice junction peptide spanning two exons
                if pep_strand == "+":
                    gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                        "ID=" + peptide]
                    gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_start_exon - 1].end, ".",
                                        pep_strand, "0", "Parent=" + peptide]
                    gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_end_exon - 1].start, pep_chr_end, ".",
                                        pep_strand, ".", "Parent=" + peptide]
                else:
                    gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                        "ID=" + peptide]
                    gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_end_exon - 1].end, ".",
                                        pep_strand, "0", "Parent=" + peptide]
                    gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_start_exon - 1].start, pep_chr_end, ".",
                                        pep_strand, ".", "Parent=" + peptide]

                output.write("\t".join(map(str, gff_format_line1)) + "\n")
                output.write("\t".join(map(str, gff_format_line2)) + "\n")
                output.write("\t".join(map(str, gff_format_line3)) + "\n")
            elif abs(pep_start_exon - pep_end_exon) > 1:  # peptide spans multiple exons, rare case
                if pep_strand == "+":
                    gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                        "ID=" + peptide]
                    gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_start_exon - 1].end, ".",
                                        pep_strand, "0", "Parent=" + peptide]
                    gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_end_exon - 1].start, pep_chr_end, ".",
                                        pep_strand, ".", "Parent=" + peptide]
                else:
                    gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                        "ID=" + peptide]
                    gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_end_exon - 1].end, ".",
                                        pep_strand, "0", "Parent=" + peptide]
                    gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_start_exon - 1].start, pep_chr_end, ".",
                                        pep_strand, ".", "Parent=" + peptide]

                output.write("\t".join(map(str, gff_format_line1)) + "\n")
                output.write("\t".join(map(str, gff_format_line2)) + "\n")
                for k in range(min(pep_start_exon, pep_end_exon) + 1, max(pep_start_exon, pep_end_exon)):
                    gff_format_line = [pep_chr, "MS", "CDS", exons[k - 1].start, exons[k - 1].end, ".", pep_strand,
                                       ".", "Parent=" + peptide]
                    output.write("\t".join(map(str, gff_format_line)) + "\n")

                output.write("\t".join(map(str, gff_format_line3)) + "\n")

    log.info("Total number of unique peptides: %d", len(pep_dic))
    log.info("Total number of unmapped peptides: %d", non_mapped_pep)


if __name__ == '__main__':
    import sys

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    if len(sys.argv) < 6:
        print("Usage: python map_peptide2genome.py <input> <gtf> <fasta> <idmap> <output> [pep_col] [prot_col]")
        sys.exit(1)

    _pep_col = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    _prot_col = int(sys.argv[7]) if len(sys.argv) > 7 else 1

    map_peptides_to_genome(
        input_file=sys.argv[1],
        gtf_file=sys.argv[2],
        fasta_file=sys.argv[3],
        idmap_file=sys.argv[4],
        output_file=sys.argv[5],
        pep_col=_pep_col,
        prot_col=_prot_col,
    )
