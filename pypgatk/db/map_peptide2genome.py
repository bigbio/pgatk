"""
    this script will map proteomics identified peptides to genome and output an peptide gff3 file with coordinates.
    The input peptide table should contain minimum two columns, peptide sequence and Ensembl protein ID.
    Three additional files are required:
    1. Ensembl protein db (same one used in database search)
    2. Ensembl GTF annotation
    3. IDmap file with three columns, first column Ensembl Gene ID, second column Ensembl transcript ID and third Ensembl protein ID.
"""

import sys
import os
import getopt
import numpy as np
import re
from Bio import SeqIO


class EXON(object):
    def __init__(self, number=0, gene=None, variant=None, chr=None, strand=None, start=0, end=0, length=0,
                 trans_start=0, trans_end=0):
        self.gene = gene
        self.variant = variant
        self.number = number
        self.start = start  # chromosome start coordinate
        self.end = end  # chromosome end coordinate
        self.strand = strand
        self.chr = chr
        self.trans_start = trans_start
        self.trans_end = trans_end

    def length(self):
        self.length = self.trans_end - self.trans_start + 1


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
    for i in range(len(exon_object_list)):
        exon = exon_object_list[i]
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
    with open(infile, "r") as infile_object:
        for line in infile_object:
            if line[0] != "#":  # skip lines commmented out
                row = line.strip().split("\t")
                if row[2] == "CDS":
                    attri_list = row[8].split(";")
                    transID = ""
                    exon = EXON(start=int(row[3]), end=int(row[4]), chr=row[0], strand=row[6])
                    for attri in attri_list:
                        if "transcript_id" in attri:
                            transID = attri.strip().replace("transcript_id ", "").replace('\"', "")

                    if transID not in dic:
                        dic[transID] = [exon]
                    else:
                        dic[transID].append(exon)
    return dic


################  Comand-line arguments ################
# default peptide sequence in first column, protein ID in second column
pep_col = 1
prot_col = 2

if len(sys.argv[1:]) <= 1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual")
    print(
        "Example: python map_peptide2genome.py --input peptideTable.txt --gtf Homo_sapiens.GRCh37.90.gtf --fasta Homo_sapiens.GRCh37.90.pep.all.fa  --IDmap Ensembl90_IDmap.txt --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:], '',
                                       ['input=', 'gtf=', 'fasta=', 'IDmap=', 'output=', 'pep_col=', 'prot_col='])
    for opt, arg in options:
        if opt == '--input':
            input_file = arg
        elif opt == '--gtf':
            gtf_file = arg
        elif opt == '--fasta':
            fasta_file = arg
        elif opt == '--IDmap':
            IDmap_file = arg
        elif opt == '--output':
            output_file = arg
        elif opt == '--pep_col':
            pep_col = int(arg)
        elif opt == '--prot_col':
            prot_col = int(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt)
            sys.exit()

print("reading GTF input file")
feature_dic = parse_gtf(gtf_file)
print("number of unique transcripts in GTF file", len(feature_dic))

IDlist_input = open(IDmap_file, "r")
id_dic = {}
for line in IDlist_input:
    row = line.strip().split("\t")
    if len(row) == 3:
        enst = row[1]
        ensp = row[2]
        if ensp not in id_dic:
            id_dic[ensp] = enst
IDlist_input.close()
print("number of unique ENSP IDs in ID table", len(id_dic))

pep_dic = {}
input = open(input_file,
             'r')  # peptide table with two columns, peptide sequence in first column, protein ID in second column

input.readline()
for line in input:
    row = line.strip().split("\t")
    pep = re.sub("[\W\d]", "", row[pep_col - 1].strip())
    acc = row[prot_col - 1].split(";")[0]  # in case there are multiple IDs
    if pep not in pep_dic:
        pep_dic[pep] = acc

input.close()

output = open(output_file, 'w')

seq_dic = SeqIO.index(fasta_file, 'fasta')
print("number of unique protein sequences in fasta file", len(seq_dic))

non_mapped_pep = 0

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

    # print pep_trans_start,pep_trans_end
    pep_chr, pep_strand, pep_chr_start, pep_chr_end, pep_start_exon, pep_end_exon = get_pep_cor(exons, pep_trans_start,
                                                                                                pep_trans_end)

    # handle exceptions
    if pep_chr_start > pep_chr_end:
        non_mapped_pep += 1
        # print peptide,ensp,enst
        continue
    if pep_chr_start <= 0:
        non_mapped_pep += 1
        # print peptide,ensp,enst,pep_trans_start,pep_trans_end
        continue

    # print pep_chr_start,pep_chr_end
    # print pep_start_exon,pep_end_exon
    pep_chr = "chr" + pep_chr.replace("MT", "M")
    if pep_start_exon == pep_end_exon:  # if peptide map to one exon
        gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".", "ID=" + peptide]
        gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, pep_chr_end, ".", pep_strand, "0", "Parent=" + peptide]
        output.write("\t".join(map(str, gff_format_line1)) + "\n")
        output.write("\t".join(map(str, gff_format_line2)) + "\n")
    elif abs(pep_start_exon - pep_end_exon) == 1:  # if it is a splice junction peptide spanning over two exons
        if pep_strand == "+":
            gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                "ID=" + peptide]
            gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_start_exon - 1].end, ".", pep_strand,
                                "0", "Parent=" + peptide]
            gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_end_exon - 1].start, pep_chr_end, ".", pep_strand, ".",
                                "Parent=" + peptide]
        else:
            gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                "ID=" + peptide]
            gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_end_exon - 1].end, ".", pep_strand, "0",
                                "Parent=" + peptide]
            gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_start_exon - 1].start, pep_chr_end, ".", pep_strand,
                                ".", "Parent=" + peptide]

        output.write("\t".join(map(str, gff_format_line1)) + "\n")
        output.write("\t".join(map(str, gff_format_line2)) + "\n")
        output.write("\t".join(map(str, gff_format_line3)) + "\n")
    elif abs(pep_start_exon - pep_end_exon) > 1:  # if peptide span multiple exons,rare case!
        if pep_strand == "+":
            gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                "ID=" + peptide]
            gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_start_exon - 1].end, ".", pep_strand,
                                "0", "Parent=" + peptide]
            gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_end_exon - 1].start, pep_chr_end, ".", pep_strand, ".",
                                "Parent=" + peptide]
        else:
            gff_format_line1 = [pep_chr, "MS", "mRNA", pep_chr_start, pep_chr_end, ".", pep_strand, ".",
                                "ID=" + peptide]
            gff_format_line2 = [pep_chr, "MS", "CDS", pep_chr_start, exons[pep_end_exon - 1].end, ".", pep_strand, "0",
                                "Parent=" + peptide]
            gff_format_line3 = [pep_chr, "MS", "CDS", exons[pep_start_exon - 1].start, pep_chr_end, ".", pep_strand,
                                ".", "Parent=" + peptide]

        output.write("\t".join(map(str, gff_format_line1)) + "\n")
        output.write("\t".join(map(str, gff_format_line2)) + "\n")
        for k in range(min(pep_start_exon, pep_end_exon) + 1, max(pep_start_exon, pep_end_exon)):
            gff_format_line = [pep_chr, "MS", "CDS", exons[k - 1].start, exons[k - 1].end, ".", pep_strand, ".",
                               "Parent=" + peptide]
            output.write("\t".join(map(str, gff_format_line)) + "\n")

        output.write("\t".join(map(str, gff_format_line3)) + "\n")

output.close()
print("total number of unique peptides", len(pep_dic))
print("total number of unmapped peptides", non_mapped_pep)
