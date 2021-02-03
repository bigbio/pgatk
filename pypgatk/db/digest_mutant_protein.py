import getopt
import re
import sys
from collections import OrderedDict

from Bio import SeqIO


def TRYPSIN(proseq, miss_cleavage):
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
        for j, pep in enumerate(peptides):
            if j + i < len(peptides):
                peptide = ''.join([x for x in (peptides[j:j + i + 1])])
                peptides_with_miss_cleavage.append(peptide)

    peptides.extend(peptides_with_miss_cleavage)

    return peptides


min = 7  # default minimum peptide length
max = 40  # default maximum peptide length
header_prefix = "Mutation"
cleavageMiss = 0  # number of misscleavage,default no misscleavage

if len(sys.argv[1:]) <= 1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command!")
    print(
        "Example: python3 digest_mutant_protein.py --input mutproteins.1.fa, mutproteins.2.fa --fa knownproteins.fa --prefix Mutation --output mutpeptides.fa --min_len 7 --max_len 40 --miss 0")
else:
    options, remainder = getopt.getopt(sys.argv[1:], '',
                                       ['input=', 'fa=', 'output=', 'prefix=', 'min_len', 'max_len', 'miss'])
    for opt, arg in options:
        if opt == '--input':
            input_file = arg
        elif opt == '--output':
            output_file = arg
        elif opt == '--fa':
            fa_file = arg
        elif opt == '--prefix':
            header_prefix = arg
        elif opt == '--min_len':
            min = int(arg)
        elif opt == '--max_len':
            max = int(arg)
        elif opt == '--miss':
            cleavageMiss = int(arg)
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt)
            sys.exit()

handle1 = SeqIO.parse(fa_file, 'fasta')  # canonical protein sequences
peptidome = {}

for record in handle1:
    aa_seq = record.seq
    peptide_list = TRYPSIN(str(aa_seq), cleavageMiss)
    for peptide in peptide_list:
        if len(peptide) in range(min, max + 1):
            if peptide not in peptidome:
                peptidome[peptide.replace("I", "L")] = 1

print("known peptides number", len(peptidome))
handle1.close()

filelist = input_file.split(",")
handle_list = []
for f in filelist:
    handle_list.append(SeqIO.parse(f, 'fasta'))  # mutproteins.fasta

output = open(output_file, 'w')

var_peptidome = OrderedDict()

for h in handle_list:
    for record in h:
        proseq = record.seq
        descrip = record.description
        if len(proseq) >= min:
            peptide_list = TRYPSIN(str(proseq), cleavageMiss)
            for peptide in peptide_list:
                if len(peptide) in range(min, max + 1):
                    peptide1 = peptide.replace("I", "L")
                    if peptide1 not in peptidome:  # check if the peptides are in canonical protein sequences
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
    print("file process done")

print("mut peptide numbers", len(var_peptidome))
for pep in var_peptidome.keys():
    acc = ";".join(set(var_peptidome[pep]))
    output.write(">%s\n%s\n" % (acc, pep))

output.close()
