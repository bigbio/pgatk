SMALL example for database and peptides to map
==============================================

- Protein database: minimal_gencode.v25.pc_translations.fa => This is the minimal fasta file that produces the same output

- Small peptide output: Testfile_small.txt => Peptide list.

- GeneCode Annotations in GTF: gencode.v25.annotation.gtf => GTF Annotations.

Steps to generate minimal fasta files:
--------------------------------------
1. From the expected "_out.gtf" file, extract the required "gene_ids" and combine them using pipe("|") as delimiter using the following command (replace 'xyz_out.gtf' with original file)

`awk '{ print $10 }' xyz_out.gtf | sort | uniq | sed s/'[";]'//g | paste -d\| -s -`

Sample output of the above cmd: `ENSG00000142798|ENSG00000011009|ENSG00000163209`

2. Now use the above generated pattern to grep only the required lines from the input fasta file
Example:
`egrep -A5 -w 'ENSG00000142798|ENSG00000011009|ENSG00000163209' original.fa > minimal.fa`

(note: If you are sure that the sequence is in only one line, then '-A1' is actually enough)

3. Check the 'minimal.fa' and clean up the unwanted lines (This is very much necessary to make sure file is in valid format)
