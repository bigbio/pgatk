.. _pepgenome

PepGenome: Mapping Peptidoforms to Genome Coordinates
------------------------------------

In proteogenomic analyses it is essential to know the loci (genome position) giving rise to **peptides** in order to improve genomic annotation and the functional characterization of protein products in their biological context. With next-generation sequencing of DNA and RNA for each sample studied by proteomic mass spectrometry integration and visualisation in a common coordinate system, i.e. the genome, is vital for systems biology. To facilitate this type of integration not only the genomic locations of modified peptides but specifically the genomic loci of associated with these modifications is required.


.. note:: The **PepGenome** tool quickly and efficiently identify genomic loci of peptides and post-translational modifications and couple these mappings with associated quantitative values. Using reference gene annotation (GTF files) and an associated transcript translations (Protein fasta files) our tool identifies the genomic loci of peptides given as input and generates output in different formats borrowed from genomics and transcriptomics which can be loaded in various genome browsers such as `UCSC Genome Browser https://genome.ucsc.edu/`_, `Ensembl Genome Browser http://www.ensembl.org/index.html`_.



Input format
~~~~~~~~~~~~

The input format required by PepGenome is a tab delimited file with four columns.

<table border="0" width="100%"><thead><tr><th scope="col">Column</th><th scope="col">Column header</th><th scope="col">Description</th></tr></thead><tbody><tr><td>1</td><td>Sample</td><td>Name of sample or experiment|<tr><td>2</td><td>Peptide</td><td>Peptide sequence with PSI-MS nodification names in round brackets following the mpdified amino acid, e.g. PEPT(Phopsho)IDE for a phosphorylated threonine</td></tr><tr><td>3</td><td>PSMs</td><td>Number of peptide-spectrum matches (PSMs) for the given peptide</td></tr><tr><td>4</td><td>Quant</td><td>Quantitative value for the given peptide in the given sample</td></tr></tbody></table>

.. note:: In addition the tool support mzTab and mzIdentML File format input.

Output formats
~~~~~~~~~~~~

BED
~~~~~~~~~~

This format contains the genomic loci for peptides, the exon-structure, the peptide sequence, as well as a colour code for uniqueness of peptides within the genome.

<table align="left" border="0" width="100%">
	<thead>
		<tr>
			<th scope="col" width="20%">Colour</th>
			<th scope="col" width="80%">Description</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td bgcolor="#F00000"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/uniquetranscript.svg" height="25px"></td>
			<td>Peptide is unique to single gene AND single transcript</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#000000"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/uniquegene.svg" height="25px"></td>
			<td>Peptide is unique to single gene BUT shared between multiple transcripts</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#808080"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/notunique.svg" height="25px"></td>
			<td>Peptide is shared between multiple genes</td>
		</tr>
	</tbody>
</table>

PTMBED
~~~~~~~~~~~

Like BED but containing the location of the post-translational modification on the genome. Thick parts of the peptide blocks indicate the position of the post-translational modification on a single amino acid (short thick block) while longer blocks indicate the occurrence of the first and last post-translational modification and residues in between. In the PTMBED the colour code is changed to indicate the type of modification.

<table border="0" width="100%">
	<thead>
		<tr>
			<th scope="col" width="20%">Colour</th>
			<th scope="col" width="80%">Post-translational Modification</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td bgcolor="#FF3333"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/phospho.svg" height="25px"></td>
			<td>Phosphorylation (phospho)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#CC6600"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/acetyl.svg" height="25px"></td>
			<td>Acetylation (acetyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#FF9933"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/amidated.svg" height="25px"></td>
			<td>Amidation (amidated)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#CCCC00"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/oxidation.svg" height="25px"></td>
			<td>Oxidation (oxidation)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#00CC00"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/methyl.svg" height="25px"></td>
			<td>Methylation (methyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#33FF33"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/glygly.svg" height="25px"></td>
			<td>Ubiquitinylation (glygly; gg)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#00CCCC"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/sulfo.svg" height="25px"></td>
			<td>Sulfation (sulfo)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#3399FF"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/palmitoyl.svg" height="25px"></td>
			<td>Palmitoylation (palmitoyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#0000CC"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/formyl.svg" height="25px"></td>
			<td>Formylation (formyl)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#3333FF"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/deamidated.svg" height="25px"></td>
			<td>Deamidation (deamidated)</td>
		</tr>
		<tr></tr>
		<tr>
			<td bgcolor="#FF3399"><img src="http://ngs.sanger.ac.uk/production/proteogenomics/trackhubs_files/color/any.svg" height="25px"></td>
			<td>Any other post-translational modification</td>
		</tr>
	</tbody>
</table>

GTF
~~~~~~~~~

This output format contains besides the genomic loci the annotated information for the genes giving rise to each peptide sequence including status and biotype. For each mapped peptide the sample, number of peptide-spectrum matches and associated quantitative value as tags.

GCT
~~~~~~~~~~~
In this format the peptide sequences are combines with the Ensembl gene identifier. It contains the genomic loci for each peptide as well as the quantitative values for each peptide in different samples as a matrix.

Usage
~~~~~~~~~~~

**Required arguments**:
<table border="0" widht="100%"><tbody><tr><td width="20%">
<pre>-fasta TRANSL</pre>
</td><td>Filepath for file containing protein sequences in FASTA format</td></tr><tr><td>
<pre>-gtf ANNO</pre>
</td><td width="80%">Gene annotation with coding sequences (CDS) in GTF format</td></tr><tr><td>
<pre>-in *.tsv</pre>
</td><td>Path to single input file or comma separated list of paths to input files containing peptides to be mapped with associated number of peptide to spectrum matches, sample name and quantitative value (see input file format)</td></tr></tbody></table>

**Optional arguments**:

<table border="0" width="100%"><tbody><tr><td width="20%">
<pre>-format OUTF</pre>
</td><td width="80%">Set output format GTF, GCT, BED, PTMBED or ALL. Comma separated combination possible. Default = ALL</td></tr><tr><td>
<pre>-merge TRUE/FALSE</pre>
</td><td>Set TRUE to merge output of multiple input files (output will be named after last input file *_merged). Default = FALSE</td></tr><tr><td><pre>-source SRC</pre></td><td>Set TRUE to merge output of multiple input files (output will be named after last input file *_merged). Default = FALSE</td></tr><tr><td>
<pre>-mm NUM</pre>
</td><td>Number of mismatches allowed in mapping (0, 1 or 2). DEFAULT = 0</td></tr><tr><td>
<pre>-mmmode TRUE/FALSE</pre>
</td><td>Set TRUE to restrict number of mismatch in kmer to 1. DEFAULT = FALSE</td></tr><tr><td>
<pre>-species SPECIES</pre></td><td>Set species using common or scientific name or taxonomy ID. Default is Human (Homo sapiens, 9606).</td></tr>
<tr><td><pre>-chr NUM</pre></td><td>Export chr prefix Allowed 0, 1.  (DEFAULT = 0)</td></tr>
</tbody></table>

Table of supported species
~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+-------------------+-----------+
| Common name       | Scientific name   | Taxon ID  |
+===================+===================+===========+
|C.intestinalis     |Ciona intestinalis |7719       |
+-------------------+-------------------+-----------+
|Cat                |Felis catus        |9685       |
+-------------------+-------------------+-----------+
|Chicken            |Gallus gallus      |9031       |
+-------------------+-------------------+-----------+
|Chimpanzee         |Pan troglodytes    |9598       |
+-------------------+-------------------+-----------+
|Cow                |Bos taurus         |9913       |
+-------------------+-------------------+-----------+
|Dog                |Canis lupus        |9615      |
+-------------------+-------------------+-----------+
|Gorilla            |Gorilla            |9595       |
+-------------------+-------------------+-----------+
|Horse              |Equus caballus     |9796       |
+-------------------+-------------------+-----------+
|Human              |Homo sapiens       |9606       |
+-------------------+-------------------+-----------+
|Macaque            |Macaca mulatta     |9544       |
+-------------------+-------------------+-----------+
|Marmoset           |Callithrix jacchus |9483       |
+-------------------+-------------------+-----------+
|Medaka             |Oryzias latipes    |8090       |
+-------------------+-------------------+-----------+
|Mouse              |Mus musculus       |10090      |
+-------------------+-------------------+-----------+
|Olive baboon       |Papio anubis       |9555       |
+-------------------+-------------------+-----------+
|Opossum            |Monodelphis        |13616      |
+-------------------+-------------------+-----------+
|Orangutan          |Pongo abelii       |9601       |
+-------------------+-------------------+-----------+
|Pig                |Sus scrofa         |9823       |
+-------------------+-------------------+-----------+
|Platypus           |Ornithorhynchus    |9258       |
+-------------------+-------------------+-----------+
|Rabbit             |Oryctolagus        |9986       |
+-------------------+-------------------+-----------+
|Rat                |Rattus norvegicus  |10116      |
+-------------------+-------------------+-----------+
|Sheep              |Ovis aries         |9940       |
+-------------------+-------------------+-----------+
|Tetraodon          |Tetraodon          |99883      |
+-------------------+-------------------+-----------+
|Turkey             |Meleagris          |9103       |
+-------------------+-------------------+-----------+
|Vervet-AGM         |Chlorocebus        |60711      |
+-------------------+-------------------+-----------+
|Zebra Finch        |Taeniopygia        |59729      |
+-------------------+-------------------+-----------+


