.. _pepgenome


PepGenome: Mapping Peptidoforms to Genome Coordinates
===========================

In proteogenomic analyses it is essential to know the loci (genome position) giving rise to **peptides** in order to improve genomic annotation and the functional characterization of protein products in their biological context. With next-generation sequencing of DNA and RNA for each sample studied by proteomic mass spectrometry integration and visualisation in a common coordinate system, i.e. the genome, is vital for systems biology. To facilitate this type of integration not only the genomic locations of modified peptides but specifically the genomic loci of associated with these modifications is required.


.. note:: The **PepGenome** tool quickly and efficiently identify genomic loci of peptides and post-translational modifications and couple these mappings with associated quantitative values. Using reference gene annotation (GTF files) and an associated transcript translations (Protein fasta files) our tool identifies the genomic loci of peptides given as input and generates output in different formats borrowed from genomics and transcriptomics which can be loaded in various genome browsers such as `UCSC Genome Browser <https://genome.ucsc.edu/>`_, `Ensembl Genome Browser <http://www.ensembl.org/index.html>`_.


Input format (Tab delimited)
------------------------

The input format required by PepGenome is a tab delimited file with four columns. It can contains the following extensions ```.tsv```, ```.pogo``` or ```.txt```.

+-------------------+------------------+-----------------------------------------+
| Column            | Column header    | Description                             |
+===================+==================+=========================================+
|1                  |Sample            |Name of sample or experiment             |
+-------------------+------------------+-----------------------------------------+
|2                  |Peptide           |Peptide sequence *                       |
+-------------------+------------------+-----------------------------------------+
|3                  |PSMs              |Number (PSMs)                            |
+-------------------+------------------+-----------------------------------------+
|4                  |Quant             |Quantitative value in the given sample   |
+-------------------+------------------+-----------------------------------------+

.. hint:: *Peptide sequence with PSI-MS modification names in round brackets following the mpdified amino acid, e.g. PEPT(Phopsho)IDE for a phosphorylated threonine.


.. note:: In addition the tool support mzTab and mzIdentML File format input.

Output formats
------------------------

BED Files
~~~~~~~~~~

This format contains the genomic loci for peptides, the exon-structure, the peptide sequence, as well as a colour code for uniqueness of peptides within the genome. Read here :ref:`bed`


Additional Files
~~~~~~~~~~~~~~

- GTF: This output format contains besides the genomic loci the annotated information for the genes giving rise to each peptide sequence including status and biotype. For each mapped peptide the sample, number of peptide-spectrum matches and associated quantitative value as tags.

- GCT: In this format the peptide sequences are combines with the Ensembl gene identifier. It contains the genomic loci for each peptide as well as the quantitative values for each peptide in different samples as a matrix.

Usage
-------------------

**Required arguments**:

- **-fasta**: Filepath for file containing protein sequences in FASTA format. (e.g. ``` -fasta gencode.v25.pc_translations.fa```)
- **-gtf**:   Gene annotation with coding sequences (CDS) in GTF format. (e.g. ``` -gtf gencode.v25.annotation.gtf```)
- **-in**:    Path to single input file or comma separated list of paths to input files containing peptides to be mapped with associated number of peptide to spectrum matches, sample name and quantitative value (see input file format). (e.g. ``` -in file.tsv```)

How to easily run the tool (**e.g. Human**)::

    $ java -jar -Xmx5G PepGenome-{version}.jar -gtf gencode-{version}.gtf
      -fasta gencode-{version}-translations.fa -in file.tsv


Optional arguments
-------------------

- **-format** : Set output format _GTF_, _GCT_, _BED_, _PTMBED_ or _ALL_. Comma separated combination possible. ```Default = ALL```
- **-merge**:   Set TRUE/FALSE to merge output of multiple input files (output will be named after last input file *_merged). ```Default = FALSE``
- **-source**:  Set TRUE/FALSE to merge output of multiple input files (output will be named after last input file *_merged). ```Default = FALSE```
- **-mm** : Number of mismatches allowed in mapping (0, 1 or 2). ```DEFAULT = 0```
- **-mmmode**: Set TRUE/FALSE to restrict number of mismatch in kmer to 1. ```DEFAULT = FALSE```.
- **-species**: Set species using common or scientific name or taxonomy ID (e.g. ```-species 9606```). Default is Human (Homo sapiens, 9606).
- **-chr**:  Export chr prefix Allowed 0, 1. (e.g. ```-chr 1```)  ```DEFAULT = 0```

Table of supported species
--------------------------

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
|Dog                |Canis lupus        |9615       |
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


