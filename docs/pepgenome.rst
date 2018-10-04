.. _pepgenome


PepGenome: Mapping Peptidoforms to Genome Coordinates
===========================

In proteogenomic analyses it is essential to know the loci (genome position) giving rise to **peptides** in order to improve genomic annotation and the functional characterization of protein products in their biological context. With next-generation sequencing of DNA and RNA for each sample studied by proteomic mass spectrometry integration and visualisation in a common coordinate system, i.e. the genome, is vital for systems biology. To facilitate this type of integration not only the genomic locations of modified peptides but specifically the genomic loci of associated with these modifications is required.


.. note:: The **PepGenome** tool quickly and efficiently identify genomic loci of peptides and post-translational modifications and couple these mappings with associated quantitative values. Using reference gene annotation (GTF files) and an associated transcript translations (Protein fasta files) our tool identifies the genomic loci of peptides given as input and generates output in different formats borrowed from genomics and transcriptomics which can be loaded in various genome browsers such as `UCSC Genome Browser https://genome.ucsc.edu/`_, `Ensembl Genome Browser http://www.ensembl.org/index.html`_.


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

.. note:: * Peptide sequence with PSI-MS nodification names in round brackets following the mpdified amino acid, e.g. PEPT(Phopsho)IDE for a phosphorylated threonine.


.. note:: In addition the tool support mzTab and mzIdentML File format input.

Output formats
------------------------

BED
~~~~~~~~~~

This format contains the genomic loci for peptides, the exon-structure, the peptide sequence, as well as a colour code for uniqueness of peptides within the genome.

+-----------------------------------------+------------------------------------------------------------------------------+
| Colour                                  | Description                                                                  |
+=========================================+==============================================================================+
|.. image:: images/uniquetranscript.svg   | Peptide is unique to single gene AND single transcript                       |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/uniquegene.svg         | Peptide is unique to single gene BUT shared between multiple transcripts     |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/notunique.svg          | Peptide is shared between multiple genes                                     |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+



PTMBED
~~~~~~~~~~~

Like BED but containing the location of the post-translational modification on the genome. Thick parts of the peptide blocks indicate the position of the post-translational modification on a single amino acid (short thick block) while longer blocks indicate the occurrence of the first and last post-translational modification and residues in between. In the PTMBED the colour code is changed to indicate the type of modification.

+-----------------------------------------+------------------------------------------------------------------------------+
| Colour                                  | Post-translational Modification                                              |
+=========================================+==============================================================================+
|.. image:: images/phospho.svg            | Phosphorylation (phospho)                                                    |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/acetyl.svg             | Acetylation (acetyl)                                                         |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/amidated.svg           | Amidation (amidated)                                                         |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/oxidation.svg          | Oxidation (oxidation)                                                        |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/oxidation.svg          | Oxidation (oxidation)                                                        |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/methyl.svg             | Methylation (methyl)                                                         |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/glygly.svg             | Ubiquitinylation (glygly; gg)                                                |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/sulfo.svg              | Sulfation (sulfo)                                                            |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/palmitoyl.svg          | Palmitoylation (palmitoyl)                                                   |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/formyl.svg             | Formylation (formyl)                                                         |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/deamidated.svg         | Deamidation (deamidated)                                                     |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+
|.. image:: images/any.svg                | Any other post-translational modification                                    |
|   :width: 25                            |                                                                              |
+-----------------------------------------+------------------------------------------------------------------------------+


GTF
~~~~~~~~~

This output format contains besides the genomic loci the annotated information for the genes giving rise to each peptide sequence including status and biotype. For each mapped peptide the sample, number of peptide-spectrum matches and associated quantitative value as tags.

GCT
~~~~~~~~~~~
In this format the peptide sequences are combines with the Ensembl gene identifier. It contains the genomic loci for each peptide as well as the quantitative values for each peptide in different samples as a matrix.

Usage
-------------------

**Required arguments**:

- **-fasta**: Filepath for file containing protein sequences in FASTA format. (e.g. ``` -fasta gencode.v25.pc_translations.fa```)
- **-gtf**:   Gene annotation with coding sequences (CDS) in GTF format. (e.g. ``` -gtf gencode.v25.annotation.gtf```)
- **-in**:    Path to single input file or comma separated list of paths to input files containing peptides to be mapped with associated number of peptide to spectrum matches, sample name and quantitative value (see input file format). (e.g. ``` -in file.tsv```)

How to easily run the tool (**e.g. Human**)::

    $ java -jar -Xmx5G PepGenome-{version}.jar -gtf gencode-{version}.gtf -fasta gencode-{version}-translations.fa -in file.tsv


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


