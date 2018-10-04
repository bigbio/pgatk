
PGATK File Formats
=====================

The ProteoGenomics Analysis ToolKit is based on standard proteomics formats developed by `HUPO-PSI <https://github.com/HUPO-PSI>`_ and Genomics Standard file formats. This section is to highlight in 10 minutes the most important features of those file formats, How they are used in PGATK and you can contribute to their development.

.. note:: It is important to notice that this Help page only provides the fundamentals of each file format used in PGATK, for more details we provide links to the original documentation of the File format.

.. _bed:

BED
-------------------

BED ***(Browser Extensible Data)** format provides a flexible way to define the data lines that are displayed in an annotation track `UCSC Bed Definition <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_. BED lines have three required fields and nine additional optional fields. The number of fields per line must be consistent throughout any single set of data in an annotation track. The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used.


.. note:: If your data set is BED-like, but it is very large (over 50MB) and you would like to keep it on your own server, you should use the `bigBed <https://genome.ucsc.edu/goldenPath/help/bigBed.html>`_ data format.

The pedBed Fields and Properties generated by PGATK:

+---------------------------+----------------------------------------------------------------------------+-------------+
|Field (bold are required)  | Description                                                                | Example     |
+===========================+============================================================================+=============+
|**chrom**                  |The name of the chromosome                                                  |chr3         |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**chromStart**             |The starting position of the feature                                        |1000         |
|                           |in the chromosome or scaffold                                               |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**chromEnd**               |The ending position of the feature in the chromosome or scaffold            |5000         |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**name**                   |Defines the label of the BED line. GPATK annotate peptide sequence          |K(Phospho)SR |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**score**                  |A score between 0 and 1000.                                                 |1000         |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**strand**                 |Defines the strand. Either "." (=no strand) or "+" or "-".                  | +           |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**thickStart**             |The starting position at which the feature is drawn thickly.                |1000         |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**thickEnd**               |The ending position at which the feature is drawn thickly                   |5000         |
+---------------------------+----------------------------------------------------------------------------+-------------+
|**itemRgb**                |An RGB value that will determine the display color of BED line :ref:`color` |(0,0,255)    |
+---------------------------+----------------------------------------------------------------------------+-------------+
|blockCount                 |The number of blocks (exons) in the BED line.                               |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|blockSizes                 |A comma-separated list of the block sizes.                                  |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|blockStarts                |A comma-separated list of block starts.                                     |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|proteinAccession           |Protein accession number                                                    |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|transcriptAccession        |Transcript Accession                                                        |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|uniqueness                 |Peptide uniqueness (See color code :ref:`color`                             |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|genomeReferenceVersion     |Genome reference version number                                             |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|psmScore                   |Best PSM score                                                              |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|fdr                        |False-discovery rate                                                        |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|modifications              |Coma separated list of Post-translational modifications                     |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|datasetAccession           |Dataset Identifier                                                          |             |
+---------------------------+----------------------------------------------------------------------------+-------------+
|uri                        |Uniform Resource Identifier                                                 |             |
+---------------------------+----------------------------------------------------------------------------+-------------+


Color
~~~~~~~~~~~~

**Uniqueness** Colors:

+-----------------------------------------+---------------------------------------------------------------------------+
| Colour                                  | Description                                                               |
+=========================================+===========================================================================+
|.. image:: images/uniquetranscript.svg   | Peptide is unique to single gene AND single transcript                    |
|   :width: 25                            |                                                                           |
+-----------------------------------------+---------------------------------------------------------------------------+
|.. image:: images/uniquegene.svg         | Peptide is unique to single gene BUT shared between multiple transcripts  |
|   :width: 25                            |                                                                           |
+-----------------------------------------+---------------------------------------------------------------------------+
|.. image:: images/notunique.svg          | Peptide is shared between multiple genes                                  |
|   :width: 25                            |                                                                           |
+-----------------------------------------+---------------------------------------------------------------------------+

**Modified** Peptides Colors:

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