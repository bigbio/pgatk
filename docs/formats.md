# PGATK File Formats

The ProteoGenomics Analysis ToolKit is based on standard proteomics formats developed by [HUPO-PSI](https://github.com/HUPO-PSI) and Genomics Standard file formats. This section highlights the most important features of those file formats, how they are used in PGATK, and how you can contribute to their development.

!!! note
    This page only provides the fundamentals of each file format used in PGATK. For more details we provide links to the original documentation of each file format.

## BED Format

BED **(Browser Extensible Data)** format provides a flexible way to define the data lines that are displayed in an annotation track ([UCSC BED Definition](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)). BED lines have three required fields and nine additional optional fields. The number of fields per line must be consistent throughout any single set of data in an annotation track. The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used.

!!! note
    If your data set is BED-like, but it is very large (over 50MB) and you would like to keep it on your own server, you should use the [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html) data format.

The **pedBed** Fields and Properties supported by PGATK:

| Field (bold are required) | Description | Example |
|---|---|---|
| **chrom** | The name of the chromosome | chr3 |
| **chromStart** | The starting position of the feature in the chromosome or scaffold | 1000 |
| **chromEnd** | The ending position of the feature in the chromosome or scaffold | 5000 |
| **name** | Defines the label of the BED line. PGATK annotates peptide sequence | K(Phospho)SR |
| **score** | A score between 0 and 1000. | 1000 |
| **strand** | Defines the strand. Either "." (=no strand) or "+" or "-". | + |
| **thickStart** | The starting position at which the feature is drawn thickly. | |
| **thickEnd** | The ending position at which the feature is drawn thickly | 5000 |
| **itemRgb** | An RGB value that will determine the display color of BED line (see [Color](#color)) | (0,0,255) |
| blockCount | The number of blocks (exons) in the BED line. | |
| blockSizes | A comma-separated list of the block sizes. | |
| blockStarts | A comma-separated list of block starts. | |
| proteinAccession | Protein accession number | |
| transcriptAccession | Transcript Accession | |
| peptideSequence | Peptide Sequence with no PTMs added | |
| proteinUniqueness | Peptide uniqueness (See color code [Color](#color)) | |
| transcriptUniqueness | Peptide uniqueness (See color code [Color](#color)) | |
| genomeReferenceVersion | Genome reference version number | |
| psmScore | Best PSM score | |
| fdr | False-discovery rate | |
| modifications | Comma separated list of Post-translational modifications | |
| peptideRepeatCount | Peptide Counting | |
| datasetAccession | Dataset Identifier | |
| uri | Uniform Resource Identifier | |

!!! tip "Hint"
    If the field content is to be empty the space should be filled with a "."

!!! note
    BED input files (and input received from stdin) are tab-delimited. The following types of BED files are supported by PGATK:

    - **BED4**: A BED file where each feature is described by chrom, start, end, and name. (e.g. `chr1  11873  14409  VLADIMIR`)
    - **BED6**: A BED file where each feature is described by chrom, start, end, name, score, and strand. (e.g. `chr1 11873 14409 VLADIMIR 0 +`)
    - **BED12**: A BED file where each feature is described by all twelve columns listed above. (Default option in all tools)
    - **BED12+11**: A complete BED file including **required** fields and optionals.

### Color

**Uniqueness** Colors:

| Colour | Description |
|---|---|
| <img src="images/uniquetranscript.svg" width="25"> | Peptide is unique to single gene AND single transcript |
| <img src="images/uniquegene.svg" width="25"> | Peptide is unique to single gene BUT shared between multiple transcripts |
| <img src="images/notunique.svg" width="25"> | Peptide is shared between multiple genes |

**Modified** Peptides Colors:

Like BED but containing the location of the post-translational modification on the genome. Thick parts of the peptide blocks indicate the position of the post-translational modification on a single amino acid (short thick block) while longer blocks indicate the occurrence of the first and last post-translational modification and residues in between. In the PTMBED the colour code is changed to indicate the type of modification.

| Colour | Post-translational Modification |
|---|---|
| <img src="images/phospho.svg" width="25"> | Phosphorylation (phospho) |
| <img src="images/acetyl.svg" width="25"> | Acetylation (acetyl) |
| <img src="images/amidated.svg" width="25"> | Amidation (amidated) |
| <img src="images/oxidation.svg" width="25"> | Oxidation (oxidation) |
| <img src="images/methyl.svg" width="25"> | Methylation (methyl) |
| <img src="images/glygly.svg" width="25"> | Ubiquitinylation (glygly; gg) |
| <img src="images/sulfo.svg" width="25"> | Sulfation (sulfo) |
| <img src="images/palmitoyl.svg" width="25"> | Palmitoylation (palmitoyl) |
| <img src="images/formyl.svg" width="25"> | Formylation (formyl) |
| <img src="images/deamidated.svg" width="25"> | Deamidation (deamidated) |
| <img src="images/any.svg" width="25"> | Any other post-translational modification |

## Additional File Formats

### Peptide Atlas Peptide List

[PeptideAtlas](http://www.peptideatlas.org) releases every month/year a list of peptides that have been found/identified by MS/MS (see the list [here](https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/defaultBuildsPepsProts?_subtab=3)).

| Column | Field | Description |
|---|---|---|
| 1 | peptide_accession | Peptide Accession (PAp06389395) |
| 2 | peptide_sequence | Peptide Sequence |
| 3 | best_probability | Best Peptide Probability |
| 4 | n_observations | Spectral counting |
| .... | | More properties not used |

!!! tip "Hint"
    For our pipelines and tools the order of the column is important.

!!! note
    A full pipeline to map the PeptideAtlas peptide evidences to Genome Coordinates can be found here.
