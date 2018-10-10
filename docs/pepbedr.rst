.. _pepbedr

PepBedR: R package to analyze peptidoforms features in Bed Files
==========================================

The PepBedR package allows users of :ref:`bed` file format to analyze and visualized their data and results. It facilitate two major things:

- Descriptive statistics of PepBed files: Number of unique peptides per chromosome, transcript, gene.
- Circle plots of peptide features by different Peptide properties (Modification, Uniqueness, Sample properties).

.. note:: The PepBedR provides set ```Rscript``` utilities to describe PepBed files.


Using R package
------------------

Parsing Bed file
~~~~~~~~~~~~~~~~~

.. code-block:: R
   :linenos:

   bed_path <- '/home/biolinux/temp/human/pride_cluster_peptides_9606_Human_pogo.bed'

   granges_peptide <- importBEDasGRange(inputFile = bed_path)

   n_features <- length(granges_peptide)
   message(c('Imported ', n_features, ' peptides...'))

Computing some basic stats from the data
~~~~~~~~~~~~~~~~~~~~~

```{r , message=FALSE, warning=FALSE}
# getting number of features(peptides) by chromosome
counts <- countsByChromosome(gr = granges_peptide, colName = 'Peptides')
counts_mod <- countsByChromosome(gr = granges_mod_peptide, colName = 'Peptides_mod')

# merging dfs
merged_counts <- merge.data.frame(counts, counts_mod, by = 'Chromosome')
# ordering by chromosome
merged_counts <- orderByChromosome(df = merged_counts, colName = 'Chromosome')

print(merged_counts)
```

\newpage

### Getting stats for unique peptides

```{r , message=FALSE, warning=FALSE}
# removing duplicated entries from original granges_peptide
unique_pep <- getUniqueFeatures(granges_peptide, colFeatures = 'name')
unique_pep_mod <- getUniqueFeatures(granges_mod_peptide, colFeatures = 'name')

# getting unique number of features(peptides) by chromosome
counts_unique <- countsByChromosome(gr = unique_pep, colName = 'Peptides')
counts_mod_unique <- countsByChromosome(gr = unique_pep_mod, colName = 'Peptides_mod')

# merging dfs
merged_counts_unique <- merge.data.frame(counts_unique,
                                         counts_mod_unique,
                                         by = 'Chromosome')

# ordering by chromosome
merged_counts_unique <- orderByChromosome(df = merged_counts_unique,
                                          colName = 'Chromosome')

print(merged_counts_unique)
```

\newpage

### Computing % coverage

```{r , message=FALSE, warning=FALSE}
## compute coverage of query (peptide evidences) on subject (transcripts) by crhomosome
data("protein_coding_transcript_hg38") # load protein coding transcript as GRanges object

coverage <- computeCoverageByChromosome(query = granges_peptide,
                                        subject = transcripts_hg38,
                                        colName = 'Coverage')

coverage_mod <- computeCoverageByChromosome(query = granges_mod_peptide,
                                            subject = transcripts_hg38,
                                            colName = 'Coverage_mod')

# merging dfs
merged_coverage <- merge.data.frame(coverage, coverage_mod, by = 'Chromosome')

# ordering by chromosome
merged_coverage <- orderByChromosome(df = merged_coverage, colName = 'Chromosome')

print(merged_coverage)
```

\newpage

### Visualizing the data

* The distribution of peptides by chromosome. (blue_track: modified peptide;  red_track: non-modified)

```{r density, echo=FALSE, fig.width=12, fig.height=12}
library(circlize)
circos.initializeWithIdeogram(species = 'hg19')
bed <- bed_df
bed_mod <- bed_mod_df
circos.genomicDensity(bed, col = c("#FF000080"), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_mod, col = c("#0000FF80"), track.height = 0.1, baseline = 0)
circos.clear()
```

\newpage

* The distribution of peptides by chromosome. (yellow_track: oxidation)

```{r density_oxidation, echo=FALSE, fig.width=12, fig.height=12}
library(circlize)
circos.initializeWithIdeogram(species = 'hg19')
bed <- bed_df
bed_oxidation <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'oxidation') # mod oxidation
circos.genomicDensity(bed, col = c(rgb(0, 0, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_oxidation, col = c(rgb(204, 204, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.clear()
```


\newpage

* The distribution of peptides by chromosome. (orange_track: acetyl)

```{r density_acetyl, echo=FALSE, fig.width=12, fig.height=12}
library(circlize)
circos.initializeWithIdeogram(species = 'hg19')
bed <- bed_df
bed_acetyl <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'acetyl') # mod acetylation
circos.genomicDensity(bed, col = c(rgb(0, 0, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_acetyl, col = c(rgb(204, 102, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.clear()
```


\newpage

* The distribution of peptides by chromosome. (red_track: phospho)

```{r density_phospho, echo=FALSE, fig.width=12, fig.height=12}
library(circlize)
circos.initializeWithIdeogram(species = 'hg19')
bed <- bed_df
bed_phospho <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'phospho') # mod phosphorylation
circos.genomicDensity(bed, col = c(rgb(0, 0, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_phospho, col = c(rgb(255, 51, 51, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.clear()
```


\newpage

* The distribution of peptides by chromosome. (all peptides vs. oxidation, acetyl, phospho)

```{r density_all_mod, echo=FALSE, fig.width=12, fig.height=12}

library(circlize)

circos.initializeWithIdeogram(species = 'hg19')

bed <- bed_df # all peptides
bed_oxidation <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'oxidation') # mod oxidation
bed_acetyl    <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'acetyl') # mod acetylation
bed_phospho   <- getModifiedSeq(bed_mod_df, colName = 'name', modPattern = 'phospho') # mod phosphorylation

circos.genomicDensity(bed, col = c(rgb(0, 0, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_oxidation, col = c(rgb(204, 204, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_acetyl, col = c(rgb(204, 102, 0, maxColorValue = 255)), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_phospho, col = c(rgb(255, 51, 51, maxColorValue = 255)), track.height = 0.1, baseline = 0)

circos.clear()

```


\newpage

* barplot with number of peptides (modified and non-modified) by chromosome

```{r counts, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- reshape2::melt(merged_counts)

plot1 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = value, fill=variable)) +
         geom_col(position = 'dodge', alpha = 0.7) +
         labs(x = 'Chromosome', y = 'Number of Peptides', fill = '') +
         scale_fill_discrete("Peptides", labels=c("non-modified", "modified")) +
         theme_bw() +
         theme(axis.text = element_text(size=12),
               axis.title = element_text(size=14),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"))
plot1
```

\newpage

* barplot with number of unique peptides (modified and non-modified) by chromosome

```{r unique_counts, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- reshape2::melt(merged_counts_unique)

plot2 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = value, fill=variable)) +
         geom_col(position = 'dodge', alpha = 0.7) +
         labs(x = 'Chromosome', y = 'Number of Peptides', fill = '') +
         scale_fill_discrete("Unique peptides", labels=c("non-modified", "modified")) +
         theme_bw() +
         theme(axis.text = element_text(size=12),
               axis.title = element_text(size=14),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"))
plot2
```

\newpage

* barplot with coverage (all peptides) by chromosome

```{r coverage, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- coverage

plot3 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = Coverage)) +
         geom_col(fill='darkgreen', alpha=0.4) +
         labs(x = 'Chromosome', y = '% coverage', fill = '') +
         theme_bw() +
         theme(axis.text = element_text(size=12),
               axis.title = element_text(size=14),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"))
plot3
```
