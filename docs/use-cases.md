# Use Cases and Recipes

This page provides end-to-end workflows for common proteogenomics scenarios.
Each recipe shows the exact commands to run, from downloading data to producing
a search-ready protein database.

---

## Cell-Type Specific Non-Canonical Peptide Discovery

!!! abstract "Featured workflow"
    This workflow reproduces the analysis from [Umer et al., *Bioinformatics* 2022](https://doi.org/10.1093/bioinformatics/btab838),
    where cell-type specific databases were generated for 65 cell lines,
    revealing non-canonical and cryptic peptides amounting to **>5% of total
    peptide identifications**.

The key insight is to build **cell-type specific** proteogenomics databases
that combine canonical proteins with non-canonical translations (pseudogenes,
lncRNAs, alternative ORFs) and variant sequences from multiple genomic sources.
Searching with a database tailored to the cell type of interest maximizes
discovery while keeping the search space focused.

### Step 1 -- Download ENSEMBL data

```bash
pgatk ensembl-downloader \
    -t 9606 \
    -o ensembl_human
```

### Step 2 -- Generate transcript sequences from GTF

```bash
gffread -F -w ensembl_human/transcripts.fa \
    -g ensembl_human/genome.fa \
    ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz
```

### Step 3 -- Canonical protein-coding sequences

Translate all protein-coding transcripts using their annotated CDS:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_human/transcripts.fa \
    --output_proteindb canonical.fa
```

### Step 4 -- Non-canonical translations

#### Pseudogenes (three-frame)

Pseudogenes were traditionally considered non-functional, but some produce
detectable peptides. Translate them in three reading frames:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_human/transcripts.fa \
    --output_proteindb pseudogene.fa \
    --var_prefix pseudo_ \
    --include_biotypes processed_pseudogene,unprocessed_pseudogene,transcribed_processed_pseudogene,transcribed_unprocessed_pseudogene,translated_processed_pseudogene \
    --num_orfs 3 \
    --skip_including_all_cds
```

#### lncRNAs (three-frame)

Long non-coding RNAs can harbor small open reading frames encoding
micropeptides:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_human/transcripts.fa \
    --output_proteindb lncrna.fa \
    --var_prefix lncrna_ \
    --include_biotypes lincRNA,antisense,sense_intronic,sense_overlapping \
    --num_orfs 3 \
    --skip_including_all_cds
```

#### Alternative ORFs (exonic out-of-frame translation)

Non-canonical reading frames of protein-coding mRNAs can produce cryptic
peptides:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_human/transcripts.fa \
    --output_proteindb altorf.fa \
    --var_prefix altorf_ \
    --include_biotypes altORFs \
    --skip_including_all_cds
```

### Step 5 -- Population variant proteins

Include common human variants from ENSEMBL:

```bash
pgatk vcf-to-proteindb \
    --vcf ensembl_human/homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta ensembl_human/transcripts.fa \
    --gene_annotations_gtf ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz \
    --output_proteindb ensembl_variants.fa
```

### Step 6 -- COSMIC somatic mutations (optional, for cancer cell lines)

For cancer cell-line studies, add cell-line-specific somatic mutations:

```bash
pgatk cosmic-downloader \
    -u your_email@example.com \
    -p your_password \
    -o cosmic_data

pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicCLP_MutantExport.tsv.gz \
    --input_genes cosmic_data/All_CellLines_Genes.fasta.gz \
    --output_db cosmic_cellline.fa \
    --filter_column "Sample name" \
    --split_by_filter_column
```

### Step 7 -- Combine and generate target-decoy database

```bash
# Combine all components
# For a specific cell line, use the matching COSMIC file if available
cat canonical.fa \
    pseudogene.fa \
    lncrna.fa \
    altorf.fa \
    ensembl_variants.fa \
    > cell_type_target.fa

# Generate decoy sequences
pgatk generate-decoy \
    --input cell_type_target.fa \
    --output cell_type_target_decoy.fa \
    --method decoypyrat \
    --decoy_prefix DECOY_
```

### Step 8 -- Extract non-canonical unique peptides

After database searching with your search engine, you can also pre-compute the
set of non-canonical peptides unique to these novel sources:

```bash
pgatk digest-mutant-protein \
    --input pseudogene.fa,lncrna.fa,altorf.fa,ensembl_variants.fa \
    --fasta canonical.fa \
    --output non_canonical_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

!!! tip "Expected results"
    In the original study, this workflow applied to 65 cell lines across six
    PRIDE datasets identified non-canonical peptides amounting to >5% of
    total identifications. The most common non-canonical sources were
    pseudogenes and lncRNAs, followed by alternative ORFs and variant
    peptides.

---

## 1. Human Variant Protein Database from ENSEMBL

Build a variant protein database using ENSEMBL population variants (common SNPs
and indels) for human proteogenomics searches. This is the most common starting
point for proteogenomics experiments -- augmenting the canonical proteome with
known variant peptides that would otherwise be missed.

### Step 1 -- Download ENSEMBL data

Download the GTF, CDS, and VCF files for *Homo sapiens* (taxonomy 9606):

```bash
pgatk ensembl-downloader \
    -t 9606 \
    -o ensembl_human \
    --skip_protein \
    --skip_ncrna \
    --skip_cdna
```

This downloads the gene annotation GTF, VCF file with known variants, and the
genome FASTA (needed by gffread to extract transcript sequences).

### Step 2 -- Generate transcript sequences

Use [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) to
extract transcript sequences from the GTF and genome FASTA:

```bash
gffread -F -w ensembl_human/transcripts.fa \
    -g ensembl_human/genome.fa \
    ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz
```

### Step 3 -- Generate the variant protein database

Translate all ENSEMBL variants that affect protein-coding transcripts:

```bash
pgatk vcf-to-proteindb \
    --vcf ensembl_human/homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta ensembl_human/transcripts.fa \
    --gene_annotations_gtf ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz \
    --output_proteindb ensembl_human/variant_proteins.fa
```

### Step 4 -- Generate the canonical protein database

Translate canonical protein-coding transcripts:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_human/transcripts.fa \
    --output_proteindb ensembl_human/canonical_proteins.fa
```

### Step 5 -- Combine and add decoy sequences

Merge canonical and variant databases, then generate decoys:

```bash
cat ensembl_human/canonical_proteins.fa \
    ensembl_human/variant_proteins.fa \
    > ensembl_human/target.fa

pgatk generate-decoy \
    --input ensembl_human/target.fa \
    --output ensembl_human/target_decoy.fa \
    --method decoypyrat \
    --decoy_prefix DECOY_
```

The file `target_decoy.fa` is ready for database searching.

---

## 2. Population-Specific Variant Database

Population-level genetic variants cause amino acid changes that are invisible
to standard reference database searches. Studies have shown that incorporating
common missense variants can increase peptide identifications by 3--5% and is
especially important when analyzing samples from underrepresented ancestries.

### Common variants above a frequency threshold

Include only variants present in at least 1% of the population:

```bash
pgatk vcf-to-proteindb \
    --vcf homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --af_field MAF \
    --af_threshold 0.01 \
    --output_proteindb common_variant_proteins.fa
```

### Missense-only variant database

Restrict to single amino acid variants (SAAVs), the most common variant type
detectable by mass spectrometry:

```bash
pgatk vcf-to-proteindb \
    --vcf homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --af_field MAF \
    --af_threshold 0.05 \
    --include_consequences missense_variant \
    --output_proteindb missense_common_proteins.fa
```

### gnomAD ancestry-stratified database

gnomAD provides allele frequencies stratified by ancestry (African, East Asian,
South Asian, European, Latino, etc.). Build a database using variants common in
a specific population:

```bash
pgatk vcf-to-proteindb \
    --vcf gnomad.exomes.v4.1.sites.vcf.bgz \
    --input_fasta gencode_transcripts.fa \
    --gene_annotations_gtf gencode.v44.annotation.gtf.gz \
    --annotation_field_name vep \
    --af_field AF_afr \
    --af_threshold 0.01 \
    --include_consequences missense_variant,inframe_insertion,inframe_deletion \
    --biotype_str transcript_type \
    --output_proteindb gnomad_afr_proteins.fa
```

!!! tip "gnomAD-specific parameters"
    - `--annotation_field_name vep` -- gnomAD uses `vep` instead of `CSQ`
    - `--af_field` -- Use population-specific AF fields: `AF_afr` (African),
      `AF_eas` (East Asian), `AF_sas` (South Asian), `AF_nfe` (Non-Finnish European),
      `AF_amr` (Latino), or `controls_AF` (all controls)
    - `--biotype_str transcript_type` -- GENCODE uses `transcript_type` instead
      of ENSEMBL's `transcript_biotype`

---

## 3. ClinVar Clinical Variant Database

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) catalogs the relationship
between human variants and clinical phenotypes. Building a ClinVar-derived
protein database enables detection of clinically annotated variant peptides by
mass spectrometry -- useful for clinical proteomics, pharmacoproteomics, and
validating pathogenic variants at the protein level.

### Step 1 -- Download NCBI / ClinVar files

```bash
pgatk ncbi-downloader -o ncbi_data
```

This downloads four files to `ncbi_data/`:

- `GRCh38_latest_genomic.gtf.gz` -- RefSeq gene annotations
- `GRCh38_latest_rna.fna.gz` -- RefSeq transcript nucleotide sequences
- `GRCh38_latest_assembly_report.txt` -- Chromosome name mapping
- `clinvar.vcf.gz` -- ClinVar variant calls

### Step 2 -- Generate the ClinVar protein database

```bash
pgatk clinvar-to-proteindb \
    --vcf ncbi_data/clinvar.vcf.gz \
    --gtf ncbi_data/GRCh38_latest_genomic.gtf.gz \
    --fasta ncbi_data/GRCh38_latest_rna.fna.gz \
    --assembly-report ncbi_data/GRCh38_latest_assembly_report.txt \
    --output clinvar_proteins.fa
```

### Step 3 -- Add decoy sequences

```bash
pgatk generate-decoy \
    --input clinvar_proteins.fa \
    --output clinvar_target_decoy.fa \
    --method decoypyrat
```

!!! tip
    The ClinVar pipeline does **not** require VEP annotations. It uses
    BedTools interval overlap to find transcripts affected by each variant,
    then applies the variant to the transcript sequence and translates it.

---

## 4. Tumor-Specific Databases for Cancer Proteogenomics

Cancer proteogenomics studies (such as those from CPTAC) build
tumor-specific protein databases to detect somatic mutant peptides, understand
therapy resistance, and prioritize neoantigen candidates. The databases combine
somatic mutations from cancer-specific sources (COSMIC, cBioPortal) and/or
patient-matched whole-exome sequencing.

### 4a. COSMIC somatic mutations by cancer type

Generate one protein database per primary tissue site. This is the standard
approach for large-scale cancer proteogenomics when patient-level WES is not
available:

```bash
# Download COSMIC data (requires account)
pgatk cosmic-downloader \
    -u your_email@example.com \
    -p your_password \
    -o cosmic_data

# Generate per-tissue databases
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_proteins.fa \
    --split_by_filter_column
```

This produces files like `cosmic_proteins_lung.fa`, `cosmic_proteins_breast.fa`,
etc. Use the tissue-matched database for your cancer type of interest.

### 4b. Single cancer type

Build a focused database for one cancer type (e.g. lung):

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_lung_proteins.fa \
    --accepted_values "lung"
```

### 4c. Cell-line proteogenomics

When analyzing cell-line proteomes, use cell-line-specific mutations. COSMIC
provides a dedicated cell-line export with mutations annotated per sample:

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicCLP_MutantExport.tsv.gz \
    --input_genes cosmic_data/All_CellLines_Genes.fasta.gz \
    --output_db cosmic_cellline_proteins.fa \
    --filter_column "Sample name" \
    --split_by_filter_column
```

### 4d. cBioPortal study-specific database

cBioPortal hosts mutation data from thousands of cancer genomics studies.
Generate a protein database from a specific study:

```bash
# List available studies
pgatk cbioportal-downloader --list_studies

# Download a study
pgatk cbioportal-downloader \
    -d brca_tcga_pan_can_atlas_2018 \
    -o cbioportal_data

# Download ENSEMBL CDS (hg19 -- cBioPortal uses GRCh37)
pgatk ensembl-downloader \
    -t 9606 --grch37 -o ensembl_hg19 \
    --skip_vcf --skip_gtf --skip_protein \
    --skip_ncrna --skip_cdna --skip_dna

# Translate mutations
pgatk cbioportal-to-proteindb \
    --input_mutation cbioportal_data/data_mutations_mskcc.txt \
    --input_cds ensembl_hg19/Homo_sapiens.GRCh37.cds.all.fa.gz \
    --output_db brca_tcga_proteins.fa
```

### 4e. Combined cancer database for immunopeptidomics

For HLA immunopeptidomics / neoantigen discovery, a broad mutation database
maximizes the chance of detecting mutant HLA-presented peptides. Studies have
shown that COSMIC-derived databases can identify 5x more mutant immunopeptides
than patient WES alone. Combine COSMIC with ClinVar:

```bash
# Generate COSMIC mutations
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_proteins.fa

# Generate ClinVar mutations
pgatk clinvar-to-proteindb \
    --vcf ncbi_data/clinvar.vcf.gz \
    --gtf ncbi_data/GRCh38_latest_genomic.gtf.gz \
    --fasta ncbi_data/GRCh38_latest_rna.fna.gz \
    --assembly-report ncbi_data/GRCh38_latest_assembly_report.txt \
    --output clinvar_proteins.fa

# Extract variant-unique peptides (undigested -- HLA peptides are not tryptic)
# For immunopeptidomics, use the full-length variant proteins directly:
cat cosmic_proteins.fa clinvar_proteins.fa > neoantigen_candidates.fa

pgatk generate-decoy \
    --input neoantigen_candidates.fa \
    --output neoantigen_target_decoy.fa \
    --method decoypyrat
```

---

## 5. Patient-Specific Database from WGS/WES

When matched whole-genome or whole-exome sequencing data is available for
a sample, build a personalized protein database from the patient's own
variants. This is the gold standard for precision proteogenomics as used in
CPTAC studies.

### From a raw VCF (no VEP annotations)

VCF files from variant callers (e.g. GATK HaplotypeCaller, Strelka, MuTect2)
typically lack VEP consequence annotations:

```bash
pgatk vcf-to-proteindb \
    --vcf patient_somatic.vcf \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --annotation_field_name '' \
    --ignore_filters \
    --output_proteindb patient_proteins.fa
```

- `--annotation_field_name ''` skips VEP/CSQ parsing
- `--ignore_filters` includes all called variants

### From a VEP-annotated VCF

If the VCF has been annotated with Ensembl VEP, you can filter by consequence
and transcript biotype:

```bash
pgatk vcf-to-proteindb \
    --vcf patient_somatic.vep.vcf \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --include_consequences missense_variant,frameshift_variant,stop_gained,inframe_insertion,inframe_deletion \
    --output_proteindb patient_proteins.fa
```

### Only PASS-filtered variants

Include only high-confidence variant calls:

```bash
pgatk vcf-to-proteindb \
    --vcf patient_somatic.vcf \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --annotation_field_name '' \
    --accepted_filters PASS \
    --output_proteindb patient_pass_proteins.fa
```

---

## 6. Novel ORF and Micropeptide Discovery

Proteogenomics is a key approach for discovering novel coding regions: small
open reading frames (smORFs), micropeptides from lncRNAs, pseudogene-encoded
proteins, and alternative reading frames. Studies have found that non-canonical
peptides can account for over 5% of total identifications.

### lincRNA-derived proteins

Long intergenic non-coding RNAs (lincRNAs) can encode small proteins. Translate
them in three reading frames:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb lincRNA_proteins.fa \
    --var_prefix lincRNA_ \
    --include_biotypes lincRNA \
    --num_orfs 3 \
    --skip_including_all_cds
```

### Pseudogene-derived proteins

Some pseudogenes are transcribed and may produce functional peptides:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb pseudogene_proteins.fa \
    --var_prefix pseudogene_ \
    --include_biotypes processed_pseudogene,transcribed_processed_pseudogene,translated_processed_pseudogene \
    --num_orfs 3 \
    --skip_including_all_cds
```

### Alternative ORFs from protein-coding transcripts

Translate non-canonical reading frames of known protein-coding genes to
discover upstream ORFs (uORFs), overlapping ORFs, and downstream ORFs:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb altorf_proteins.fa \
    --var_prefix altorf_ \
    --include_biotypes altORFs \
    --skip_including_all_cds
```

### Antisense transcripts

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb antisense_proteins.fa \
    --var_prefix antisense_ \
    --include_biotypes antisense,antisense_RNA \
    --num_orfs 3 \
    --skip_including_all_cds
```

### Combined non-canonical database

Combine all non-canonical sources with the canonical proteome:

```bash
cat canonical_proteins.fa \
    lincRNA_proteins.fa \
    pseudogene_proteins.fa \
    altorf_proteins.fa \
    antisense_proteins.fa \
    > novel_orf_target.fa

pgatk generate-decoy \
    --input novel_orf_target.fa \
    --output novel_orf_target_decoy.fa \
    --method decoypyrat
```

!!! tip
    After searching, use `digest-mutant-protein` to extract only the novel
    peptides that do not appear in the canonical proteome:

    ```bash
    pgatk digest-mutant-protein \
        --input lincRNA_proteins.fa,pseudogene_proteins.fa,altorf_proteins.fa \
        --fasta canonical_proteins.fa \
        --output novel_unique_peptides.fa
    ```

---

## 7. Genome Annotation Refinement

Proteogenomics is used to validate and refine genome annotations by mapping
proteomics-identified peptides back to genomic coordinates. This has been
applied to improve annotations for organisms from wheat (validating 33,612
gene models) to human (discovering novel splice junctions and exons).

### Six-frame genome translation

Translate the entire genome in all six reading frames to create an unbiased
search space:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta genome.fa \
    --output_proteindb genome_6frame.fa \
    --biotype_str '' \
    --num_orfs 3 \
    --num_orfs_complement 3
```

!!! warning
    Six-frame translation of a full mammalian genome produces a very large
    database. Consider translating individual chromosomes or using it in
    combination with a smaller targeted database for a two-pass search strategy.

### Three-frame transcriptome translation

A more focused approach: translate all annotated transcripts in three reading
frames to capture alternative ORFs without the full genome search space:

```bash
pgatk threeframe-translation \
    --input_fasta transcripts.fa \
    --output translated_transcripts.fa
```

### Map identified peptides back to the genome

After database searching, map novel peptide identifications to genomic
coordinates for annotation and genome browser visualization:

```bash
pgatk map-peptide2genome \
    --input novel_peptide_identifications.tsv \
    --gtf genes.gtf \
    --fasta proteins.fa \
    --idmap protein_to_transcript.tsv \
    --output novel_peptides.gff3
```

The GFF3 output can be loaded into UCSC Genome Browser, IGV, or Ensembl for
visual inspection alongside gene models, RNA-seq coverage, and conservation
tracks.

---

## 8. Metaproteomics: Microbiome and Environmental Samples

Metaproteomics characterizes the functional activity of microbial communities.
The central challenge is building the protein search database, since communities
can contain hundreds to thousands of species. When matched metagenomic
sequencing data is available, six-frame translation of assembled contigs
provides a sample-matched database.

### Six-frame translation of metagenome assembly

After assembling metagenomic reads into contigs (e.g. with megahit or
metaSPAdes), translate them:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta metagenome_contigs.fa \
    --output_proteindb metagenome_6frame.fa \
    --biotype_str '' \
    --num_orfs 3 \
    --num_orfs_complement 3 \
    --protein_prefix meta_
```

### Three-frame translation of predicted ORFs

If gene prediction has been run on the metagenome assembly (e.g. Prodigal,
MetaGeneAnnotator), translate the predicted CDS regions:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta predicted_genes.fna \
    --output_proteindb metagenome_proteins.fa \
    --biotype_str '' \
    --protein_prefix meta_
```

### Add decoy for metaproteomics search

```bash
pgatk generate-decoy \
    --input metagenome_proteins.fa \
    --output metagenome_target_decoy.fa \
    --method decoypyrat
```

!!! note
    For very large metaproteomics databases (millions of sequences),
    consider using `--memory_save` with `generate-decoy` to reduce memory
    usage, or splitting the database by contig/bin before decoy generation.

---

## 9. Long-Read Transcriptomics + Proteogenomics

Long-read RNA sequencing (PacBio Iso-Seq, Oxford Nanopore) captures full-length
transcript isoforms without assembly, enabling unambiguous protein isoform
prediction. Studies have shown that 25% of long-read-detected protein isoforms
may be absent from reference annotations.

### Translate long-read transcript sequences

After processing long reads through an isoform classification pipeline (e.g.
SQANTI3), translate the high-quality transcript sequences:

```bash
# Translate all long-read transcripts, including novel isoforms
pgatk dnaseq-to-proteindb \
    --input_fasta long_read_transcripts.fa \
    --output_proteindb long_read_proteins.fa \
    --biotype_str '' \
    --protein_prefix lr_
```

### Combine with canonical proteome

```bash
cat canonical_proteins.fa long_read_proteins.fa > lr_target.fa

pgatk generate-decoy \
    --input lr_target.fa \
    --output lr_target_decoy.fa \
    --method decoypyrat
```

### Extract novel isoform peptides

After searching, identify peptides unique to long-read isoforms that are not
present in the reference:

```bash
pgatk digest-mutant-protein \
    --input long_read_proteins.fa \
    --fasta canonical_proteins.fa \
    --output novel_isoform_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

---

## 10. Plant and Non-Model Organism Proteogenomics

For non-model organisms with draft or incomplete genome annotations,
proteogenomics improves gene model quality. This approach has been particularly
impactful for crop species -- for example, a wheat proteogenomics study
validated over 33,000 gene models and recommended promoting 3,700
low-confidence genes.

### Step 1 -- Download ENSEMBL data for the species

pgatk supports any species available in ENSEMBL. For example, rice
(*Oryza sativa*, taxonomy 39947):

```bash
pgatk ensembl-downloader \
    -t 39947 \
    -o ensembl_rice
```

### Step 2 -- Generate transcript sequences and canonical proteome

```bash
gffread -F -w ensembl_rice/transcripts.fa \
    -g ensembl_rice/genome.fa \
    ensembl_rice/Oryza_sativa.IRGSP-1.0.*.gtf.gz

pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_rice/transcripts.fa \
    --output_proteindb rice_canonical.fa
```

### Step 3 -- Three-frame translation for novel gene discovery

```bash
pgatk threeframe-translation \
    --input_fasta ensembl_rice/transcripts.fa \
    --output rice_3frame.fa
```

### Step 4 -- Six-frame genome translation (unbiased)

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_rice/genome.fa \
    --output_proteindb rice_genome_6frame.fa \
    --biotype_str '' \
    --num_orfs 3 \
    --num_orfs_complement 3
```

### Step 5 -- Combine databases

```bash
cat rice_canonical.fa rice_3frame.fa > rice_target.fa

pgatk generate-decoy \
    --input rice_target.fa \
    --output rice_target_decoy.fa \
    --method decoypyrat
```

### Step 6 -- Map identified novel peptides to the genome

```bash
pgatk map-peptide2genome \
    --input rice_novel_peptides.tsv \
    --gtf ensembl_rice/Oryza_sativa.IRGSP-1.0.*.gtf.gz \
    --fasta rice_target.fa \
    --idmap protein_to_transcript.tsv \
    --output rice_novel_peptides.gff3
```

### Mouse example

For model organisms like mouse (*Mus musculus*, taxonomy 10090):

```bash
# Download
pgatk ensembl-downloader -t 10090 -o ensembl_mouse

# Generate transcript sequences
gffread -F -w ensembl_mouse/transcripts.fa \
    -g ensembl_mouse/genome.fa \
    ensembl_mouse/Mus_musculus.GRCm39.*.gtf.gz

# Canonical + variant proteins
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_mouse/transcripts.fa \
    --output_proteindb mouse_canonical.fa

pgatk vcf-to-proteindb \
    --vcf ensembl_mouse/mus_musculus_incl_consequences.vcf.gz \
    --input_fasta ensembl_mouse/transcripts.fa \
    --gene_annotations_gtf ensembl_mouse/Mus_musculus.GRCm39.*.gtf.gz \
    --output_proteindb mouse_variants.fa

# Combine and generate decoy
cat mouse_canonical.fa mouse_variants.fa > mouse_target.fa
pgatk generate-decoy \
    --input mouse_target.fa \
    --output mouse_target_decoy.fa \
    --method decoypyrat
```

---

## 11. Digest and Filter Variant Peptides

After generating variant protein databases from any source, *in silico*
digestion and filtering against the canonical proteome extracts only the
peptides that are unique to the variant sequences. This compact peptide list
is useful for targeted validation, spectral library building, and
quantification of variant peptides.

### Single source

```bash
pgatk digest-mutant-protein \
    --input variant_proteins.fa \
    --fasta canonical_proteins.fa \
    --output unique_variant_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

### Combine multiple variant sources

```bash
pgatk digest-mutant-protein \
    --input cosmic_proteins.fa,clinvar_proteins.fa,ensembl_variants.fa \
    --fasta canonical_proteins.fa \
    --output combined_variant_peptides.fa
```

---

## 12. Complete Multi-Source Proteogenomics Workflow

A full end-to-end workflow combining canonical proteins, population variants,
clinical variants, cancer mutations, and non-coding translations into a single
search-ready database. This is the approach used in large-scale studies that aim
to maximize the discovery of non-canonical peptides.

### Step 1 -- Download all data sources

```bash
# ENSEMBL (canonical + population variants)
pgatk ensembl-downloader -t 9606 -o ensembl_data

# NCBI / ClinVar
pgatk ncbi-downloader -o ncbi_data

# COSMIC (requires account)
pgatk cosmic-downloader -u user@example.com -p password -o cosmic_data
```

### Step 2 -- Prepare transcript sequences

```bash
gffread -F -w ensembl_data/transcripts.fa \
    -g ensembl_data/genome.fa \
    ensembl_data/Homo_sapiens.GRCh38.*.gtf.gz
```

### Step 3 -- Generate all protein databases

```bash
# Canonical proteins
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_data/transcripts.fa \
    --output_proteindb canonical.fa

# ENSEMBL population variants
pgatk vcf-to-proteindb \
    --vcf ensembl_data/homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta ensembl_data/transcripts.fa \
    --gene_annotations_gtf ensembl_data/Homo_sapiens.GRCh38.*.gtf.gz \
    --output_proteindb ensembl_variants.fa

# ClinVar clinical variants
pgatk clinvar-to-proteindb \
    --vcf ncbi_data/clinvar.vcf.gz \
    --gtf ncbi_data/GRCh38_latest_genomic.gtf.gz \
    --fasta ncbi_data/GRCh38_latest_rna.fna.gz \
    --assembly-report ncbi_data/GRCh38_latest_assembly_report.txt \
    --output clinvar_variants.fa

# COSMIC somatic mutations
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_variants.fa

# Non-coding RNA (lincRNA)
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_data/transcripts.fa \
    --output_proteindb lincRNA.fa \
    --var_prefix lincRNA_ \
    --include_biotypes lincRNA \
    --num_orfs 3 \
    --skip_including_all_cds

# Pseudogenes
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_data/transcripts.fa \
    --output_proteindb pseudogene.fa \
    --var_prefix pseudogene_ \
    --include_biotypes processed_pseudogene,transcribed_processed_pseudogene \
    --num_orfs 3 \
    --skip_including_all_cds
```

### Step 4 -- Combine all databases

```bash
cat canonical.fa \
    ensembl_variants.fa \
    clinvar_variants.fa \
    cosmic_variants.fa \
    lincRNA.fa \
    pseudogene.fa \
    > combined_target.fa
```

### Step 5 -- Quality check and decoy generation

```bash
# Filter short sequences and validate
pgatk ensembl-check \
    --input_fasta combined_target.fa \
    --output validated_target.fa \
    --num_aa 6

# Add decoy sequences
pgatk generate-decoy \
    --input validated_target.fa \
    --output proteogenomics_target_decoy.fa \
    --method decoypyrat \
    --decoy_prefix DECOY_
```

### Step 6 (optional) -- Extract unique variant peptides

```bash
pgatk digest-mutant-protein \
    --input ensembl_variants.fa,clinvar_variants.fa,cosmic_variants.fa,lincRNA.fa,pseudogene.fa \
    --fasta canonical.fa \
    --output unique_variant_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

The file `proteogenomics_target_decoy.fa` is a comprehensive search database,
and `unique_variant_peptides.fa` provides a focused list of non-canonical
peptides for validation.

!!! tip "Database size considerations"
    Larger proteogenomic databases improve coverage but increase the FDR
    penalty. Consider a two-pass strategy: first search against a focused
    database (canonical + variants), then a broader database (adding
    non-coding ORFs and six-frame translations) for discovery.

---

## 13. Quality Check: Validate Protein Database

Before using a protein database for searching, validate it for internal stop
codons and short sequences:

```bash
pgatk ensembl-check \
    --input_fasta protein_database.fa \
    --output validated_database.fa \
    --num_aa 6
```

This filters out sequences shorter than 6 amino acids. Use `--add_stop_codons`
to include proteins that contain internal stop codons (useful for
proteogenomics where frameshifts may introduce premature stops).

---

## References

If you use pgatk in your research, please cite:

> Husen M Umer, Enrique Audain, Yafeng Zhu, Julianus Pfeuffer, Timo
> Sachsenberg, Janne Lehtiö, Rui M Branca, Yasset Perez-Riverol.
> **Generation of ENSEMBL-based proteogenomics databases boosts the
> identification of non-canonical peptides.**
> *Bioinformatics*, Volume 38, Issue 5, 1 March 2022, Pages 1470--1472.
> [https://doi.org/10.1093/bioinformatics/btab838](https://doi.org/10.1093/bioinformatics/btab838)
