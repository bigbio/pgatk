# Use Cases and Recipes

This page provides end-to-end workflows for common proteogenomics scenarios.
Each recipe shows the exact commands to run, from downloading data to producing
a search-ready protein database.

---

## 1. Human Variant Protein Database from ENSEMBL

Build a variant protein database using ENSEMBL population variants (common SNPs
and indels) for human proteogenomics searches.

### Step 1 -- Download ENSEMBL data

Download the GTF, CDS, and VCF files for *Homo sapiens* (taxonomy 9606):

```bash
pgatk ensembl-downloader \
    -t 9606 \
    -o ensembl_human \
    --skip_protein \
    --skip_ncrna \
    --skip_cdna \
    --skip_dna
```

This downloads the gene annotation GTF and the VCF file with known variants.

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

## 2. High-Frequency Variant Database (AF Filtering)

Build a protein database that only includes common population variants above a
given allele frequency threshold.

```bash
pgatk vcf-to-proteindb \
    --vcf homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --af_field MAF \
    --af_threshold 0.01 \
    --output_proteindb common_variant_proteins.fa
```

To restrict to specific consequence types (e.g. only missense variants):

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

---

## 3. ClinVar Clinical Variant Database

Build a protein database from ClinVar pathogenic and likely pathogenic variants
using the NCBI RefSeq gene models.

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

## 4. Cancer-Specific Database from COSMIC

Generate tissue-specific protein databases from COSMIC somatic mutations.

### Step 1 -- Download COSMIC data

```bash
pgatk cosmic-downloader \
    -u your_email@example.com \
    -p your_password \
    -o cosmic_data
```

!!! note
    A COSMIC account is required. Register at
    [https://cancer.sanger.ac.uk/cosmic/register](https://cancer.sanger.ac.uk/cosmic/register).

### Step 2a -- Single database from all mutations

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_all_proteins.fa
```

### Step 2b -- One database per cancer type

Split mutations by primary tissue site, creating one FASTA per cancer type:

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_proteins.fa \
    --split_by_filter_column
```

This produces files named `cosmic_proteins_<tissue>.fa` for each tissue type
(e.g. `cosmic_proteins_lung.fa`, `cosmic_proteins_breast.fa`).

### Step 2c -- Database for a specific cancer type

Limit to a single tissue type:

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicMutantExport.tsv.gz \
    --input_genes cosmic_data/All_COSMIC_Genes.fasta.gz \
    --output_db cosmic_lung_proteins.fa \
    --accepted_values "lung"
```

### Step 2d -- Cell-line specific databases

Use cell-line data split by sample name:

```bash
pgatk cosmic-to-proteindb \
    --input_mutation cosmic_data/CosmicCLP_MutantExport.tsv.gz \
    --input_genes cosmic_data/All_CellLines_Genes.fasta.gz \
    --output_db cosmic_cellline_proteins.fa \
    --filter_column "Sample name" \
    --split_by_filter_column
```

---

## 5. Cancer Study from cBioPortal

Generate a protein database from cancer mutation data hosted in cBioPortal.

### Step 1 -- List available studies

```bash
pgatk cbioportal-downloader --list_studies
```

### Step 2 -- Download a specific study

```bash
pgatk cbioportal-downloader \
    -d blca_mskcc_solit_2014 \
    -o cbioportal_data
```

### Step 3 -- Download ENSEMBL CDS (hg19 assembly)

cBioPortal mutations are aligned to hg19 (GRCh37):

```bash
pgatk ensembl-downloader \
    -t 9606 \
    --grch37 \
    -o ensembl_hg19 \
    --skip_vcf \
    --skip_gtf \
    --skip_protein \
    --skip_ncrna \
    --skip_cdna \
    --skip_dna
```

### Step 4 -- Translate mutations to protein sequences

```bash
pgatk cbioportal-to-proteindb \
    --input_mutation cbioportal_data/data_mutations_mskcc.txt \
    --input_cds ensembl_hg19/Homo_sapiens.GRCh37.cds.all.fa.gz \
    --output_db bladder_cancer_proteins.fa
```

To generate tissue-type-specific databases using the clinical sample file:

```bash
pgatk cbioportal-to-proteindb \
    --input_mutation cbioportal_data/data_mutations_mskcc.txt \
    --input_cds ensembl_hg19/Homo_sapiens.GRCh37.cds.all.fa.gz \
    --clinical_sample_file cbioportal_data/data_clinical_sample.txt \
    --output_db cbioportal_proteins.fa \
    --split_by_filter_column
```

---

## 6. Non-Coding RNA Protein Database

Translate non-coding RNA transcripts to search for novel peptides from
supposedly non-coding regions.

### lincRNA database

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb lincRNA_proteins.fa \
    --var_prefix lincRNA_ \
    --include_biotypes lincRNA \
    --num_orfs 3 \
    --skip_including_all_cds
```

### Pseudogene database

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

Translate non-canonical reading frames of protein-coding genes:

```bash
pgatk dnaseq-to-proteindb \
    --input_fasta transcripts.fa \
    --output_proteindb altorf_proteins.fa \
    --var_prefix altorf_ \
    --include_biotypes altORFs \
    --skip_including_all_cds
```

---

## 7. Six-Frame Genome Translation

Perform a full six-frame translation of a genome assembly. This is useful as
an unbiased search space but produces very large databases.

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
    combination with a smaller targeted database search.

---

## 8. gnomAD Population Variant Database

Build a protein database from gnomAD variants using GENCODE annotations.

```bash
pgatk vcf-to-proteindb \
    --vcf gnomad_genome.vcf.gz \
    --input_fasta gencode_transcripts.fa \
    --gene_annotations_gtf gencode.v38.annotation.gtf.gz \
    --annotation_field_name vep \
    --af_field controls_AF \
    --af_threshold 0.01 \
    --include_consequences missense_variant,inframe_insertion,inframe_deletion,frameshift_variant \
    --biotype_str transcript_type \
    --output_proteindb gnomad_proteins.fa
```

!!! tip "gnomAD-specific parameters"
    - `--annotation_field_name vep` -- gnomAD uses `vep` instead of `CSQ`
    - `--af_field controls_AF` -- Use control allele frequency
    - `--biotype_str transcript_type` -- GENCODE uses `transcript_type` instead
      of ENSEMBL's `transcript_biotype`

---

## 9. Sample-Specific Database from WGS/WES

Translate variants from a whole-genome or whole-exome sequencing VCF of an
individual sample. These VCFs typically lack VEP annotations.

```bash
pgatk vcf-to-proteindb \
    --vcf sample_variants.vcf \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --annotation_field_name '' \
    --ignore_filters \
    --output_proteindb sample_proteins.fa
```

- `--annotation_field_name ''` tells the tool to skip parsing VEP/CSQ
  annotations (not present in raw variant callers like GATK HaplotypeCaller)
- `--ignore_filters` includes all variants regardless of the FILTER column

To include only variants that passed quality filtering:

```bash
pgatk vcf-to-proteindb \
    --vcf sample_variants.vcf \
    --input_fasta transcripts.fa \
    --gene_annotations_gtf genes.gtf \
    --annotation_field_name '' \
    --accepted_filters PASS \
    --output_proteindb sample_proteins_pass.fa
```

---

## 10. Digest and Filter Variant Peptides

After generating a variant protein database you may want to digest the proteins
*in silico* and keep only peptides that differ from the canonical proteome.
This produces a compact FASTA of variant-specific peptides.

```bash
pgatk digest-mutant-protein \
    --input variant_proteins.fa \
    --fasta canonical_proteins.fa \
    --output unique_variant_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

To combine multiple variant sources (e.g. COSMIC + ClinVar):

```bash
pgatk digest-mutant-protein \
    --input cosmic_proteins.fa,clinvar_proteins.fa \
    --fasta canonical_proteins.fa \
    --output combined_variant_peptides.fa
```

---

## 11. Map Identified Peptides to Genomic Coordinates

After a proteomics search, map the identified peptides back to genomic
coordinates for visualization in genome browsers (e.g. UCSC, IGV).

```bash
pgatk map-peptide2genome \
    --input peptide_identifications.tsv \
    --gtf genes.gtf \
    --fasta proteins.fa \
    --idmap protein_to_transcript.tsv \
    --output peptides.gff3
```

The input TSV should contain at least peptide sequence and protein accession
columns (configurable with `--pep-col` and `--prot-col`).

---

## 12. Complete Proteogenomics Workflow

A full end-to-end workflow combining canonical proteins, population variants,
clinical variants, cancer mutations, and non-coding translations into a single
search-ready database.

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
```

### Step 4 -- Combine all databases

```bash
cat canonical.fa \
    ensembl_variants.fa \
    clinvar_variants.fa \
    cosmic_variants.fa \
    lincRNA.fa \
    > combined_target.fa
```

### Step 5 -- Add decoy sequences

```bash
pgatk generate-decoy \
    --input combined_target.fa \
    --output proteogenomics_target_decoy.fa \
    --method decoypyrat \
    --decoy_prefix DECOY_
```

### Step 6 (optional) -- Extract unique variant peptides

```bash
pgatk digest-mutant-protein \
    --input ensembl_variants.fa,clinvar_variants.fa,cosmic_variants.fa \
    --fasta canonical.fa \
    --output unique_variant_peptides.fa \
    --min-len 7 \
    --max-len 40 \
    --missed-cleavages 2
```

The file `proteogenomics_target_decoy.fa` is a comprehensive search database,
and `unique_variant_peptides.fa` provides a focused list of variant-specific
peptides for validation.

---

## 13. Non-Human Species

pgatk supports any species available in ENSEMBL. Here is an example for mouse
(*Mus musculus*, taxonomy 10090):

```bash
# Download
pgatk ensembl-downloader -t 10090 -o ensembl_mouse --skip_dna

# Generate transcript sequences
gffread -F -w ensembl_mouse/transcripts.fa \
    -g ensembl_mouse/genome.fa \
    ensembl_mouse/Mus_musculus.GRCm39.*.gtf.gz

# Canonical proteins
pgatk dnaseq-to-proteindb \
    --input_fasta ensembl_mouse/transcripts.fa \
    --output_proteindb mouse_canonical.fa

# Variant proteins
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

## 14. Three-Frame Translation

Perform a simple three-frame translation of any nucleotide FASTA file:

```bash
pgatk threeframe-translation \
    --input_fasta input_sequences.fa \
    --output translated_proteins.fa
```

---

## 15. Quality Check: Validate Protein Database

Before using a protein database, check it for internal stop codons and short
sequences:

```bash
pgatk ensembl-check \
    --input_fasta protein_database.fa \
    --output validated_database.fa \
    --num_aa 6
```

This filters out sequences shorter than 6 amino acids. Use `--add_stop_codons`
to include proteins that contain internal stop codons (useful for
proteogenomics where frameshifts may introduce premature stops).
