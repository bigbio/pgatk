# pgatk -- ProteoGenomics Analysis Toolkit

![Python application](https://github.com/bigbio/pgatk/workflows/Python%20application/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgatk/README.html)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6d030fd7d69413987f7265a01193324)](https://www.codacy.com/gh/bigbio/pgatk/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bigbio/pgatk&amp;utm_campaign=Badge_Grade)
[![PyPI version](https://badge.fury.io/py/pgatk.svg)](https://badge.fury.io/py/pgatk)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pgatk)

**pgatk** is a Python toolkit for building proteogenomics protein sequence databases. It downloads, translates, and combines variant and non-canonical sequences from multiple genomic sources into search-ready FASTA databases compatible with all major proteomics search engines.

## Key Features

- **Multi-source variant integration** -- Translate variants from ENSEMBL, VCF files, COSMIC, cBioPortal, ClinVar, and gnomAD into protein sequences
- **Non-canonical ORF discovery** -- Three-frame and six-frame translation of lncRNAs, pseudogenes, antisense transcripts, and alternative reading frames
- **Any species** -- Supports all organisms available in ENSEMBL (human, mouse, rice, wheat, etc.)
- **Search engine compatible** -- Output FASTA files work with MaxQuant, SearchGUI, MSFragger, Comet, DIA-NN, and Proteome Discoverer
- **Decoy generation** -- Multiple target-decoy strategies (DecoyPYrat, protein-reverse, protein-shuffle)
- **Peptide-to-genome mapping** -- Map identified peptides back to genomic coordinates (GFF3) for genome browser visualization
- **ClinVar without VEP** -- ClinVar pipeline uses BedTools interval overlap, no VEP annotation required

## Installation

### pip (recommended)

```bash
pip install pgatk
```

### Bioconda

```bash
conda install -c bioconda pgatk
```

### From source

```bash
git clone https://github.com/bigbio/pgatk.git
cd pgatk
pip install .
```

## Quick Start

Build a human variant protein database in four commands:

```bash
# 1. Download ENSEMBL data for human
pgatk ensembl-downloader -t 9606 -o ensembl_human

# 2. Extract transcript sequences (requires gffread)
gffread -F -w ensembl_human/transcripts.fa \
    -g ensembl_human/genome.fa \
    ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz

# 3. Translate variants to protein sequences
pgatk vcf-to-proteindb \
    --vcf ensembl_human/homo_sapiens_incl_consequences.vcf.gz \
    --input_fasta ensembl_human/transcripts.fa \
    --gene_annotations_gtf ensembl_human/Homo_sapiens.GRCh38.*.gtf.gz \
    --output_proteindb variant_proteins.fa

# 4. Generate target-decoy database
pgatk generate-decoy \
    --input variant_proteins.fa \
    --output target_decoy.fa \
    --method decoypyrat
```

## Commands

### Data Downloaders

| Command | Description |
|---------|-------------|
| `ensembl-downloader` | Download ENSEMBL reference data (GTF, FASTA, VCF) for any species by taxonomy ID |
| `ncbi-downloader` | Download NCBI RefSeq annotations and ClinVar VCF |
| `cosmic-downloader` | Download COSMIC somatic mutation data (requires account) |
| `cbioportal-downloader` | Download cBioPortal cancer genomics studies |

### Variant-to-Protein Translation

| Command | Description |
|---------|-------------|
| `vcf-to-proteindb` | Translate VCF variants (ENSEMBL, gnomAD, patient WES/WGS) to protein sequences |
| `clinvar-to-proteindb` | Translate ClinVar clinical variants (no VEP required) |
| `cosmic-to-proteindb` | Translate COSMIC somatic mutations, with optional tissue-type splitting |
| `cbioportal-to-proteindb` | Translate cBioPortal study mutations to protein sequences |

### Sequence Translation

| Command | Description |
|---------|-------------|
| `dnaseq-to-proteindb` | Translate DNA sequences with biotype filtering, multi-frame ORFs, and expression thresholds |
| `threeframe-translation` | Three-frame translation of transcript sequences |

### Database Processing

| Command | Description |
|---------|-------------|
| `generate-decoy` | Generate decoy sequences (methods: `decoypyrat`, `protein-reverse`, `protein-shuffle`, `pgdbdeep`) |
| `ensembl-check` | Validate protein database -- filter short sequences, handle stop codons |

### Post-Processing

| Command | Description |
|---------|-------------|
| `digest-mutant-protein` | In silico digest of variant proteins, filter against canonical proteome to extract unique peptides |
| `map-peptide2genome` | Map identified peptides to genomic coordinates (GFF3 output) |
| `spectrumai` | Inspect MS2 spectra of peptide identifications |
| `blast_get_position` | BLAST peptides against a reference database |

## Supported Variant Sources

| Source | Command | Description |
|--------|---------|-------------|
| ENSEMBL | `vcf-to-proteindb` | Population variants (SNPs, indels) for any ENSEMBL species |
| gnomAD | `vcf-to-proteindb` | Ancestry-stratified population variants (AF_afr, AF_eas, AF_nfe, etc.) |
| ClinVar | `clinvar-to-proteindb` | Clinically annotated pathogenic/benign variants |
| COSMIC | `cosmic-to-proteindb` | Somatic cancer mutations, per tissue type or cell line |
| cBioPortal | `cbioportal-to-proteindb` | Cancer study mutations from TCGA, METABRIC, etc. |
| Custom VCF | `vcf-to-proteindb` | Patient WGS/WES variants from any variant caller (GATK, Strelka, MuTect2) |

## Use Cases

Detailed end-to-end workflows are available in [docs/use-cases.md](docs/use-cases.md):

1. **Cell-type specific non-canonical peptide discovery** -- Reproduce the analysis from Umer et al. 2022
2. **Human variant protein database** -- Standard ENSEMBL-based variant proteogenomics
3. **Population-specific databases** -- gnomAD ancestry-stratified variant databases
4. **ClinVar clinical variants** -- Clinical variant detection at the protein level
5. **Cancer proteogenomics** -- COSMIC, cBioPortal, and patient-specific tumor databases
6. **Novel ORF and micropeptide discovery** -- lncRNA, pseudogene, and alternative ORF translation
7. **Genome annotation refinement** -- Six-frame translation and peptide-to-genome mapping
8. **Metaproteomics** -- Six-frame translation of metagenome assemblies
9. **Long-read transcriptomics** -- Isoform-resolved protein databases from PacBio/ONT data
10. **Plant and non-model organisms** -- Proteogenomics for any ENSEMBL species

## Project Structure

```
pgatk/
├── commands/           # CLI command definitions (Click)
├── ensembl/            # ENSEMBL data download and VCF translation
├── cgenomes/           # COSMIC and cBioPortal handling
├── clinvar/            # ClinVar variant translation
├── proteogenomics/     # Spectral validation tools
├── proteomics/         # Protein database utilities (decoy generation)
├── db/                 # Peptide digestion and genome mapping
├── config/             # YAML configuration files
└── toolbox/            # Shared utilities
```

## Full Documentation

[https://pgatk.quantms.org](https://pgatk.quantms.org)

## Cite

If you use pgatk in your research, please cite:

> Husen M Umer, Enrique Audain, Yafeng Zhu, Julianus Pfeuffer, Timo Sachsenberg, Janne Lehtiö, Rui M Branca, Yasset Perez-Riverol.
> **Generation of ENSEMBL-based proteogenomics databases boosts the identification of non-canonical peptides.**
> *Bioinformatics*, Volume 38, Issue 5, 1 March 2022, Pages 1470--1472.
> [https://doi.org/10.1093/bioinformatics/btab838](https://doi.org/10.1093/bioinformatics/btab838)

## Contributing

```bash
git clone https://github.com/bigbio/pgatk.git
cd pgatk
pip install -e ".[dev]"
pytest
```

## License

Apache License 2.0
