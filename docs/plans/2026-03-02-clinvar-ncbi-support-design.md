# ClinVar/NCBI RefSeq Protein Database Generation — Design

**Goal:** Add support for generating variant protein databases from ClinVar VCF files using NCBI RefSeq GTF annotations, independent of VEP/Ensembl.

**Architecture:** New `pgatk/clinvar/` module with ClinVar-specific pipeline, chromosome mapper, and NCBI data downloader. Shared algorithms (`get_altseq`, `get_orfs_vcf`) are reused from existing code.

**Tech Stack:** Python, pysam (VCF), gffutils (GTF), pyBedTools (interval overlap), pyfaidx (FASTA), tqdm, PyYAML.

---

## 1. Architecture Overview

### New module: `pgatk/clinvar/`

```
pgatk/clinvar/
    __init__.py
    clinvar_service.py        # Main pipeline: ClinVar VCF + RefSeq GTF -> protein FASTA
    chromosome_mapper.py      # NC_000001.11 <-> 1 <-> chr1 mapping
    data_downloader.py        # Download NCBI RefSeq + ClinVar files
```

### New CLI commands

```
pgatk/commands/
    clinvar_to_proteindb.py   # clinvar-to-proteindb command
    ncbi_downloader.py        # ncbi-downloader command
```

### Shared utilities

Extract reusable VCF-to-protein functions into `pgatk/toolbox/vcf_utils.py`:
- `get_altseq()` — apply variant to transcript sequence
- `get_orfs_vcf()` — translate modified transcript to protein
- `three_frame_translation()` — 3-frame translation helper

These are currently in `pgatk/ensembl/ensembl.py` and will be imported by both the Ensembl and ClinVar pipelines.

### New configuration

```
pgatk/config/clinvar_config.yaml
```

---

## 2. Pipeline Data Flow

### Stage 1: Chromosome Mapping

Parse NCBI assembly report (`GRCh38_latest_assembly_report.txt`) to build bidirectional chromosome name maps:
- `NC_000001.11` (RefSeq accession in GTF)
- `1` (numeric, used in ClinVar VCF)
- `chr1` (UCSC-style)

The `chromosome_mapper.py` module exposes `map_chrom(name, target_convention)` for converting between conventions.

### Stage 2: GTF Indexing

Use gffutils to build a database from the RefSeq GTF:
- Filter to `protein_coding` biotype transcripts (`NM_*` and optionally `XM_*`)
- Index CDS features per transcript
- Build transcript-to-gene mapping from `gene_id` attribute (gene symbol in RefSeq)
- Extract transcript CDS sequences from protein FASTA (or reconstruct from genome FASTA)

### Stage 3: ClinVar VCF Processing

For each VCF record:
1. **Filter by CLNSIG:** Exclude Benign, Likely_benign, Benign/Likely_benign
2. **Filter by MC (molecular consequence):** Keep missense, nonsense, frameshift, indels, splice variants
3. **Map chromosome:** Convert VCF chromosome (numeric) to GTF chromosome (NC_) using chromosome mapper
4. **Find overlapping transcripts:** Use BedTools/interval tree to find transcripts whose CDS region overlaps the variant position
5. **Apply variant:** For each overlapping transcript, use `get_altseq()` to modify the transcript sequence
6. **Translate:** Use `get_orfs_vcf()` to translate modified sequence to protein
7. **Collect unique variant proteins**

### Stage 4: Output

Write FASTA with ClinVar metadata headers:
```
>clinvar|RCV000012345|NM_000123|BRCA1|Pathogenic|missense_variant VARIANT_DESCRIPTION
MPROTEINSEQUENCE...
```

---

## 3. NCBI Downloader

### Files downloaded

| File | Source | Purpose |
|------|--------|---------|
| `GRCh38_latest_genomic.gtf.gz` | RefSeq FTP | Gene annotations |
| `GRCh38_latest_protein.faa.gz` | RefSeq FTP | Protein sequences |
| `GRCh38_latest_assembly_report.txt` | RefSeq FTP | Chromosome name mapping |
| `clinvar.vcf.gz` + `.tbi` | ClinVar FTP | Variant calls |

### Implementation (`pgatk/clinvar/data_downloader.py`)

- Uses `urllib.request` for FTP/HTTPS downloads (no extra dependencies)
- Downloads to user-specified `--output-dir`
- Skips existing files unless `--force` is passed
- Parses assembly report to produce chromosome mapping JSON
- Progress reporting via `tqdm`

### CLI

```
pgatk ncbi-downloader --output-dir ./ncbi_data [--genome-build GRCh38] [--force]
```

### Chromosome Mapper (`pgatk/clinvar/chromosome_mapper.py`)

- Reads assembly report TSV (columns: Sequence-Name, Sequence-Role, Assigned-Molecule, GenBank-Accn, RefSeq-Accn, UCSC-style-name)
- Builds bidirectional maps: `NC_000001.11 <-> 1 <-> chr1`
- Exposes `map_chrom(name, target_convention)` for any format conversion
- Handles unknown chromosomes gracefully (returns input unchanged + warning)

---

## 4. Configuration

### `pgatk/config/clinvar_config.yaml`

```yaml
# ClinVar-specific defaults
clinical_significance_exclude:
  - "Benign"
  - "Likely_benign"
  - "Benign/Likely_benign"

# NCBI RefSeq biotypes to include
biotypes:
  - "protein_coding"

# Molecular consequences to include (from ClinVar MC field)
consequences:
  - "missense_variant"
  - "nonsense"
  - "frameshift_variant"
  - "inframe_insertion"
  - "inframe_deletion"
  - "stop_gained"
  - "stop_lost"
  - "start_lost"
  - "splice_donor_variant"
  - "splice_acceptor_variant"

# FTP URLs
ncbi_ftp_base: "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/"
clinvar_ftp_base: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
```

### CLI parameters for `clinvar-to-proteindb`

```
pgatk clinvar-to-proteindb \
  --vcf clinvar.vcf.gz \
  --gtf GRCh38_latest_genomic.gtf \
  --fasta GRCh38_latest_protein.faa \
  --assembly-report GRCh38_latest_assembly_report.txt \
  --output clinvar_proteins.fa \
  [--clnsig-exclude "Benign,Likely_benign"] \
  [--consequences "missense_variant,stop_gained,..."] \
  [--biotypes "protein_coding"] \
  [--af-field AF_EXAC]
```

All filtering parameters have sensible defaults from the config file and can be overridden via CLI.

---

## 5. Testing Strategy

### Test data (committed to repo)

- `testdata/clinvar/mini_clinvar.vcf` — ~20 hand-picked ClinVar records: SNPs (missense, nonsense, stop_lost), indels, Benign + Pathogenic entries, multiple chromosomes
- `testdata/clinvar/mini_refseq.gtf` — GTF subset with ~5 transcripts matching the VCF variants
- `testdata/clinvar/mini_refseq_protein.faa` — Matching protein FASTA for the 5 transcripts
- `testdata/clinvar/mini_assembly_report.txt` — Assembly report subset for tested chromosomes
- `testdata/clinvar/expected_output.fa` — Expected output protein FASTA for validation

### Unit tests (`pgatk/tests/test_clinvar/`)

1. **`test_chromosome_mapper.py`** — NC_ <-> numeric <-> chr mapping, unknown chromosome handling, assembly report parsing
2. **`test_clinvar_service.py`** — CLNSIG filtering, transcript overlap detection, variant application (reusing `get_altseq`/`get_orfs_vcf`), end-to-end mini pipeline
3. **`test_data_downloader.py`** — Assembly report parsing (no actual FTP downloads), chromosome map JSON generation

### Integration test

- **`test_clinvar_integration.py`** — Full pipeline run with mini test data, compare output to expected FASTA

### Reuse validation

Shared functions (`get_altseq`, `get_orfs_vcf`) are already tested in `test_ensembl_core.py`. ClinVar tests focus on new code paths: chromosome mapping, CLNSIG filtering, transcript overlap, and CLI wiring.

---

## Design Decisions

1. **New module vs. extending Ensembl:** ClinVar/RefSeq is fundamentally different from the Ensembl+VEP pipeline (no CSQ annotations, different chromosome naming, different GTF structure). A separate module avoids coupling and keeps both pipelines clean.

2. **Internal transcript mapping vs. VEP:** ClinVar VCF lacks transcript-level annotations. Rather than requiring users to run VEP first, we use BedTools interval overlap to find transcripts internally. This is simpler for users and keeps the pipeline self-contained.

3. **Shared algorithm extraction:** `get_altseq` and `get_orfs_vcf` are generic enough to work with any VCF+GTF combination. Extracting them to `toolbox/vcf_utils.py` enables reuse without code duplication.

4. **Assembly report for chromosome mapping:** NCBI's assembly report is the authoritative source for chromosome name mappings. Parsing it once is more reliable than hardcoding mappings.
