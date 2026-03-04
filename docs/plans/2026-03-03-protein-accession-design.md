# Protein Accession and FASTA Header Design

Issue: https://github.com/bigbio/pgatk/issues/18
Branch: `feature/protein-accession-design`
Date: 2026-03-03

## Problem

Current pgatk FASTA headers are inconsistent across variant sources (VCF, COSMIC, ClinVar) and incompatible with major search engines like SearchGUI, which cannot parse ENSEMBL-style IDs.

## Design

### Two protein categories, two prefix strategies

| Category | Prefix | Accession | Description |
|----------|--------|-----------|-------------|
| Canonical (reference) | Keep original (`sp\|`, `tr\|`, `ensp\|`) | Original accession | Untouched from source database |
| Variant (mutated) | `pgvar\|` | `{TRANSCRIPT_ID}-{INDEX}` | pgatk-generated variant protein |

### Variant header format

```
>pgvar|{TRANSCRIPT_ID}-{INDEX}|{GENE_SYMBOL} {key=value metadata}
```

**Fields:**

- `pgvar` -- database tag identifying pgatk-generated variant proteins.
- `{TRANSCRIPT_ID}-{INDEX}` -- accession composed of parent transcript ID and a dash-separated 1-based index (per transcript, per run). Mirrors UniProt isoform convention (`P12345-2`).
- `{GENE_SYMBOL}` -- gene name, first token after the second pipe.
- Metadata key=value pairs in the description field:

| Key | Description | Example |
|-----|-------------|---------|
| `VariantSource` | Origin database | `COSMIC`, `ClinVar`, `gnomAD`, `dbSNP` |
| `GenomicCoord` | `chr:pos:ref:alt` | `12:25245347:C:G` |
| `AAChange` | HGVS protein notation | `p.G13R` |
| `MutationType` | SO term or short label | `missense_variant` |
| `dbSNP` | rsID if available | `rs121913529` |
| `ORF` | Reading frame number (only when multi-ORF) | `1`, `2`, `3` |

### Examples

```fasta
# Canonical proteins -- untouched from source databases
>sp|P01112|RASH_HUMAN GTPase HRas OS=Homo sapiens
>ensp|ENSP00000309845|BRCA1

# Variant proteins -- unified pgvar| prefix regardless of source
>pgvar|ENST00000311189-1|HRAS VariantSource=COSMIC AAChange=p.G13R MutationType=missense_variant GenomicCoord=12:25245347:C:G
>pgvar|ENST00000311189-2|HRAS VariantSource=COSMIC AAChange=p.Q61L MutationType=missense_variant GenomicCoord=12:25245350:A:T
>pgvar|ENST00000357654-1|BRCA1 VariantSource=ClinVar AAChange=p.R1699Q MutationType=missense_variant GenomicCoord=17:43094464:G:A dbSNP=rs41293455

# Multiple ORFs -- each ORF gets its own index, ORF number in metadata
>pgvar|ENST00000311189-3|HRAS VariantSource=COSMIC AAChange=p.G13R ORF=1
>pgvar|ENST00000311189-4|HRAS VariantSource=COSMIC AAChange=p.G13R ORF=2
>pgvar|ENST00000311189-5|HRAS VariantSource=COSMIC AAChange=p.G13R ORF=3
```

### Indexing logic

The `-{index}` is per-transcript, per-file generation run:

- First variant on `ENST00000311189` gets `-1`
- Second variant on same transcript gets `-2`
- Multi-ORF outputs each consume an index (3 ORFs = 3 indices)
- First variant on a different transcript resets to `-1`

### Search engine compatibility

| Engine | Compatible | Notes |
|--------|-----------|-------|
| SearchGUI / PeptideShaker | Yes | Matches UniProt-like `db\|acc\|name` pattern |
| MaxQuant | Yes | Default UniProt parse rule works |
| MSFragger / FragPipe | Yes | Reads full header, splits on first whitespace |
| Comet | Yes | Parses `>db\|acc\|` natively |
| DIA-NN | Yes | Follows UniProt-style parsing |
| Proteome Discoverer | Yes | Supports pipe-delimited headers |

### Files to modify

| File | Change |
|------|--------|
| `pgatk/ensembl/ensembl.py` | Refactor `vcf_to_proteindb()` header construction (lines 661-664) |
| `pgatk/clinvar/clinvar_service.py` | Refactor header construction (lines 554-560) |
| `pgatk/cgenomes/cgenomes_proteindb.py` | Refactor COSMIC header (line 317), cBioPortal header |
| `pgatk/toolbox/vcf_utils.py` | Update `write_output()` to handle new format cleanly |
| `pgatk/config/` | Add constants for `PGVAR_PREFIX`, metadata keys |

### Design decisions

1. **Dash separator** (`-`) between transcript and index, consistent with UniProt isoform convention.
2. **No ORF suffix in accession** -- ORF number is metadata (`ORF=N`), each ORF gets its own index.
3. **Canonical proteins are pass-through** -- pgatk does not reformat existing database headers.
4. **Unified format across all sources** -- COSMIC, ClinVar, VCF variants all use `pgvar|` regardless of origin.
5. **Key=value metadata** in description field for structured downstream parsing.
