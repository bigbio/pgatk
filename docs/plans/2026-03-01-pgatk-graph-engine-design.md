# pgatk Evolution: Graph-Based Transcript Modeling & Feature Parity

**Date:** 2026-03-01
**Status:** Approved
**Authors:** Yasset Perez-Riverol, Claude (AI assistant)

---

## 1. Executive Summary

This document describes the research and implementation plan to evolve pgatk from a linear, one-variant-at-a-time proteogenomics database generator into a graph-based transcript modeling engine capable of multi-variant co-occurrence, gene fusions, RNA editing, and circular RNA support. The plan includes a prerequisite infrastructure cleanup phase and a phased feature rollout with performance optimization.

---

## 2. Current State Assessment

### What pgatk does well

- Only Python tool with native one-command downloaders for ENSEMBL, COSMIC, and cBioPortal
- Integrated decoy generation (4 methods: protein-reverse, protein-shuffle, DecoyPyrat, pgdbdeep)
- Spectrum-level validation via SpectrumAI (unique among database generators)
- Aho-Corasick exact matching + sliding-window mismatch search (BlastGetPosition)
- Published benchmark: 43,501 non-canonical peptides across 64 human cell lines (Umer et al. 2022)
- Pangenome proteogenomics demonstrated (Wang et al. 2024 preprint): 4,991 novel peptides from HPRC

### Critical gaps vs. the field (2024-2026)

| Gap | Impact | Solved by |
|---|---|---|
| No multi-variant co-occurrence | Transcript with N SNPs produces 1 sequence instead of up to 2^N | moPepGen (DAG), vcf2prot (phased VCF) |
| No gene fusions | Major cancer event type missed entirely | moPepGen |
| No RNA editing | A-to-I editing creates novel peptides invisible to DNA-only | moPepGen |
| No circular RNA | Emerging class of translated non-coding RNA | moPepGen |
| No ClinVar/NCBI support | Open issue #24 since 2019 | Manual workaround only |
| pandas iterrows in VCF processing | ~100x slower than vectorized operations | Internal bottleneck |
| No transcript feature caching | Redundant gffutils DB queries per variant | Internal bottleneck |

### Competitive landscape

**moPepGen** (Nature Biotechnology, June 2025): Graph-based (DAG per transcript), handles all event types, ~4x more non-canonical peptides than prior methods. The algorithmic benchmark to match.

**vcf2prot** (ElAbd et al. 2022, bioRxiv): Rust-based Sequence Intermediate Representation (SIR) approach. ~1000x faster than PrecisionProDB. Processes 99,254 variants across 8,192 patients in ~11 minutes. Demonstrates the performance ceiling achievable with systems-level optimization. Limited to phased VCFs with BCFtools/csq annotation.

**PG2** (J. Proteome Research 2023): Snakemake pipeline integrating genome + transcriptome. Handles splicing and fusions but requires heavy infrastructure.

**pXg** (MCP 2024): RNA-Seq + proteomics + de novo sequencing for immunopeptidomics.

**NeoDisc** (Nature Biotechnology 2024): End-to-end clinical neoantigen pipeline.

**PepCentric** (bioRxiv 2025): Repository-scale validation against 2.3B spectra. Complementary to database generators.

### Key code quality findings

- `get_altseq()` (ensembl.py): 60-line heart of variant processing, handles SNP/ins/del but only one variant at a time
- `vcf_to_proteindb()` (ensembl.py): ~260 lines, main bottleneck is `vcf_reader.iterrows()` + per-row gffutils DB lookups
- No `pathlib`, no `dataclasses`, inconsistent type hints
- `db/` standalone scripts use `sys.argv` at import time (not importable as modules)
- Variable shadowing bug in `protein_database_decoy.py` (`decoy_sequence` variable shadows the imported `pyteomics.fasta.decoy_sequence` function)
- Broad `except Exception` in several places suppresses unexpected errors
- Mixed `print()` and logger usage for error reporting
- Repetitive 3-layer config fallback pattern across all service classes

---

## 3. Architecture: TranscriptGraph Engine

### Core Concept

Each transcript is modeled as a directed acyclic graph (DAG):

- **Nodes** represent sequence segments (reference sequence between variant positions)
- **Edges** represent either the reference path or an alternative (variant) path
- **Traversal** of all root-to-leaf paths generates all combinatorial protein sequences

```
Reference:  ──[seg1]──[seg2]──[seg3]──[seg4]──
                        │                │
Variant A:              └─[alt_A]─┘      │
Variant B:                               └─[alt_B]─┘

Paths: ref, A only, B only, A+B  →  4 protein sequences
```

### Module Structure

```
pgatk/
  graph/                          # NEW module
    __init__.py
    transcript_graph.py           # TranscriptGraph class (DAG builder + traverser)
    variant_nodes.py              # Node types: SNP, Insertion, Deletion, Fusion, RNAEdit, CircRNA
    graph_enumerator.py           # Path enumeration with combinatorial explosion controls
    event_parsers/                # Input parsers for each event type
      __init__.py
      vcf_parser.py              # VCF → graph edges (replaces current get_altseq)
      fusion_parser.py           # Gene fusion calls → cross-transcript edges
      rnaedit_parser.py          # RNA editing sites → substitution edges
      circrna_parser.py          # Back-splice junctions → circular edges
      clinvar_parser.py          # ClinVar VCF → graph edges (issue #24)
    filters/
      __init__.py
      allele_frequency.py        # AF-based pruning
      consequence.py             # VEP consequence filtering
      expression.py              # Expression-level filtering
      max_variants.py            # Cap combinatorial explosion
```

### Key Classes

```python
@dataclass
class VariantNode:
    """Base class for all genomic events."""
    position: int          # CDS position
    ref_allele: str
    alt_allele: str
    event_type: EventType  # SNP, INS, DEL, FUSION, RNA_EDIT, CIRC_RNA
    metadata: dict         # AF, consequence, source, etc.


class TranscriptGraph:
    """DAG representing a transcript with all its variants."""

    def __init__(self, transcript_id: str, reference_seq: str, strand: str) -> None: ...
    def add_variant(self, variant: VariantNode) -> None: ...
    def add_fusion(self, partner_graph: 'TranscriptGraph', breakpoints: tuple) -> None: ...
    def enumerate_paths(self, max_paths: int = 1000) -> Iterator[ProteinSequence]: ...
    def to_fasta(self, output: TextIO) -> int: ...  # returns count written
```

### Combinatorial Explosion Controls

With N variants, naive enumeration produces 2^N sequences. Controls:

1. **Max variants per transcript** (default: 10) — skip transcripts exceeding this
2. **Max paths per transcript** (default: 1000) — stop enumeration after limit
3. **Allele frequency pruning** — only include variants above a threshold
4. **Phasing support** — if phased VCF, only enumerate haplotype-consistent paths
5. **Consequence filtering** — only include protein-altering consequences

### Integration Strategy

The graph engine sits alongside existing modules (no replacement):

- **Downloaders** remain unchanged
- **New CLI commands** (`graph-vcf-to-proteindb`, etc.) use the graph engine
- **Existing commands** continue working for backward compatibility
- **Shared infrastructure** (`toolbox/`, `config/`) is reused

---

## 4. Phased Implementation Plan

### Phase 0 — Infrastructure & Code Quality

**Goal:** Fix existing issues so the graph engine is built on solid ground.

**Deliverables:**

1. Fix variable shadowing bug in `protein_database_decoy.py`
2. Convert `db/digest_mutant_protein.py` and `db/map_peptide2genome.py` from `sys.argv`/`getopt` import-time execution to proper Click CLI commands
3. Add `conftest.py` with proper fixtures; remove relative path dependencies in tests
4. Replace `print()` error reporting with proper logger usage across all modules
5. Add type hints to all public APIs
6. Use `dataclasses` for data models (`SNP` in `cgenomes/models.py`, new models)
7. Use `pathlib.Path` for file operations (replace string concatenation)
8. Replace broad `except Exception` with specific exception handling
9. Standardize the 3-layer config fallback pattern (DRY up repetitive `get_*_parameters()` methods)
10. Expand unit test coverage for core functions: `get_altseq()`, `get_mut_pro_seq()`, `revswitch()`
11. Replace `pathos` with `concurrent.futures.ProcessPoolExecutor` (stdlib)

---

### Phase 1 — Graph Core + SNP/Indel Support

**Goal:** Working graph engine for basic variants, producing a superset of current output.

**Deliverables:**

1. `TranscriptGraph` class with add/enumerate/serialize
2. `VariantNode` dataclass hierarchy (SNP, Insertion, Deletion)
3. `vcf_parser.py` — reads VCF (VEP-annotated or unannotated) into graph edges
4. `graph_enumerator.py` — path enumeration with explosion controls
5. New CLI command: `graph-vcf-to-proteindb`
6. Transcript feature caching (eliminate redundant gffutils queries)
7. Vectorized VCF grouping (replace `iterrows()`)
8. Unit tests for graph operations + integration tests
9. Benchmark against current `vcf_to_proteindb()` on test VCFs

**Validation:** Graph engine must produce a superset of sequences from the linear approach.

---

### Phase 2 — Cancer Genomics + ClinVar

**Goal:** Bring cancer mutation sources into the graph model and add ClinVar (issue #24).

**Deliverables:**

1. `clinvar_parser.py` — parse ClinVar VCF into graph edges
2. Refactor `CancerGenomesService.get_mut_pro_seq()` to produce `VariantNode` objects
3. Graph-aware COSMIC processing (multiple mutations on same gene → combinatorial)
4. Graph-aware cBioPortal processing
5. New CLI commands: `clinvar-to-proteindb`, `graph-cosmic-to-proteindb`
6. Tests with real COSMIC data (multiple mutations per transcript)

---

### Phase 3 — Novel Event Types: Fusions, RNA Editing, Circular RNA

**Goal:** Add the event types that differentiate moPepGen from other tools.

**Deliverables:**

1. `fusion_parser.py` — parse gene fusion calls (STAR-Fusion, Arriba, FusionCatcher formats)
2. `rnaedit_parser.py` — parse RNA editing databases (REDIportal, DARNED) or user sites
3. `circrna_parser.py` — parse back-splice junction files (CIRCexplorer, CIRI)
4. Extended `TranscriptGraph` for cross-transcript fusions and circular paths
5. New CLI commands: `fusion-to-proteindb`, `rnaedit-to-proteindb`, `circrna-to-proteindb`
6. Combined command: `multi-event-to-proteindb` (merges all event types per transcript)

---

### Phase 4 — Performance Optimization

**Goal:** Make the graph engine fast enough for genome-scale datasets.

**Deliverables:**

1. Profile Phases 1-3 on real-world datasets (gnomAD, TCGA)
2. Identify hot paths via `cProfile` / `py-spy`
3. Optimize graph traversal (topological sort, memoization of shared subpaths)
4. Replace remaining `pathos` usage with `concurrent.futures`
5. Optional: C extension (via `cffi` or `cython`) for innermost graph traversal — only if profiling shows it's the bottleneck
6. Benchmark against moPepGen and vcf2prot on same datasets

---

### Phase Dependencies

```
Phase 0 (Infrastructure & Quality)
    └── Phase 1 (Graph Core + SNP/Indel)
            ├── Phase 2 (Cancer + ClinVar)
            │       └── Phase 3 (Fusions, RNA Edit, CircRNA)
            │               └── Phase 4 (Performance)
            └── (ongoing quality improvements)
```

---

## 5. References

### Primary pgatk publication
- Umer HM, Audain E, Zhu Y, Pfeuffer J, Sachsenberg T, Lehtiö J, Branca RM, Perez-Riverol Y. Generation of ENSEMBL-based proteogenomics databases boosts the identification of non-canonical peptides. *Bioinformatics*. 2022;38(5):1470-1472. doi:10.1093/bioinformatics/btab838

### Pangenome proteogenomics
- Wang D, Bouwmeester R, Zheng P, Dai C, Sanchez A, Shu K, Bai M, Umer HM, Perez-Riverol Y. Proteogenomics analysis of human tissues using pangenomes. *bioRxiv*. 2024. doi:10.1101/2024.05.24.595489

### Competing tools
- moPepGen: Graph-based proteogenomic database generation. *Nature Biotechnology*. June 2025. doi:10.1038/s41587-025-02701-0
- vcf2prot: ElAbd H, Degenhardt F, Lenz TL, Franke A, Wendorff M. VCF2Prot: An efficient and parallel tool for generating personalized proteomes from VCF files. *bioRxiv*. 2022. doi:10.1101/2022.01.21.477084
- ProteomeGenerator2: *J. Proteome Research*. 2023. doi:10.1021/acs.jproteome.3c00005
- pXg: *Molecular & Cellular Proteomics*. 2024. doi:10.1016/j.mcpro.2024.100733
- NeoDisc: *Nature Biotechnology*. October 2024. doi:10.1038/s41587-024-02420-y
- PepCentric: *bioRxiv*. February 2025. doi:10.1101/2025.02.24.639867

### Related ecosystem
- quantms: Dai C et al. *Nature Methods*. 2024. doi:10.1038/s41592-024-02343-1
- PRIDE 2025: *Nucleic Acids Research*. 2025;53(D1):D543. doi:10.1093/nar/gkae1011
