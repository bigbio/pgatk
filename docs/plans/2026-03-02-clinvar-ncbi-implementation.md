# ClinVar/NCBI RefSeq Protein Database Generation — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add CLI commands to generate variant protein databases from ClinVar VCF files using NCBI RefSeq GTF annotations, without requiring VEP.

**Architecture:** New `pgatk/clinvar/` module with chromosome mapper, ClinVar VCF→protein pipeline, NCBI data downloader. Shared VCF→protein functions (`get_altseq`, `get_orfs_vcf`, etc.) are extracted from `pgatk/ensembl/ensembl.py` into `pgatk/toolbox/vcf_utils.py` and imported by both pipelines.

**Tech Stack:** Python 3.9+, gffutils (GTF indexing), pybedtools (interval overlap), BioPython (Seq, SeqIO), pandas, Click (CLI), PyYAML (config), tqdm (progress), urllib (downloads).

---

## Task 1: Extract Shared VCF→Protein Functions to `pgatk/toolbox/vcf_utils.py`

These static methods in `pgatk/ensembl/ensembl.py` are generic (no Ensembl-specific logic) and will be reused by the ClinVar pipeline. Extract them to a shared module.

**Files:**
- Create: `pgatk/toolbox/vcf_utils.py`
- Modify: `pgatk/ensembl/ensembl.py`
- Test: `pgatk/tests/test_vcf_utils.py`

**Step 1: Write the failing test**

Create `pgatk/tests/test_vcf_utils.py`:

```python
"""Tests for shared VCF utility functions extracted from EnsemblDataService."""

from Bio.Seq import Seq

from pgatk.toolbox.vcf_utils import check_overlap, get_altseq, get_orfs_vcf, write_output


class TestCheckOverlap:
    """Tests for check_overlap()."""

    def test_variant_fully_within_feature(self):
        features = [[100, 200, 'CDS']]
        assert check_overlap(150, 160, features) is True

    def test_variant_outside_feature(self):
        features = [[100, 200, 'CDS']]
        assert check_overlap(50, 60, features) is False

    def test_variant_start_neg1_always_overlaps(self):
        assert check_overlap(-1, 10, [[100, 200, 'CDS']]) is True

    def test_variant_partial_overlap_start(self):
        features = [[100, 200, 'CDS']]
        assert check_overlap(90, 110, features) is True


class TestGetAltseq:
    """Tests for get_altseq()."""

    def test_snp_plus_strand(self):
        ref_seq = Seq("ATGCCCGGG")
        features = [[10, 18, 'exon']]
        ref, alt = get_altseq(ref_seq, Seq("C"), Seq("T"), 13, '+', features)
        assert len(alt) > 0
        assert ref != alt

    def test_returns_tuple(self):
        ref_seq = Seq("ATGCCC")
        features = [[10, 15, 'exon']]
        result = get_altseq(ref_seq, Seq("C"), Seq("T"), 12, '+', features)
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestGetOrfsVcf:
    """Tests for get_orfs_vcf()."""

    def test_single_orf(self):
        ref = Seq("ATGCCCGGG")
        alt = Seq("ATGTCCGGG")
        ref_orfs, alt_orfs = get_orfs_vcf(ref, alt, 1, num_orfs=1)
        assert len(ref_orfs) == 1
        assert len(alt_orfs) == 1

    def test_three_orfs(self):
        ref = Seq("ATGCCCGGG")
        alt = Seq("ATGTCCGGG")
        ref_orfs, alt_orfs = get_orfs_vcf(ref, alt, 1, num_orfs=3)
        assert len(ref_orfs) == 3
        assert len(alt_orfs) == 3


class TestWriteOutput:
    """Tests for write_output()."""

    def test_writes_fasta(self, tmp_path):
        out_file = tmp_path / "test.fa"
        with open(out_file, 'w') as f:
            write_output("seq1", "description", ["MPKF"], f)
        content = out_file.read_text()
        assert ">seq1 description" in content
        assert "MPKF" in content
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_vcf_utils.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgatk.toolbox.vcf_utils'`

**Step 3: Create `pgatk/toolbox/vcf_utils.py`**

Extract these 4 static methods from `pgatk/ensembl/ensembl.py` (lines 142-166, 168-227, 280-296, 826-849) as standalone functions:

```python
"""Shared VCF-to-protein utility functions.

Extracted from pgatk.ensembl.ensembl.EnsemblDataService so that both the
Ensembl and ClinVar pipelines can reuse the same core algorithms.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

from Bio.Seq import Seq

logger = logging.getLogger(__name__)


def check_overlap(var_start: int, var_end: int, features_info: Optional[list] = None) -> bool:
    """Return True when the variant overlaps any of the features.

    :param var_start: Variant start position (genomic). -1 means always overlap.
    :param var_end: Variant end position (genomic).
    :param features_info: List of [start, end, type] for each feature.
    """
    if features_info is None:
        features_info = [[0, 1, 'type']]
    if var_start == -1:
        return True
    for feature_pos in features_info:
        pep_start = feature_pos[0]
        pep_end = feature_pos[1]
        if var_start <= pep_start <= var_end:
            return True
        elif var_start <= pep_end <= var_end:
            return True
        elif pep_start <= var_start and pep_end >= var_end:
            return True
    return False


def get_altseq(
    ref_seq: Seq,
    ref_allele: Seq,
    var_allele: Seq,
    var_pos: int,
    strand: str,
    features_info: list,
    cds_info: Optional[list] = None,
) -> tuple:
    """Modify a reference sequence based on a variant allele.

    Handles strand orientation, CDS trimming, and exon-based coordinate
    calculation.  Works identically with Ensembl and RefSeq annotations.

    :param ref_seq: Full transcript sequence (all exons concatenated).
    :param ref_allele: Reference allele (Seq object).
    :param var_allele: Alternate allele (Seq object).
    :param var_pos: Genomic position of the variant.
    :param strand: '+' or '-'.
    :param features_info: List of [start, end, type] per exon/CDS.
    :param cds_info: Optional [cds_start, cds_end] for CDS trimming.
    :return: (ref_coding_seq, alt_coding_seq) tuple.
    """
    if cds_info is None:
        cds_info = []
    alt_seq = ""
    if len(cds_info) == 2:
        start_coding_index = cds_info[0] - 1
        stop_coding_index = cds_info[1]
    else:
        start_coding_index = 0
        total_len = 0
        for x in features_info:
            total_len += x[1] - x[0] + 1
        stop_coding_index = total_len

    if strand == '-':
        ref_seq = ref_seq[::-1]
        ref_allele = ref_allele.complement()
        var_allele = var_allele.complement()

    if strand == '-' and len(cds_info) == 2:
        n = len(ref_seq)
        ref_seq = ref_seq[n - stop_coding_index:n - start_coding_index]
    else:
        ref_seq = ref_seq[start_coding_index:stop_coding_index]

    nc_index = 0
    if len(ref_allele) == len(var_allele) or ref_allele[0] == var_allele[0]:
        for feature in features_info:
            if var_pos in range(feature[0], feature[1] + 1):
                var_index_in_cds = nc_index + (var_pos - feature[0])
                c = len(ref_allele)
                alt_seq = ref_seq[0:var_index_in_cds] + var_allele + ref_seq[var_index_in_cds + c::]
                if strand == '-':
                    return ref_seq[::-1], alt_seq[::-1]
                else:
                    return ref_seq, alt_seq
            nc_index += (feature[1] - feature[0] + 1)

    return ref_seq, alt_seq


def get_orfs_vcf(
    ref_seq: Seq,
    alt_seq: Seq,
    translation_table: int,
    num_orfs: int = 1,
) -> tuple[list, list]:
    """Translate reference and alternate sequences into ORFs.

    :param ref_seq: Reference coding sequence.
    :param alt_seq: Alternate coding sequence with variant applied.
    :param translation_table: NCBI translation table number.
    :param num_orfs: Number of reading frames to generate (1 or 3).
    :return: (ref_orfs, alt_orfs) lists of translated protein sequences.
    """
    ref_orfs = []
    alt_orfs = []
    for n in range(0, num_orfs):
        ref_orfs.append(ref_seq[n::].translate(translation_table))
        alt_orfs.append(alt_seq[n::].translate(translation_table))
    return ref_orfs, alt_orfs


def write_output(
    seq_id: str,
    desc: str,
    seqs: list,
    prots_fn: Any,
    seqs_filter: Optional[list] = None,
) -> None:
    """Write ORFs to a FASTA output file.

    :param seq_id: Sequence accession/ID for the FASTA header.
    :param desc: Description string for the FASTA header.
    :param seqs: List of protein sequences (ORFs) to write.
    :param prots_fn: Open file handle for writing.
    :param seqs_filter: If provided, skip ORFs that match reference (used to
        write only variant-specific proteins).
    """
    if seqs_filter is None:
        seqs_filter = []
    write_i = False
    if len(seqs) > 1:
        write_i = True

    formatted_desc = " " + desc if desc else ""
    for i, orf in enumerate(seqs):
        if orf in seqs_filter:
            continue
        if write_i:
            prots_fn.write('>{}{}\n{}\n'.format(seq_id + "_" + str(i + 1), formatted_desc, orf))
        else:
            prots_fn.write('>{}{}\n{}\n'.format(seq_id, formatted_desc, orf))
```

**Step 4: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_vcf_utils.py -v`
Expected: All PASS

**Step 5: Update `pgatk/ensembl/ensembl.py` to import from shared module**

Replace the 4 static method bodies in `EnsemblDataService` with thin wrappers that delegate to the shared functions. This preserves backward compatibility for any code calling `EnsemblDataService.check_overlap(...)` etc.

In `pgatk/ensembl/ensembl.py`:
- Add import at the top: `from pgatk.toolbox.vcf_utils import check_overlap as _check_overlap, get_altseq as _get_altseq, get_orfs_vcf as _get_orfs_vcf, write_output as _write_output`
- Replace method bodies of `check_overlap` (lines 142-165), `get_altseq` (lines 168-227), `get_orfs_vcf` (lines 280-296), `write_output` (lines 826-849) with delegation calls:

```python
@staticmethod
def check_overlap(var_start, var_end, features_info=None):
    return _check_overlap(var_start, var_end, features_info)

@staticmethod
def get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info=None):
    return _get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info)

@staticmethod
def get_orfs_vcf(ref_seq, alt_seq, translation_table, num_orfs=1):
    return _get_orfs_vcf(ref_seq, alt_seq, translation_table, num_orfs)

@staticmethod
def write_output(seq_id, desc, seqs, prots_fn, seqs_filter=None):
    _write_output(seq_id, desc, seqs, prots_fn, seqs_filter)
```

**Step 6: Run ALL existing tests to verify no regressions**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v`
Expected: All 66+ tests PASS (existing tests still call `EnsemblDataService.get_altseq(...)` etc., which now delegate to the shared module)

**Step 7: Commit**

```bash
git add pgatk/toolbox/vcf_utils.py pgatk/tests/test_vcf_utils.py pgatk/ensembl/ensembl.py
git commit -m "refactor: extract shared VCF utility functions to toolbox/vcf_utils.py

Extract check_overlap, get_altseq, get_orfs_vcf, write_output from
EnsemblDataService to pgatk.toolbox.vcf_utils for reuse by ClinVar pipeline.
EnsemblDataService methods now delegate to shared functions."
```

---

## Task 2: Chromosome Mapper

Build the chromosome name mapper that converts between NC_ accessions (RefSeq GTF), numeric names (ClinVar VCF), and UCSC chr-prefixed names.

**Files:**
- Create: `pgatk/clinvar/__init__.py`
- Create: `pgatk/clinvar/chromosome_mapper.py`
- Create: `pgatk/testdata/clinvar/mini_assembly_report.txt`
- Test: `pgatk/tests/test_clinvar/__init__.py`
- Test: `pgatk/tests/test_clinvar/test_chromosome_mapper.py`

**Step 1: Create test data**

Create `pgatk/testdata/clinvar/mini_assembly_report.txt` — a subset of the real NCBI assembly report format. The file has comment lines starting with `#`, then a header line, then TSV data. The columns we care about are: col 0 = Sequence-Name (e.g. "1"), col 6 = RefSeq-Accn (e.g. "NC_000001.11"), col 9 = UCSC-style-name (e.g. "chr1"):

```
# Assembly name:  GRCh38.p14
# Description:    Genome Reference Consortium Human Build 38 patch release 14
# Organism name:  Homo sapiens
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-loc/type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
2	assembled-molecule	2	Chromosome	CM000664.2	=	NC_000002.12	Primary Assembly	242193529	chr2
X	assembled-molecule	X	Chromosome	CM000685.2	=	NC_000023.11	Primary Assembly	156040895	chrX
Y	assembled-molecule	Y	Chromosome	CM000686.2	=	NC_000024.10	Primary Assembly	57227415	chrY
MT	assembled-molecule	MT	Mitochondria	J01415.2	=	NC_012920.1	non-nuclear	16569	chrM
```

**Step 2: Write the failing test**

Create `pgatk/tests/test_clinvar/__init__.py` (empty file).

Create `pgatk/tests/test_clinvar/test_chromosome_mapper.py`:

```python
"""Tests for pgatk.clinvar.chromosome_mapper."""

from pathlib import Path

import pytest

from pgatk.clinvar.chromosome_mapper import ChromosomeMapper


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"
ASSEMBLY_REPORT = TESTDATA_DIR / "mini_assembly_report.txt"


class TestChromosomeMapper:
    """Tests for ChromosomeMapper."""

    @pytest.fixture
    def mapper(self):
        return ChromosomeMapper.from_assembly_report(str(ASSEMBLY_REPORT))

    def test_numeric_to_refseq(self, mapper):
        assert mapper.map_chrom("1", "refseq") == "NC_000001.11"

    def test_refseq_to_numeric(self, mapper):
        assert mapper.map_chrom("NC_000001.11", "numeric") == "1"

    def test_numeric_to_ucsc(self, mapper):
        assert mapper.map_chrom("1", "ucsc") == "chr1"

    def test_ucsc_to_numeric(self, mapper):
        assert mapper.map_chrom("chr1", "numeric") == "1"

    def test_refseq_to_ucsc(self, mapper):
        assert mapper.map_chrom("NC_000001.11", "ucsc") == "chr1"

    def test_x_chromosome(self, mapper):
        assert mapper.map_chrom("X", "refseq") == "NC_000023.11"
        assert mapper.map_chrom("NC_000023.11", "numeric") == "X"

    def test_mt_chromosome(self, mapper):
        assert mapper.map_chrom("MT", "refseq") == "NC_012920.1"
        assert mapper.map_chrom("NC_012920.1", "ucsc") == "chrM"

    def test_unknown_chromosome_returns_input(self, mapper):
        assert mapper.map_chrom("UNKNOWN_CONTIG", "refseq") == "UNKNOWN_CONTIG"

    def test_invalid_convention_raises(self, mapper):
        with pytest.raises(ValueError, match="convention"):
            mapper.map_chrom("1", "invalid_convention")

    def test_chr_prefixed_input(self, mapper):
        """Input with chr prefix should be recognized as UCSC convention."""
        assert mapper.map_chrom("chrX", "numeric") == "X"
        assert mapper.map_chrom("chrX", "refseq") == "NC_000023.11"
```

**Step 3: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_chromosome_mapper.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgatk.clinvar'`

**Step 4: Implement chromosome mapper**

Create `pgatk/clinvar/__init__.py` (empty file).

Create `pgatk/clinvar/chromosome_mapper.py`:

```python
"""Bidirectional chromosome name mapping: RefSeq (NC_) <-> numeric <-> UCSC (chr).

Parses an NCBI assembly report to build the mapping tables.
"""
from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Valid target convention names
_VALID_CONVENTIONS = {"numeric", "refseq", "ucsc"}


class ChromosomeMapper:
    """Maps chromosome names between naming conventions.

    Supports three conventions:
    - 'numeric': e.g. '1', '2', 'X', 'Y', 'MT'
    - 'refseq': e.g. 'NC_000001.11'
    - 'ucsc': e.g. 'chr1', 'chrX', 'chrM'
    """

    def __init__(
        self,
        numeric_to_refseq: dict[str, str],
        numeric_to_ucsc: dict[str, str],
    ) -> None:
        self._num_to_ref = numeric_to_refseq
        self._num_to_ucsc = numeric_to_ucsc

        # Build reverse maps
        self._ref_to_num = {v: k for k, v in numeric_to_refseq.items()}
        self._ucsc_to_num = {v: k for k, v in numeric_to_ucsc.items()}

    @classmethod
    def from_assembly_report(cls, report_path: str) -> ChromosomeMapper:
        """Parse an NCBI assembly report file and build chromosome maps.

        The assembly report is a TSV with comment lines starting with '#'.
        Relevant columns (0-indexed):
        - 0: Sequence-Name (e.g. '1', 'X', 'MT')
        - 6: RefSeq-Accn (e.g. 'NC_000001.11')
        - 9: UCSC-style-name (e.g. 'chr1')
        """
        numeric_to_refseq: dict[str, str] = {}
        numeric_to_ucsc: dict[str, str] = {}

        with open(report_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                cols = line.split("\t")
                if len(cols) < 10:
                    continue
                seq_name = cols[0]       # numeric name
                refseq_accn = cols[6]    # NC_ accession
                ucsc_name = cols[9]      # chr-prefixed name

                if refseq_accn and refseq_accn != "na":
                    numeric_to_refseq[seq_name] = refseq_accn
                if ucsc_name and ucsc_name != "na":
                    numeric_to_ucsc[seq_name] = ucsc_name

        logger.info(
            "Loaded chromosome mapping: %d entries from %s",
            len(numeric_to_refseq),
            Path(report_path).name,
        )
        return cls(numeric_to_refseq, numeric_to_ucsc)

    def _to_numeric(self, name: str) -> str | None:
        """Convert any convention to numeric. Returns None if unknown."""
        # Already numeric?
        if name in self._num_to_ref or name in self._num_to_ucsc:
            return name
        # RefSeq?
        if name in self._ref_to_num:
            return self._ref_to_num[name]
        # UCSC?
        if name in self._ucsc_to_num:
            return self._ucsc_to_num[name]
        return None

    def map_chrom(self, name: str, target_convention: str) -> str:
        """Map a chromosome name to the target convention.

        :param name: Input chromosome name in any convention.
        :param target_convention: One of 'numeric', 'refseq', 'ucsc'.
        :return: Mapped name, or the original name if mapping is unknown.
        :raises ValueError: If target_convention is not valid.
        """
        if target_convention not in _VALID_CONVENTIONS:
            raise ValueError(
                f"Invalid convention '{target_convention}'. "
                f"Must be one of: {', '.join(sorted(_VALID_CONVENTIONS))}"
            )

        numeric = self._to_numeric(name)
        if numeric is None:
            logger.warning("Unknown chromosome name: %s — returning as-is", name)
            return name

        if target_convention == "numeric":
            return numeric
        elif target_convention == "refseq":
            return self._num_to_ref.get(numeric, name)
        elif target_convention == "ucsc":
            return self._num_to_ucsc.get(numeric, name)
        return name  # unreachable but satisfies linters
```

**Step 5: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_chromosome_mapper.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/clinvar/__init__.py pgatk/clinvar/chromosome_mapper.py \
       pgatk/testdata/clinvar/mini_assembly_report.txt \
       pgatk/tests/test_clinvar/__init__.py \
       pgatk/tests/test_clinvar/test_chromosome_mapper.py
git commit -m "feat: add ChromosomeMapper for NC_/numeric/UCSC chromosome name conversion

Parses NCBI assembly report to build bidirectional chromosome name maps.
Supports RefSeq (NC_000001.11), numeric (1), and UCSC (chr1) conventions."
```

---

## Task 3: ClinVar Config and Config Registry Integration

Add ClinVar-specific configuration defaults and register them in the config registry.

**Files:**
- Create: `pgatk/config/clinvar_config.yaml`
- Modify: `pgatk/config/registry.py`

**Step 1: Create `pgatk/config/clinvar_config.yaml`**

```yaml
clinvar_translation:
  translation_table: 1
  mito_translation_table: 2
  protein_prefix: "clinvar_"
  proteindb_output_file: "clinvar-peptide-database.fa"
  report_ref_seq: false
  num_orfs: 1

  # ClinVar-specific filtering
  clinical_significance_exclude:
    - "Benign"
    - "Likely_benign"
    - "Benign/Likely_benign"

  # NCBI RefSeq biotypes to include
  include_biotypes: "protein_coding"

  # Molecular consequences to include (from ClinVar MC field)
  include_consequences:
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

  # NCBI FTP URLs
  ncbi_refseq_ftp: "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/"
  clinvar_ftp: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"

  logger:
    formatters:
      DEBUG: "%(asctime)s [%(levelname)7s][%(name)48s][%(module)32s, %(lineno)4s] %(message)s"
      INFO: "%(asctime)s [%(levelname)7s][%(name)48s] %(message)s"
    loglevel: DEBUG
```

**Step 2: Register in `pgatk/config/registry.py`**

Add `"clinvar": "clinvar_config.yaml"` to the `COMMAND_CONFIGS` dict at line 16 of `pgatk/config/registry.py`:

```python
COMMAND_CONFIGS = {
    "ensembl_downloader": "ensembl_downloader_config.yaml",
    "ensembl_config": "ensembl_config.yaml",
    "cosmic": "cosmic_config.yaml",
    "cbioportal": "cbioportal_config.yaml",
    "protein_decoy": "protein_decoy.yaml",
    "clinvar": "clinvar_config.yaml",
}
```

**Step 3: Verify config loads**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -c "from pgatk.config.registry import load_config; c = load_config('clinvar'); print(c['clinvar_translation']['clinical_significance_exclude'])"`
Expected: `['Benign', 'Likely_benign', 'Benign/Likely_benign']`

**Step 4: Commit**

```bash
git add pgatk/config/clinvar_config.yaml pgatk/config/registry.py
git commit -m "feat: add ClinVar config with CLNSIG filtering and NCBI FTP URLs"
```

---

## Task 4: ClinVar Service — Core Pipeline

The main ClinVar VCF-to-protein pipeline. This is the largest task.

**Files:**
- Create: `pgatk/clinvar/clinvar_service.py`
- Create: `pgatk/testdata/clinvar/mini_clinvar.vcf`
- Create: `pgatk/testdata/clinvar/mini_refseq.gtf`
- Create: `pgatk/testdata/clinvar/mini_refseq_protein.faa`
- Test: `pgatk/tests/test_clinvar/test_clinvar_service.py`

**Step 1: Create test data files**

These are minimal hand-crafted files that exercise the key code paths.

Create `pgatk/testdata/clinvar/mini_clinvar.vcf`:
```
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant">
##INFO=<ID=MC,Number=.,Type=String,Description="Molecular consequence">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100015	rs001	C	T	.	.	CLNSIG=Pathogenic;GENEINFO=TESTGENE1:1234;MC=SO:0001583|missense_variant
1	100025	rs002	G	A	.	.	CLNSIG=Benign;GENEINFO=TESTGENE1:1234;MC=SO:0001583|missense_variant
2	200010	rs003	A	G	.	.	CLNSIG=Likely_pathogenic;GENEINFO=TESTGENE2:5678;MC=SO:0001587|stop_gained
```

Create `pgatk/testdata/clinvar/mini_refseq.gtf` — RefSeq-style GTF with NC_ chromosome accessions. Note: features must overlap the VCF variant positions after chromosome mapping.

```
NC_000001.11	BestRefSeq	transcript	100001	100090	.	+	.	gene_id "TESTGENE1"; transcript_id "NM_000001.1"; gene_biotype "protein_coding";
NC_000001.11	BestRefSeq	exon	100001	100030	.	+	.	gene_id "TESTGENE1"; transcript_id "NM_000001.1";
NC_000001.11	BestRefSeq	CDS	100001	100030	.	+	0	gene_id "TESTGENE1"; transcript_id "NM_000001.1";
NC_000001.11	BestRefSeq	exon	100051	100090	.	+	.	gene_id "TESTGENE1"; transcript_id "NM_000001.1";
NC_000001.11	BestRefSeq	CDS	100051	100090	.	+	2	gene_id "TESTGENE1"; transcript_id "NM_000001.1";
NC_000002.12	BestRefSeq	transcript	200001	200060	.	+	.	gene_id "TESTGENE2"; transcript_id "NM_000002.1"; gene_biotype "protein_coding";
NC_000002.12	BestRefSeq	exon	200001	200060	.	+	.	gene_id "TESTGENE2"; transcript_id "NM_000002.1";
NC_000002.12	BestRefSeq	CDS	200001	200060	.	+	0	gene_id "TESTGENE2"; transcript_id "NM_000002.1";
```

Create `pgatk/testdata/clinvar/mini_refseq_protein.faa` — transcript sequences for the transcripts referenced in the GTF. Sequences must be nucleotide, matching the CDS lengths:

```
>NM_000001.1 CDS=1-70
ATGCCCGGGAAATTTCCCGGGAAATTTCCCATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCC
>NM_000002.1 CDS=1-60
ATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAA
```

**Step 2: Write the failing test**

Create `pgatk/tests/test_clinvar/test_clinvar_service.py`:

```python
"""Tests for ClinVar VCF-to-protein pipeline."""

from pathlib import Path

import pytest

from pgatk.clinvar.clinvar_service import ClinVarService


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"


class TestClinSigFiltering:
    """Tests for CLNSIG filtering logic."""

    def test_pathogenic_passes(self):
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Pathogenic", exclude) is True

    def test_benign_excluded(self):
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Benign", exclude) is False

    def test_likely_benign_excluded(self):
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Likely_benign", exclude) is False

    def test_uncertain_passes(self):
        exclude = ["Benign", "Likely_benign", "Benign/Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("Uncertain_significance", exclude) is True

    def test_empty_clnsig_passes(self):
        exclude = ["Benign", "Likely_benign"]
        assert ClinVarService.passes_clnsig_filter("", exclude) is True


class TestMolecularConsequenceParser:
    """Tests for MC field parsing."""

    def test_parse_missense(self):
        mc_field = "SO:0001583|missense_variant"
        assert ClinVarService.parse_mc_consequence(mc_field) == "missense_variant"

    def test_parse_multiple_consequences(self):
        mc_field = "SO:0001583|missense_variant,SO:0001587|stop_gained"
        result = ClinVarService.parse_mc_consequences(mc_field)
        assert "missense_variant" in result
        assert "stop_gained" in result

    def test_parse_empty_mc(self):
        assert ClinVarService.parse_mc_consequence("") == ""


class TestGeneInfoParser:
    """Tests for GENEINFO field parsing."""

    def test_parse_single_gene(self):
        gene_symbol, gene_id = ClinVarService.parse_geneinfo("BRCA1:672")
        assert gene_symbol == "BRCA1"
        assert gene_id == "672"

    def test_parse_multi_gene(self):
        gene_symbol, gene_id = ClinVarService.parse_geneinfo("BRCA1:672|TP53:7157")
        assert gene_symbol == "BRCA1"  # first gene


class TestClinVarPipeline:
    """Integration test for ClinVar-to-proteindb pipeline."""

    def test_pipeline_produces_output(self, tmp_path):
        output_file = tmp_path / "output.fa"
        service = ClinVarService(
            vcf_file=str(TESTDATA_DIR / "mini_clinvar.vcf"),
            gtf_file=str(TESTDATA_DIR / "mini_refseq.gtf"),
            fasta_file=str(TESTDATA_DIR / "mini_refseq_protein.faa"),
            assembly_report=str(TESTDATA_DIR / "mini_assembly_report.txt"),
            output_file=str(output_file),
        )
        service.run()
        assert output_file.exists()
        content = output_file.read_text()
        # Should contain at least one protein entry
        assert ">" in content

    def test_benign_variants_excluded(self, tmp_path):
        output_file = tmp_path / "output.fa"
        service = ClinVarService(
            vcf_file=str(TESTDATA_DIR / "mini_clinvar.vcf"),
            gtf_file=str(TESTDATA_DIR / "mini_refseq.gtf"),
            fasta_file=str(TESTDATA_DIR / "mini_refseq_protein.faa"),
            assembly_report=str(TESTDATA_DIR / "mini_assembly_report.txt"),
            output_file=str(output_file),
        )
        service.run()
        content = output_file.read_text()
        # rs002 is Benign and should NOT appear
        assert "rs002" not in content
```

**Step 3: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_service.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'pgatk.clinvar.clinvar_service'`

**Step 4: Implement `pgatk/clinvar/clinvar_service.py`**

```python
"""ClinVar VCF-to-protein database pipeline.

Processes ClinVar VCF records against NCBI RefSeq GTF annotations to produce
a FASTA file of variant protein sequences. Does NOT require VEP — uses
BedTools interval overlap to find affected transcripts internally.
"""
from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Optional

import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pybedtools import BedTool

from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
from pgatk.toolbox.vcf_utils import check_overlap, get_altseq, get_orfs_vcf, write_output

logger = logging.getLogger(__name__)

# Default CLNSIG values to exclude
_DEFAULT_CLNSIG_EXCLUDE = ["Benign", "Likely_benign", "Benign/Likely_benign"]

# Default molecular consequences to include
_DEFAULT_CONSEQUENCES = [
    "missense_variant", "nonsense", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "stop_gained",
    "stop_lost", "start_lost", "splice_donor_variant", "splice_acceptor_variant",
]


class ClinVarService:
    """Generate variant protein databases from ClinVar VCF + RefSeq GTF."""

    def __init__(
        self,
        vcf_file: str,
        gtf_file: str,
        fasta_file: str,
        assembly_report: str,
        output_file: str,
        clnsig_exclude: Optional[list[str]] = None,
        consequences: Optional[list[str]] = None,
        translation_table: int = 1,
        mito_translation_table: int = 2,
        protein_prefix: str = "clinvar_",
        num_orfs: int = 1,
        report_ref_seq: bool = False,
    ) -> None:
        self._vcf_file = vcf_file
        self._gtf_file = gtf_file
        self._fasta_file = fasta_file
        self._assembly_report = assembly_report
        self._output_file = output_file
        self._clnsig_exclude = clnsig_exclude or _DEFAULT_CLNSIG_EXCLUDE
        self._consequences = consequences or _DEFAULT_CONSEQUENCES
        self._translation_table = translation_table
        self._mito_translation_table = mito_translation_table
        self._protein_prefix = protein_prefix
        self._num_orfs = num_orfs
        self._report_ref_seq = report_ref_seq

    # ------------------------------------------------------------------
    # Static helper methods for parsing ClinVar INFO fields
    # ------------------------------------------------------------------

    @staticmethod
    def passes_clnsig_filter(clnsig: str, exclude_list: list[str]) -> bool:
        """Return True if the clinical significance is NOT in the exclude list."""
        if not clnsig:
            return True
        return clnsig not in exclude_list

    @staticmethod
    def parse_mc_consequence(mc_field: str) -> str:
        """Parse first molecular consequence from the MC INFO field.

        MC format: 'SO:0001583|missense_variant'
        """
        if not mc_field:
            return ""
        first_mc = mc_field.split(",")[0]
        parts = first_mc.split("|")
        return parts[1] if len(parts) > 1 else ""

    @staticmethod
    def parse_mc_consequences(mc_field: str) -> list[str]:
        """Parse all molecular consequences from the MC field."""
        if not mc_field:
            return []
        consequences = []
        for entry in mc_field.split(","):
            parts = entry.strip().split("|")
            if len(parts) > 1:
                consequences.append(parts[1])
        return consequences

    @staticmethod
    def parse_geneinfo(geneinfo: str) -> tuple[str, str]:
        """Parse GENEINFO field. Format: 'SYMBOL:GENEID[|SYMBOL2:GENEID2]'.

        Returns (gene_symbol, gene_id) for the first gene.
        """
        if not geneinfo:
            return "", ""
        first_gene = geneinfo.split("|")[0]
        parts = first_gene.split(":")
        return (parts[0], parts[1]) if len(parts) == 2 else (geneinfo, "")

    @staticmethod
    def _get_info_field(info_str: str, key: str) -> str:
        """Extract a value from the VCF INFO column by key."""
        for field in info_str.split(";"):
            if field.startswith(key + "="):
                return field.split("=", 1)[1]
        return ""

    # ------------------------------------------------------------------
    # GTF parsing (reuses gffutils, same pattern as EnsemblDataService)
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_gtf(gtf_file: str) -> gffutils.FeatureDB:
        """Create or load a gffutils database from a GTF file."""
        db_file = str(Path(gtf_file).with_suffix(".db"))
        try:
            gffutils.create_db(
                gtf_file, db_file,
                merge_strategy="create_unique",
                keep_order=True,
                disable_infer_transcripts=True,
                disable_infer_genes=True,
                verbose=False,
                force=False,
            )
        except (ValueError, sqlite3.OperationalError):
            logger.info("GTF database already exists: %s", db_file)

        return gffutils.FeatureDB(db_file)

    @staticmethod
    def _get_features(
        db: gffutils.FeatureDB,
        feature_id: str,
        feature_types: Optional[list[str]] = None,
    ) -> tuple:
        """Retrieve chromosome, strand, and coding features for a transcript."""
        if feature_types is None:
            feature_types = ["exon"]
        try:
            feature = db[feature_id]
        except gffutils.exceptions.FeatureNotFoundError:
            try:
                feature = db[feature_id.split(".")[0]]
            except gffutils.exceptions.FeatureNotFoundError:
                logger.warning("Feature %s not found in GTF database", feature_id)
                return None, None, None

        coding_features = []
        for f in db.children(feature, featuretype=feature_types, order_by="end"):
            coding_features.append([f.start, f.end, f.featuretype])
        return feature.chrom, feature.strand, coding_features

    # ------------------------------------------------------------------
    # VCF reading (same pandas approach as EnsemblDataService)
    # ------------------------------------------------------------------

    @staticmethod
    def _read_vcf(vcf_file: str) -> tuple[list[str], pd.DataFrame]:
        """Read a VCF file into metadata headers + DataFrame."""
        metadata = []
        records = []
        with open(vcf_file, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#"):
                    metadata.append(line)
                    continue
                cols = line.split("\t")
                if len(cols) < 8:
                    continue
                records.append({
                    "CHROM": cols[0],
                    "POS": cols[1],
                    "ID": cols[2],
                    "REF": cols[3],
                    "ALT": cols[4],
                    "QUAL": cols[5],
                    "FILTER": cols[6],
                    "INFO": cols[7],
                })
        df = pd.DataFrame(records) if records else pd.DataFrame(
            columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        )
        return metadata, df

    # ------------------------------------------------------------------
    # Transcript overlap via BedTools
    # ------------------------------------------------------------------

    @staticmethod
    def _annotate_vcf_with_transcripts(
        vcf_file: str,
        gtf_file: str,
        chrom_mapper: ChromosomeMapper,
    ) -> dict[str, list[str]]:
        """Find overlapping transcripts for each VCF record using BedTools.

        Returns a dict mapping 'CHROM:POS:REF:ALT' -> [transcript_id, ...].
        Chromosome names are normalized to RefSeq (NC_) for GTF matching.
        """
        overlaps: dict[str, list[str]] = {}

        # Read VCF records and convert to BED format in RefSeq chromosome names
        vcf_bed_entries = []
        _, vcf_df = ClinVarService._read_vcf(vcf_file)

        for _, row in vcf_df.iterrows():
            chrom_refseq = chrom_mapper.map_chrom(str(row.CHROM), "refseq")
            start = int(row.POS) - 1  # BED is 0-based
            end = start + len(str(row.REF))
            key = f"{row.CHROM}:{row.POS}:{row.REF}:{row.ALT}"
            vcf_bed_entries.append(f"{chrom_refseq}\t{start}\t{end}\t{key}")

        if not vcf_bed_entries:
            return overlaps

        vcf_bed = BedTool("\n".join(vcf_bed_entries), from_string=True)

        # Intersect with GTF — extract CDS features only
        try:
            gtf_bed = BedTool(gtf_file)
            intersection = gtf_bed.intersect(vcf_bed, wo=True)
        except Exception as e:
            logger.error("BedTools intersection failed: %s", e)
            return overlaps

        for feature in intersection:
            fields = str(feature).strip().split("\t")
            # GTF feature type is at index 2
            if fields[2] != "CDS":
                continue
            # Parse transcript_id from GTF attributes (index 8)
            attrs = fields[8]
            transcript_id = ""
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    transcript_id = attr.split('"')[1] if '"' in attr else attr.split(" ")[1]
                    break
            if not transcript_id:
                continue

            # The key is in the last few columns (from the VCF BED name field)
            # BedTool intersect with wo appends the second file fields
            variant_key = fields[-2]  # the name column from the VCF BED
            overlaps.setdefault(variant_key, [])
            if transcript_id not in overlaps[variant_key]:
                overlaps[variant_key].append(transcript_id)

        return overlaps

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def run(self) -> str:
        """Execute the ClinVar-to-proteindb pipeline. Returns output path."""

        # 1. Load chromosome mapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(self._assembly_report)

        # 2. Parse GTF and index transcripts
        db = self._parse_gtf(self._gtf_file)

        # 3. Load transcript FASTA sequences
        transcripts_dict = SeqIO.index(self._fasta_file, "fasta",
                                        key_function=lambda h: h.split("|")[0].split(" ")[0])
        transcript_id_mapping = {k.split(".")[0]: k for k in transcripts_dict.keys()}

        # 4. Find overlapping transcripts for all VCF records
        variant_transcripts = self._annotate_vcf_with_transcripts(
            self._vcf_file, self._gtf_file, chrom_mapper,
        )

        # 5. Process VCF records
        _, vcf_df = self._read_vcf(self._vcf_file)

        stats = {
            "variants_processed": 0,
            "variants_excluded_clnsig": 0,
            "variants_no_overlap": 0,
            "variants_translated": 0,
            "transcript_not_found": 0,
        }

        with open(self._output_file, "w", encoding="utf-8") as out_fh:
            for _, record in vcf_df.iterrows():
                stats["variants_processed"] += 1

                # Validate REF/ALT
                ref = str(record.REF)
                alts = [a for a in str(record.ALT).split(",")
                        if a and all(c in "ACGT" for c in a)]
                if not all(c in "ACGT" for c in ref) or not alts:
                    continue

                # Filter by CLNSIG
                clnsig = self._get_info_field(str(record.INFO), "CLNSIG")
                if not self.passes_clnsig_filter(clnsig, self._clnsig_exclude):
                    stats["variants_excluded_clnsig"] += 1
                    continue

                # Translation table
                trans_table = self._translation_table
                chrom_str = str(record.CHROM).lstrip("chr").upper()
                if chrom_str in ("M", "MT"):
                    trans_table = self._mito_translation_table

                # Get overlapping transcripts
                variant_key = f"{record.CHROM}:{record.POS}:{record.REF}:{record.ALT}"
                transcript_ids = variant_transcripts.get(variant_key, [])
                if not transcript_ids:
                    stats["variants_no_overlap"] += 1
                    continue

                for transcript_id in transcript_ids:
                    # Resolve transcript version
                    try:
                        transcript_id_v = transcript_id_mapping.get(
                            transcript_id.split(".")[0], transcript_id
                        )
                    except (AttributeError, KeyError):
                        transcript_id_v = transcript_id

                    # Get sequence from FASTA
                    try:
                        row = transcripts_dict[transcript_id_v]
                        ref_seq = row.seq
                        desc = str(row.description)
                    except KeyError:
                        stats["transcript_not_found"] += 1
                        continue

                    # Determine feature types and CDS info
                    feature_types = ["exon"]
                    cds_info: list[int] = []
                    num_orfs = 3
                    if "CDS=" in desc:
                        try:
                            cds_info = [int(x) for x in desc.split(" ")[1].split("=")[1].split("-")]
                            feature_types = ["CDS", "stop_codon"]
                            num_orfs = 1
                        except (ValueError, IndexError):
                            pass

                    # Override with configured num_orfs
                    if self._num_orfs:
                        num_orfs = self._num_orfs

                    # Get genomic features from GTF
                    chrom, strand, features_info = self._get_features(
                        db, transcript_id_v, feature_types,
                    )
                    if chrom is None:
                        continue

                    # Verify chromosome match (after mapping)
                    record_chrom_refseq = chrom_mapper.map_chrom(str(record.CHROM), "refseq")
                    if chrom != record_chrom_refseq:
                        continue

                    for alt in alts:
                        # Check overlap
                        var_pos = int(record.POS)
                        if not check_overlap(var_pos, var_pos + len(ref) - 1, features_info):
                            continue

                        # Apply variant
                        coding_ref, coding_alt = get_altseq(
                            ref_seq, Seq(ref), Seq(alt), var_pos,
                            strand, features_info, cds_info,
                        )

                        if coding_alt != "":
                            ref_orfs, alt_orfs = get_orfs_vcf(
                                coding_ref, coding_alt, trans_table, num_orfs,
                            )

                            # Build header with ClinVar metadata
                            gene_symbol, _ = self.parse_geneinfo(
                                self._get_info_field(str(record.INFO), "GENEINFO"),
                            )
                            record_id = str(record.ID) if record.ID and record.ID != "." else ""
                            seq_id = "_".join(filter(None, [
                                self._protein_prefix + record_id,
                                f"{record.CHROM}.{record.POS}.{ref}.{alt}",
                                transcript_id_v,
                            ]))

                            desc_str = f"{clnsig}|{gene_symbol}" if gene_symbol else clnsig

                            write_output(
                                seq_id=seq_id,
                                desc=desc_str,
                                seqs=alt_orfs,
                                prots_fn=out_fh,
                                seqs_filter=ref_orfs,
                            )
                            stats["variants_translated"] += 1

                            if self._report_ref_seq:
                                write_output(
                                    seq_id=transcript_id_v,
                                    desc="",
                                    seqs=ref_orfs,
                                    prots_fn=out_fh,
                                )

        logger.info("ClinVar pipeline summary: %s", stats)
        return self._output_file
```

**Step 5: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_service.py -v`
Expected: All PASS

Note: The integration tests (`TestClinVarPipeline`) may need test data adjustments. The test data sequences and positions must be consistent — CDS positions in the GTF must overlap VCF variant positions, and FASTA sequences must match CDS lengths. Debug and adjust the test data files if positions don't align. The important thing is: Pathogenic/Likely_pathogenic variants get translated, Benign variants are excluded.

**Step 6: Run all tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v`
Expected: All tests PASS

**Step 7: Commit**

```bash
git add pgatk/clinvar/clinvar_service.py \
       pgatk/testdata/clinvar/mini_clinvar.vcf \
       pgatk/testdata/clinvar/mini_refseq.gtf \
       pgatk/testdata/clinvar/mini_refseq_protein.faa \
       pgatk/tests/test_clinvar/test_clinvar_service.py
git commit -m "feat: add ClinVarService pipeline for VCF-to-protein generation

Processes ClinVar VCF against RefSeq GTF using BedTools overlap detection.
Filters by CLNSIG, reuses shared get_altseq/get_orfs_vcf for variant
application and translation."
```

---

## Task 5: NCBI Data Downloader

Download NCBI RefSeq and ClinVar files from FTP.

**Files:**
- Create: `pgatk/clinvar/data_downloader.py`
- Test: `pgatk/tests/test_clinvar/test_data_downloader.py`

**Step 1: Write the failing test**

Create `pgatk/tests/test_clinvar/test_data_downloader.py`:

```python
"""Tests for NCBI data downloader (no actual downloads)."""

from pathlib import Path

import pytest

from pgatk.clinvar.data_downloader import NcbiDataDownloader


class TestNcbiDataDownloader:
    """Tests for NcbiDataDownloader."""

    def test_build_refseq_urls(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        urls = downloader.get_refseq_urls()
        assert any("GRCh38_latest_genomic.gtf.gz" in u for u in urls)
        assert any("GRCh38_latest_protein.faa.gz" in u for u in urls)
        assert any("assembly_report.txt" in u for u in urls)

    def test_build_clinvar_urls(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        urls = downloader.get_clinvar_urls()
        assert any("clinvar.vcf.gz" in u for u in urls)

    def test_expected_files_list(self):
        downloader = NcbiDataDownloader(output_dir="/tmp/test")
        files = downloader.expected_files()
        assert len(files) >= 4  # gtf, protein fasta, assembly report, clinvar vcf

    def test_output_dir_created(self, tmp_path):
        out_dir = tmp_path / "ncbi_data"
        downloader = NcbiDataDownloader(output_dir=str(out_dir))
        downloader.ensure_output_dir()
        assert out_dir.exists()
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_data_downloader.py -v`
Expected: FAIL with `ModuleNotFoundError`

**Step 3: Implement `pgatk/clinvar/data_downloader.py`**

```python
"""Download NCBI RefSeq and ClinVar reference files.

Uses urllib for FTP/HTTPS downloads. Skips existing files unless force=True.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path

from pgatk.toolbox.general import download_file, check_create_folders

logger = logging.getLogger(__name__)

_DEFAULT_REFSEQ_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/"
    "GRCh38_latest/refseq_identifiers/"
)
_DEFAULT_CLINVAR_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
)

_REFSEQ_FILES = [
    "GRCh38_latest_genomic.gtf.gz",
    "GRCh38_latest_protein.faa.gz",
    "GRCh38_latest_assembly_report.txt",
]
_CLINVAR_FILES = [
    "clinvar.vcf.gz",
    "clinvar.vcf.gz.tbi",
]


class NcbiDataDownloader:
    """Download NCBI RefSeq and ClinVar reference data."""

    def __init__(
        self,
        output_dir: str,
        refseq_base_url: str = _DEFAULT_REFSEQ_BASE,
        clinvar_base_url: str = _DEFAULT_CLINVAR_BASE,
    ) -> None:
        self._output_dir = output_dir
        self._refseq_base = refseq_base_url
        self._clinvar_base = clinvar_base_url

    def ensure_output_dir(self) -> None:
        """Create output directory if it doesn't exist."""
        check_create_folders([self._output_dir])

    def get_refseq_urls(self) -> list[str]:
        """Return list of RefSeq file URLs to download."""
        return [self._refseq_base + f for f in _REFSEQ_FILES]

    def get_clinvar_urls(self) -> list[str]:
        """Return list of ClinVar file URLs to download."""
        return [self._clinvar_base + f for f in _CLINVAR_FILES]

    def expected_files(self) -> list[str]:
        """Return list of expected local file paths after download."""
        all_files = _REFSEQ_FILES + _CLINVAR_FILES
        return [os.path.join(self._output_dir, f) for f in all_files]

    def download_all(self, force: bool = False) -> list[str]:
        """Download all required files. Returns list of downloaded file paths.

        :param force: If True, re-download even if files exist.
        """
        self.ensure_output_dir()
        downloaded = []

        all_urls = self.get_refseq_urls() + self.get_clinvar_urls()
        all_names = _REFSEQ_FILES + _CLINVAR_FILES

        for url, name in zip(all_urls, all_names):
            local_path = os.path.join(self._output_dir, name)

            if os.path.exists(local_path) and not force:
                logger.info("File already exists, skipping: %s", local_path)
                downloaded.append(local_path)
                continue

            logger.info("Downloading %s -> %s", url, local_path)
            result = download_file(url, local_path, logger)
            if result:
                downloaded.append(result)
            else:
                logger.error("Failed to download: %s", url)

        return downloaded
```

**Step 4: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_data_downloader.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add pgatk/clinvar/data_downloader.py pgatk/tests/test_clinvar/test_data_downloader.py
git commit -m "feat: add NcbiDataDownloader for RefSeq/ClinVar file downloads

Downloads GTF, protein FASTA, assembly report, and ClinVar VCF from NCBI FTP.
Skips existing files, supports --force re-download."
```

---

## Task 6: CLI Commands and CLI Registration

Wire up the two new CLI commands and register them in `cli.py`.

**Files:**
- Create: `pgatk/commands/clinvar_to_proteindb.py`
- Create: `pgatk/commands/ncbi_downloader.py`
- Modify: `pgatk/cli.py`

**Step 1: Create `pgatk/commands/clinvar_to_proteindb.py`**

```python
import click

from pgatk.clinvar.clinvar_service import ClinVarService
from pgatk.config.registry import load_config


@click.command("clinvar-to-proteindb", short_help="Generate protein database from ClinVar VCF + RefSeq GTF")
@click.option("-c", "--config_file", help="Configuration YAML file (optional)")
@click.option("-v", "--vcf", required=True, help="ClinVar VCF file path")
@click.option("-g", "--gtf", required=True, help="NCBI RefSeq GTF file path")
@click.option("-f", "--fasta", required=True, help="RefSeq protein/transcript FASTA file path")
@click.option("-a", "--assembly-report", required=True, help="NCBI assembly report file path")
@click.option("-o", "--output", required=True, help="Output protein FASTA file path")
@click.option("--clnsig-exclude", default=None,
              help="Comma-separated CLNSIG values to exclude (default: Benign,Likely_benign,Benign/Likely_benign)")
@click.option("--consequences", default=None,
              help="Comma-separated molecular consequences to include")
@click.option("-t", "--translation-table", type=int, default=None, help="Translation table (default: 1)")
@click.option("-p", "--protein-prefix", default=None, help="Prefix for variant protein IDs (default: clinvar_)")
@click.option("--report-ref-seq", is_flag=True, default=False,
              help="Also report reference protein sequences")
@click.pass_context
def clinvar_to_proteindb(ctx, config_file, vcf, gtf, fasta, assembly_report, output,
                          clnsig_exclude, consequences, translation_table, protein_prefix,
                          report_ref_seq):
    """Generate a variant protein database from ClinVar VCF and NCBI RefSeq GTF.

    This command does NOT require VEP annotations. It uses BedTools interval
    overlap to find transcripts affected by each ClinVar variant, then applies
    the variant and translates to protein.
    """
    config = load_config("clinvar", config_file)
    clinvar_cfg = config.get("clinvar_translation", {})

    # Build kwargs with config defaults, CLI overrides
    kwargs = {
        "vcf_file": vcf,
        "gtf_file": gtf,
        "fasta_file": fasta,
        "assembly_report": assembly_report,
        "output_file": output,
        "translation_table": translation_table or clinvar_cfg.get("translation_table", 1),
        "mito_translation_table": clinvar_cfg.get("mito_translation_table", 2),
        "protein_prefix": protein_prefix or clinvar_cfg.get("protein_prefix", "clinvar_"),
        "num_orfs": clinvar_cfg.get("num_orfs", 1),
        "report_ref_seq": report_ref_seq or clinvar_cfg.get("report_ref_seq", False),
    }

    if clnsig_exclude:
        kwargs["clnsig_exclude"] = [x.strip() for x in clnsig_exclude.split(",")]
    elif "clinical_significance_exclude" in clinvar_cfg:
        kwargs["clnsig_exclude"] = clinvar_cfg["clinical_significance_exclude"]

    if consequences:
        kwargs["consequences"] = [x.strip() for x in consequences.split(",")]
    elif "include_consequences" in clinvar_cfg:
        kwargs["consequences"] = clinvar_cfg["include_consequences"]

    service = ClinVarService(**kwargs)
    service.run()
```

**Step 2: Create `pgatk/commands/ncbi_downloader.py`**

```python
import click

from pgatk.clinvar.data_downloader import NcbiDataDownloader
from pgatk.config.registry import load_config


@click.command("ncbi-downloader", short_help="Download NCBI RefSeq and ClinVar reference files")
@click.option("-c", "--config_file", help="Configuration YAML file (optional)")
@click.option("-o", "--output-dir", required=True, help="Output directory for downloaded files")
@click.option("--force", is_flag=True, default=False, help="Re-download files even if they exist")
@click.pass_context
def ncbi_downloader(ctx, config_file, output_dir, force):
    """Download NCBI RefSeq GTF, protein FASTA, assembly report, and ClinVar VCF.

    Files are downloaded to the specified output directory. Existing files
    are skipped unless --force is used.
    """
    config = load_config("clinvar", config_file)
    clinvar_cfg = config.get("clinvar_translation", {})

    refseq_base = clinvar_cfg.get("ncbi_refseq_ftp",
                                   "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/")
    clinvar_base = clinvar_cfg.get("clinvar_ftp",
                                    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/")

    downloader = NcbiDataDownloader(
        output_dir=output_dir,
        refseq_base_url=refseq_base,
        clinvar_base_url=clinvar_base,
    )
    downloaded = downloader.download_all(force=force)
    click.echo(f"Downloaded {len(downloaded)} files to {output_dir}")
```

**Step 3: Register commands in `pgatk/cli.py`**

Add at line 30 (after the existing imports):

```python
from pgatk.commands import clinvar_to_proteindb as clinvar_to_proteindb_cmd
from pgatk.commands import ncbi_downloader as ncbi_downloader_cmd
```

Add at line 56 (after the last `cli.add_command`):

```python
cli.add_command(clinvar_to_proteindb_cmd.clinvar_to_proteindb)
cli.add_command(ncbi_downloader_cmd.ncbi_downloader)
```

**Step 4: Verify CLI help works**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pgatk.cli --help`
Expected: Output should include `clinvar-to-proteindb` and `ncbi-downloader` commands

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pgatk.cli clinvar-to-proteindb --help`
Expected: Shows help text with all options

**Step 5: Run all tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add pgatk/commands/clinvar_to_proteindb.py \
       pgatk/commands/ncbi_downloader.py \
       pgatk/cli.py
git commit -m "feat: add clinvar-to-proteindb and ncbi-downloader CLI commands

Wire up ClinVar pipeline and NCBI downloader as Click commands.
Register both in the main CLI entry point."
```

---

## Task 7: Final Integration Test and Cleanup

End-to-end test running the full CLI pipeline with mini test data.

**Files:**
- Create: `pgatk/tests/test_clinvar/test_clinvar_integration.py`

**Step 1: Write integration test**

Create `pgatk/tests/test_clinvar/test_clinvar_integration.py`:

```python
"""End-to-end integration test for the ClinVar pipeline."""

from pathlib import Path

from click.testing import CliRunner

from pgatk.cli import cli


TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"


class TestClinVarCLIIntegration:
    """Test the clinvar-to-proteindb CLI command end-to-end."""

    def test_cli_produces_output(self, tmp_path):
        output_file = tmp_path / "clinvar_output.fa"
        runner = CliRunner()
        result = runner.invoke(cli, [
            "clinvar-to-proteindb",
            "--vcf", str(TESTDATA_DIR / "mini_clinvar.vcf"),
            "--gtf", str(TESTDATA_DIR / "mini_refseq.gtf"),
            "--fasta", str(TESTDATA_DIR / "mini_refseq_protein.faa"),
            "--assembly-report", str(TESTDATA_DIR / "mini_assembly_report.txt"),
            "--output", str(output_file),
        ])
        assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception}"
        assert output_file.exists()

    def test_cli_help_shows_options(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["clinvar-to-proteindb", "--help"])
        assert result.exit_code == 0
        assert "clinvar" in result.output.lower()
        assert "--vcf" in result.output
        assert "--assembly-report" in result.output

    def test_ncbi_downloader_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["ncbi-downloader", "--help"])
        assert result.exit_code == 0
        assert "--output-dir" in result.output
        assert "--force" in result.output
```

**Step 2: Run integration test**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_integration.py -v`
Expected: All PASS

**Step 3: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v`
Expected: All tests PASS (existing + new)

**Step 4: Commit**

```bash
git add pgatk/tests/test_clinvar/test_clinvar_integration.py
git commit -m "test: add end-to-end CLI integration tests for ClinVar pipeline"
```

**Step 5: Clean up any gffutils .db files created during tests**

If gffutils created `.db` files in testdata during test runs, add them to `.gitignore`:

```bash
echo "*.db" >> pgatk/testdata/clinvar/.gitignore
git add pgatk/testdata/clinvar/.gitignore
git commit -m "chore: gitignore gffutils .db files in ClinVar testdata"
```
