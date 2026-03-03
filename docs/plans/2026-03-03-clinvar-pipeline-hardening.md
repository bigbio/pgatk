# ClinVar Pipeline Hardening Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the three weaknesses in the ClinVar pipeline: activate MC consequence filtering, eliminate double VCF read, and add biotype filtering.

**Architecture:** Refactor `ClinVarService.run()` to read the VCF once into a DataFrame, build the BedTools overlap map from the DataFrame (not a second file read), apply MC and biotype filters early to skip irrelevant variants, and add a duplicate variant-transcript guard.

**Tech Stack:** Python, pandas, pybedtools, gffutils, pytest

---

### Task 1: Activate MC Consequence Filtering

The `include_consequences` config is defined in `clinvar_config.yaml` but never used.
The ClinVar VCF already has an `MC` field (e.g., `MC=SO:0001583|missense_variant`) that
tells us the variant consequence. We need to parse it and filter variants whose MC
does not match any entry in `include_consequences`.

**Files:**
- Modify: `pgatk/clinvar/clinvar_service.py:52-80` (add `_include_consequences` to `__init__`)
- Modify: `pgatk/clinvar/clinvar_service.py:340-532` (add MC filter in `run()`)
- Modify: `pgatk/config/clinvar_config.yaml:18-29` (update defaults — remove splice variants we can't model)
- Test: `pgatk/tests/test_clinvar/test_clinvar_service.py`
- Modify: `pgatk/testdata/clinvar/mini_clinvar.vcf` (add test records with different MC values)

**Step 1: Write the failing tests**

Add to `pgatk/tests/test_clinvar/test_clinvar_service.py`:

```python
class TestMCFiltering:
    """Tests for MC consequence filtering."""

    def test_passes_mc_filter_matching_consequence(self):
        """A variant with missense_variant should pass when it's in the include list."""
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001583|missense_variant"
        assert ClinVarService.passes_mc_filter(mc, include) is True

    def test_passes_mc_filter_no_match(self):
        """A variant with synonymous_variant should NOT pass."""
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001819|synonymous_variant"
        assert ClinVarService.passes_mc_filter(mc, include) is False

    def test_passes_mc_filter_multi_consequence_any_match(self):
        """If ANY consequence matches the include list, the variant passes."""
        include = ["missense_variant", "stop_gained"]
        mc = "SO:0001819|synonymous_variant,SO:0001587|stop_gained"
        assert ClinVarService.passes_mc_filter(mc, include) is True

    def test_passes_mc_filter_empty_mc_passes(self):
        """Variants without MC field should pass (no info to filter on)."""
        include = ["missense_variant"]
        assert ClinVarService.passes_mc_filter("", include) is True

    def test_passes_mc_filter_all_keyword(self):
        """When include list is ['all'], everything passes."""
        mc = "SO:0001819|synonymous_variant"
        assert ClinVarService.passes_mc_filter(mc, ["all"]) is True
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_service.py::TestMCFiltering -v`
Expected: FAIL with `AttributeError: type object 'ClinVarService' has no attribute 'passes_mc_filter'`

**Step 3: Implement `passes_mc_filter` and wire into `__init__` and `run()`**

In `pgatk/clinvar/clinvar_service.py`:

1. Add `passes_mc_filter` static method (after `passes_clnsig_filter`):

```python
@staticmethod
def passes_mc_filter(mc_field: str, include_list: list[str]) -> bool:
    """Return True when *mc_field* contains at least one consequence in *include_list*.

    An empty MC field always passes (no information to filter on).
    The special value ``'all'`` in include_list disables filtering.
    """
    if not mc_field:
        return True
    if "all" in include_list:
        return True
    consequences = ClinVarService.parse_mc_consequences(mc_field)
    return any(c in include_list for c in consequences)
```

2. In `__init__`, after `self._clnsig_exclude = ...` (line 80), add:

```python
self._include_consequences = self._cfg.get(
    "include_consequences", ["all"]
)
```

3. In `run()`, after the CLNSIG filter block (after line 401), add:

```python
# --- MC consequence filter ---
mc_field = self._get_info_field(info, "MC")
if not self.passes_mc_filter(mc_field, self._include_consequences):
    stats["variants_filtered_mc"] += 1
    continue
```

4. Add `"variants_filtered_mc": 0` to the `stats` dict (line 373).

**Step 4: Update config defaults**

In `pgatk/config/clinvar_config.yaml`, update the `include_consequences` list.
Remove splice variants (we can't model their protein consequence — that would
require exon-skipping prediction). Add a comment explaining why:

```yaml
  # Molecular consequences to include (from ClinVar MC field).
  # Only consequences that produce modelable protein changes are included.
  # Splice variants are excluded because the pipeline cannot predict
  # exon-skipping outcomes — use VEP for those.
  include_consequences:
    - "missense_variant"
    - "nonsense"
    - "frameshift_variant"
    - "inframe_insertion"
    - "inframe_deletion"
    - "stop_gained"
    - "stop_lost"
    - "start_lost"
```

**Step 5: Update test data**

Add a synonymous variant to `pgatk/testdata/clinvar/mini_clinvar.vcf` (after line 9):

```
1	1054	rs00004	G	A	.	.	GENEINFO=GeneA:1234;CLNSIG=Pathogenic;MC=SO:0001819|synonymous_variant
```

This is Pathogenic but synonymous — it should pass CLNSIG but fail MC filter.

**Step 6: Add integration test for MC filtering**

Add to `TestClinVarPipeline` in `test_clinvar_service.py`:

```python
def test_synonymous_variant_excluded(self, tmp_path):
    """Pathogenic synonymous variant (rs00004) should not appear in output."""
    output_file = str(tmp_path / "output.fa")
    service = ClinVarService(
        vcf_file=MINI_VCF,
        gtf_file=MINI_GTF,
        fasta_file=MINI_FASTA,
        assembly_report=ASSEMBLY_REPORT,
        output_file=output_file,
    )
    service.run()
    with open(output_file, "r") as f:
        content = f.read()
    assert "rs00004" not in content
```

**Step 7: Run all tests to verify**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/ -v`
Expected: ALL PASS

**Step 8: Commit**

```bash
git add pgatk/clinvar/clinvar_service.py pgatk/config/clinvar_config.yaml pgatk/tests/test_clinvar/test_clinvar_service.py pgatk/testdata/clinvar/mini_clinvar.vcf
git commit -m "feat(clinvar): activate MC consequence filtering from config"
```

---

### Task 2: Single-Pass VCF Reading

Currently the VCF is read twice: once in `_annotate_vcf_with_transcripts()` (to build
the BED for BedTools) and once in `_read_vcf()` (to get the DataFrame for processing).
Refactor so the VCF is read once into a DataFrame, and the BED is built from the
DataFrame rows.

**Files:**
- Modify: `pgatk/clinvar/clinvar_service.py:246-334` (refactor `_annotate_vcf_with_transcripts`)
- Modify: `pgatk/clinvar/clinvar_service.py:340-380` (refactor `run()` to read once)
- Test: `pgatk/tests/test_clinvar/test_clinvar_service.py`

**Step 1: Write the failing test**

Add to `test_clinvar_service.py`:

```python
class TestBuildOverlapMap:
    """Tests for _build_overlap_map (DataFrame-based BedTools annotation)."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_build_overlap_map_from_dataframe(self):
        """_build_overlap_map should accept a DataFrame and return overlap dict."""
        from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(ASSEMBLY_REPORT)
        _meta, vcf_df = ClinVarService._read_vcf(MINI_VCF)
        overlap_map = ClinVarService._build_overlap_map(vcf_df, MINI_GTF, chrom_mapper)
        assert isinstance(overlap_map, dict)
        # rs00001 at chr1:1006 should overlap NM_000001.1 CDS (1000-1299)
        assert any("1:1006:" in k for k in overlap_map)

    def test_build_overlap_map_matches_old_method(self):
        """_build_overlap_map should produce same results as _annotate_vcf_with_transcripts."""
        from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(ASSEMBLY_REPORT)
        _meta, vcf_df = ClinVarService._read_vcf(MINI_VCF)
        new_result = ClinVarService._build_overlap_map(vcf_df, MINI_GTF, chrom_mapper)
        old_result = ClinVarService._annotate_vcf_with_transcripts(MINI_VCF, MINI_GTF, chrom_mapper)
        assert new_result == old_result
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_service.py::TestBuildOverlapMap -v`
Expected: FAIL with `AttributeError: ... has no attribute '_build_overlap_map'`

**Step 3: Implement `_build_overlap_map`**

Add new method to `ClinVarService` (replace or supplement `_annotate_vcf_with_transcripts`):

```python
@staticmethod
def _build_overlap_map(
    vcf_df: pd.DataFrame,
    gtf_file: str,
    chrom_mapper: ChromosomeMapper,
) -> dict[str, list[str]]:
    """Find transcripts overlapping each VCF variant via BedTools.

    Unlike ``_annotate_vcf_with_transcripts``, this method builds the BED
    from an already-loaded DataFrame, avoiding a second file read.

    Returns a dict mapping ``"CHROM:POS:REF:ALT"`` variant keys to lists
    of overlapping transcript IDs.
    """
    bed_lines: list[str] = []
    for _, row in vcf_df.iterrows():
        ref = str(row.REF)
        if any(c not in "ACGT" for c in ref):
            continue
        chrom_numeric = str(row.CHROM)
        pos = int(row.POS)
        alt_field = str(row.ALT)
        chrom_refseq = chrom_mapper.map_chrom(chrom_numeric, "refseq")
        start = pos - 1  # BED is 0-based half-open
        end = start + len(ref)
        for alt in alt_field.split(","):
            alt = alt.strip()
            if not alt or not all(c in "ACGT" for c in alt):
                continue
            variant_key = f"{chrom_numeric}:{pos}:{ref}:{alt}"
            bed_lines.append(f"{chrom_refseq}\t{start}\t{end}\t{variant_key}\n")

    if not bed_lines:
        return {}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
        tmp.writelines(bed_lines)
        tmp_bed_path = tmp.name

    try:
        vcf_bed = BedTool(tmp_bed_path)
        gtf_bed = BedTool(gtf_file)
        intersection = vcf_bed.intersect(gtf_bed, wo=True)

        result: dict[str, list[str]] = {}
        for feature in intersection:
            fields = str(feature).strip().split("\t")
            variant_key = fields[3]
            gtf_type_idx = 4 + 2
            if len(fields) <= gtf_type_idx:
                continue
            if fields[gtf_type_idx] != "CDS":
                continue
            gtf_attrs_idx = 4 + 8
            if len(fields) <= gtf_attrs_idx:
                continue
            transcript_id = _extract_transcript_id(fields[gtf_attrs_idx])
            if transcript_id:
                result.setdefault(variant_key, [])
                if transcript_id not in result[variant_key]:
                    result[variant_key].append(transcript_id)
        return result
    finally:
        Path(tmp_bed_path).unlink(missing_ok=True)
```

**Step 4: Refactor `run()` to use single-pass reading**

In `run()`, replace the two calls:

```python
# OLD (two reads):
overlap_map = self._annotate_vcf_with_transcripts(
    self._vcf_file, self._gtf_file, chrom_mapper
)
_metadata, vcf_df = self._read_vcf(self._vcf_file)

# NEW (one read):
_metadata, vcf_df = self._read_vcf(self._vcf_file)
overlap_map = self._build_overlap_map(vcf_df, self._gtf_file, chrom_mapper)
```

**Step 5: Remove `_annotate_vcf_with_transcripts`**

Delete the old method entirely (lines 246-334). It is no longer called.

**Step 6: Run all tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/ -v`
Expected: ALL PASS

**Step 7: Commit**

```bash
git add pgatk/clinvar/clinvar_service.py pgatk/tests/test_clinvar/test_clinvar_service.py
git commit -m "refactor(clinvar): single-pass VCF reading — build overlap map from DataFrame"
```

---

### Task 3: Biotype Filtering from GTF

The config has `include_biotypes: "protein_coding"` but it's never applied. When
processing a transcript, check its `gene_biotype` attribute from the gffutils DB.
Skip transcripts that don't match the include list.

**Files:**
- Modify: `pgatk/clinvar/clinvar_service.py:52-80` (add `_include_biotypes` to `__init__`)
- Modify: `pgatk/clinvar/clinvar_service.py` in `run()` (add biotype check in transcript loop)
- Test: `pgatk/tests/test_clinvar/test_clinvar_service.py`
- Modify: `pgatk/testdata/clinvar/mini_refseq.gtf` (add `gene_biotype` attribute to existing records)

**Step 1: Update test data**

Add `gene_biotype "protein_coding";` to all GTF lines in `mini_refseq.gtf`.
For example, the first transcript line becomes:

```
NC_000001.11	BestRefSeq	transcript	1000	1299	.	+	.	gene_id "GeneA"; transcript_id "NM_000001.1"; gene_biotype "protein_coding";
```

Do this for ALL lines (transcript, exon, CDS) in the GTF.

**Step 2: Write the failing test**

Add to `test_clinvar_service.py`:

```python
class TestBiotypeFiltering:
    """Tests for biotype filtering from GTF."""

    def test_get_biotype_from_db(self):
        """_get_transcript_biotype should extract gene_biotype from gffutils DB."""
        from pgatk.clinvar.clinvar_service import ClinVarService
        db = ClinVarService._parse_gtf(MINI_GTF)
        biotype = ClinVarService._get_transcript_biotype(db, "NM_000001.1")
        assert biotype == "protein_coding"

    def test_get_biotype_missing_returns_empty(self):
        """Missing biotype returns empty string (passes filter)."""
        from pgatk.clinvar.clinvar_service import ClinVarService
        db = ClinVarService._parse_gtf(MINI_GTF)
        biotype = ClinVarService._get_transcript_biotype(db, "NONEXISTENT")
        assert biotype == ""
```

**Step 3: Run tests to verify they fail**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/test_clinvar_service.py::TestBiotypeFiltering -v`
Expected: FAIL with `AttributeError`

**Step 4: Implement biotype filtering**

1. Add `_get_transcript_biotype` static method:

```python
@staticmethod
def _get_transcript_biotype(db: gffutils.FeatureDB, transcript_id: str) -> str:
    """Extract gene_biotype from a gffutils transcript feature.

    Returns empty string if the transcript or attribute is not found.
    """
    try:
        feature = db[transcript_id]
    except gffutils.exceptions.FeatureNotFoundError:
        try:
            feature = db[transcript_id.split(".")[0]]
        except gffutils.exceptions.FeatureNotFoundError:
            return ""
    try:
        return feature.attributes.get("gene_biotype", [""])[0]
    except (IndexError, AttributeError):
        return ""
```

2. In `__init__`, load the biotype config:

```python
biotypes_raw = self._cfg.get("include_biotypes", "all")
if isinstance(biotypes_raw, str):
    self._include_biotypes = [b.strip() for b in biotypes_raw.split(",")]
else:
    self._include_biotypes = list(biotypes_raw)
```

3. In `run()`, inside the transcript loop (after looking up the transcript in
FASTA, before getting features), add:

```python
# --- Biotype filter ---
if self._include_biotypes != ["all"]:
    biotype = self._get_transcript_biotype(db, tid)
    if biotype and biotype not in self._include_biotypes:
        stats["transcripts_filtered_biotype"] += 1
        continue
```

4. Add `"transcripts_filtered_biotype": 0` to `stats`.

**Step 5: Run all tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/ -v`
Expected: ALL PASS

**Step 6: Commit**

```bash
git add pgatk/clinvar/clinvar_service.py pgatk/tests/test_clinvar/test_clinvar_service.py pgatk/testdata/clinvar/mini_refseq.gtf
git commit -m "feat(clinvar): add biotype filtering from GTF gene_biotype attribute"
```

---

### Task 4: Duplicate Variant-Transcript Guard and Stats Cleanup

The Ensembl pipeline tracks processed `(transcript_id, ref, alt)` combinations to
avoid processing the same variant-transcript pair twice. The ClinVar pipeline lacks
this. Also clean up the `_read_vcf` HEADERS dict (deferred I3 from code review).

**Files:**
- Modify: `pgatk/clinvar/clinvar_service.py` in `run()` (add duplicate guard)
- Modify: `pgatk/clinvar/clinvar_service.py:213-240` (clean up `_read_vcf`)
- Test: `pgatk/tests/test_clinvar/test_clinvar_service.py`

**Step 1: Write the failing test**

```python
class TestDuplicateGuard:
    """Pipeline should not process the same variant-transcript pair twice."""

    @pytest.fixture(autouse=True)
    def _cleanup_db(self):
        db_path = Path(MINI_GTF).with_suffix(".db")
        if db_path.exists():
            db_path.unlink()
        yield
        if db_path.exists():
            db_path.unlink()

    def test_no_duplicate_fasta_entries(self, tmp_path):
        """Each variant-transcript combo should appear at most once in output."""
        output_file = str(tmp_path / "output.fa")
        service = ClinVarService(
            vcf_file=MINI_VCF,
            gtf_file=MINI_GTF,
            fasta_file=MINI_FASTA,
            assembly_report=ASSEMBLY_REPORT,
            output_file=output_file,
        )
        service.run()
        with open(output_file, "r") as f:
            headers = [line.strip() for line in f if line.startswith(">")]
        # No duplicate headers
        assert len(headers) == len(set(headers)), (
            f"Duplicate FASTA entries found: {headers}"
        )
```

**Step 2: Implement duplicate guard in `run()`**

Before the main `for _, record in vcf_df.iterrows():` loop, add:

```python
processed_pairs: set[str] = set()
```

Inside the transcript loop, before the `# Resolve transcript in FASTA` block, add:

```python
pair_key = f"{variant_key}:{transcript_id}"
if pair_key in processed_pairs:
    continue
processed_pairs.add(pair_key)
```

**Step 3: Clean up `_read_vcf` HEADERS**

Replace the `HEADERS` dict with a simple column name list (remove the misleading
type values that are never applied):

```python
@staticmethod
def _read_vcf(vcf_file: str) -> tuple[list, pd.DataFrame]:
    """Read a VCF file and return metadata lines and a DataFrame of records."""
    COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    metadata: list[str] = []
    data: list[list[str]] = []
    with open(vcf_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                metadata.append(line)
            else:
                data.append(line.split("\t")[0:8])

    vcf_df = pd.DataFrame(data, columns=COLUMNS)
    return metadata, vcf_df
```

**Step 4: Run all tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_clinvar/ -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add pgatk/clinvar/clinvar_service.py pgatk/tests/test_clinvar/test_clinvar_service.py
git commit -m "fix(clinvar): add duplicate variant-transcript guard, clean up _read_vcf"
```

---

### Task 5: Comprehensive Integration Test Updates

Verify the full pipeline with the new filters active. Ensure test data exercises
all filter paths (CLNSIG, MC, biotype). Verify the single-pass refactor doesn't
change output.

**Files:**
- Modify: `pgatk/tests/test_clinvar/test_clinvar_integration.py`
- Modify: `pgatk/tests/test_clinvar/test_clinvar_service.py`

**Step 1: Add integration tests**

Add to `TestClinVarCLIIntegration` in `test_clinvar_integration.py`:

```python
def test_cli_filters_synonymous_variants(self, tmp_path, clinvar_testdata):
    """Synonymous variants should be filtered by MC even if Pathogenic."""
    output_file = tmp_path / "clinvar_output.fa"
    runner = CliRunner()
    result = runner.invoke(cli, [
        "clinvar-to-proteindb",
        "--vcf", str(clinvar_testdata / "mini_clinvar.vcf"),
        "--gtf", str(clinvar_testdata / "mini_refseq.gtf"),
        "--fasta", str(clinvar_testdata / "mini_refseq_protein.faa"),
        "--assembly-report", str(clinvar_testdata / "mini_assembly_report.txt"),
        "--output", str(output_file),
    ])
    assert result.exit_code == 0
    content = output_file.read_text()
    # rs00004 is Pathogenic but synonymous — should NOT appear
    assert "rs00004" not in content
    # rs00001 is Pathogenic missense — SHOULD appear
    assert "rs00001" in content
```

Add a test that existing variants still produce correct output:

```python
def test_pipeline_output_contains_clnsig_in_header(self, tmp_path, clinvar_testdata):
    """Output FASTA headers should contain CLNSIG and gene symbol."""
    output_file = tmp_path / "clinvar_output.fa"
    runner = CliRunner()
    result = runner.invoke(cli, [
        "clinvar-to-proteindb",
        "--vcf", str(clinvar_testdata / "mini_clinvar.vcf"),
        "--gtf", str(clinvar_testdata / "mini_refseq.gtf"),
        "--fasta", str(clinvar_testdata / "mini_refseq_protein.faa"),
        "--assembly-report", str(clinvar_testdata / "mini_assembly_report.txt"),
        "--output", str(output_file),
    ])
    assert result.exit_code == 0
    content = output_file.read_text()
    assert "Pathogenic" in content
```

**Step 2: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v`
Expected: ALL PASS (including existing Ensembl tests and VCF utils tests)

**Step 3: Commit**

```bash
git add pgatk/tests/test_clinvar/
git commit -m "test(clinvar): add integration tests for MC filtering and output validation"
```
