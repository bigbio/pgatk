# Test Validations

This page documents the validation strategy for PGATK, describing what the automated tests verify and how to independently confirm their correctness.
# 1. Variant Translation Tests - COSMIC
## 1.1 Validation per variant type (`test_variant_types_hgvs.py`)

These tests verify that `CancerGenomesService.get_mut_pro_seq()` correctly processes all 18 Sequence Ontology (SO) variant types present in COSMIC v103 mutation files, following [HGVS nomenclature rules](https://hgvs-nomenclature.org/stable/).

Each test uses a real mutation from `Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.tsv` paired with the matching CDS from `Cosmic_Genes_v103_GRCh38.fasta`.

---

### Actionable Variant Types

These types are expected to produce a non-empty mutant protein sequence. Tests assert biologically specific properties, not just that the result is non-empty.

| Variant Type | Gene | Mutation | HGVS Rule Verified |
|---|---|---|---|
| `missense_variant` | EZH2 | c.656C>T p.S219F | Protein length unchanged; position 218 (0-indexed) is `F`; all other residues identical to WT |
| `stop_gained` | FCGR3B | c.394A>T p.K132* | Position 131 (0-indexed) is `*` (stop codon) in the translated output |
| `inframe_deletion` | CTNNB1 | c.104_124del p.T35_G41del | Protein shortened by exactly 7 AA (21 nt deleted ÷ 3); residues before position 35 unchanged |
| `inframe_insertion` (dup) | SMAP1 | c.1226_1228dup p.I409dup | Protein lengthened by exactly 1 AA (3 nt duplicated ÷ 3); residues before the dup site unchanged |
| `protein_altering_variant` | NF2 | c.252_262delinsCT p.L85_K88delins* | Result is non-empty and differs from the wild-type sequence |
| `frameshift_variant` | MEN1 | c.292del p.R98Efs*21 | Non-empty translated sequence produced from the frame-shifted CDS |
| `start_lost` | INPP4B | c.3G>T p.M1? | First AA is not Met; it is Ile (`ATG→ATT` via the DNA substitution path) |
| `stop_lost` | MYD88 | c.611_613delinsGCGGCCCCC p.D204_*205delinsGGPR | Mutant protein is longer than WT (stop codon removed; read-through occurs) |
| `stop_retained_variant` | EPHA3 | c.2756G>A p.*919= | Length unchanged; last residue is still `*` (`TGA→TAA`, both stop codons) |

---

### Filtered / Skipped Variant Types

These types must produce an empty string `""`. There are two mechanisms:

- **Pre-filter** (`cosmic_to_proteindb`): `synonymous_variant` rows are discarded before `get_mut_pro_seq` is called.
- **`p.?` guard** (`get_mut_pro_seq`): When `aa_mut == "p.?"`, the function immediately returns `""`. All intronic, UTR, and splice variants in COSMIC v103 carry `p.?`.

| Variant Type | Mechanism | Source Mutation |
|---|---|---|
| `synonymous_variant` | Pre-filter in `cosmic_to_proteindb` | Filter logic validated directly on the SO term strings |
| `intron_variant` | `p.?` guard | RAD51C c.145+706C>T p.? |
| `3_prime_UTR_variant` | `p.?` guard | SUFU c.*2816G>A p.? |
| `5_prime_UTR_variant` | `p.?` guard | PSMD13 c.-20C>T p.? |
| `coding_sequence_variant` | `p.?` guard (also has intronic offset `+`) | OAS2 c.2157+2A>G p.? |
| `splice_acceptor_variant` | `p.?` guard | TSHR c.171-1G>A p.? |
| `splice_donor_variant` | `p.?` guard | PTCH1 c.3306+2T>G p.? |
| `splice_region_variant` (composite) | Fires missense branch when combined with `missense_variant` | EZH2 c.656C>T p.S219F with type `missense_variant,splice_region_variant` |
| `incomplete_terminal_codon_variant` | `p.?` guard (absent from v103 targeted-screen data; tested by convention) | c.2213_2214insA p.? |

!!! note
    `splice_region_variant` only appears in COSMIC v103 as part of a composite type string (e.g. `missense_variant,splice_region_variant`). When combined with an actionable type, the actionable type's branch fires first and the mutation is processed normally.

---

## 1.2 Independent Validation of the Tests

### 1. Re-derive CDS sequences from the FASTA

The CDS strings embedded in the test file came from `Cosmic_Genes_v103_GRCh38.fasta`. To verify a gene (e.g. EZH2):

```bash
python3 - <<'EOF'
from Bio import SeqIO
fasta = "use-cases/cosmic_data/Cosmic_Genes_v103_GRCh38.fasta"
for rec in SeqIO.parse(fasta, "fasta"):
    if "EZH2" in rec.id and "ENST00000483967" in rec.id:
        print(f"len={len(rec.seq)}")
        print(f"prot[218]={rec.seq.translate(to_stop=False)[218]}")
        print(rec.seq[:50])
EOF
```

Expected: `len=2211`, `prot[218]=S` (which becomes `F` after the c.656C>T substitution).

### 2. Confirm mutations exist in the COSMIC TSV

Every test documents its source mutation. To verify a mutation is real and not fabricated:

```bash
grep "EZH2" use-cases/cosmic_data/Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.tsv | \
  grep "c.656C>T" | \
  cut -f1-10
```

### 3. Manually compute expected protein output

For the missense test, manually apply the substitution and translate:

```python
from Bio.Seq import Seq

cds = '''ATGGGCCAGACTGGGAAGAAATCTGAGAAGGGACCAGTTTGTTGGCGGAAGCGTGTAAAATCAGAGTACATGCGACTGAGACAGCTCAAGAGGTTCAGACGAGCTGATGAAGTAAAGAGTATGTTTAGTTCCAATCGTCAGAAAATT
TTGGAAAGAACGGAAATCTTAAACCAAGAATGGAAACAGCGAAGGATACAGCCTGTGCACATCCTGACTTCTTGTTCGGTGACCAGTGACTTGGATTTTCCAACACAAGTCATCCCATTAAAGACTCTGAATGCAGTTGCTTCAGTACCCATAATG
TATTCTTGGTCTCCCCTACAGCAGAATTTTATGGTGGAAGATGAAACTGTTTTACATAACATTCCTTATATGGGAGATGAAGTTTTAGATCAGGATGGTACTTTCATTGAAGAACTAATAAAAAATTATGATGGGAAAGTACACGGGGATAGAGAA
TGTGGGTTTATAAATGATGAAATTTTTGTGGAGTTGGTGAATGCCCTTGGTCAATATAATGATGATGACGATGATGATGATGGAGACGATCCTGAAGAAAGAGAAGAAAAGCAGAAAGATCTGGAGGATCACCGAGATGATAAAGAAAGCCGCCCA
CCTCGGAAATTTCCTTCTGATAAAATTTTTGAAGCCATTTCCTCAATGTTTCCAGATAAGGGCACAGCAGAAGAACTAAAGGAAAAATATAAAGAACTCACCGAACAGCAGCTCCCAGGCGCACTTCCTCCTGAATGTACCCCCAACATAGATGGA
CCAAATGCTAAATCTGTTCAGAGAGAGCAAAGCTTACACTCCTTTCATACGCTTTTCTGTAGGCGATGTTTTAAATATGACTGCTTCCTACATCCTTTTCATGCAACACCCAACACTTATAAGCGGAAGAACACAGAAACAGCTCTAGACAACAAA
CCTTGTGGACCACAGTGTTACCAGCATTTGGAGGGAGCAAAGGAGTTTGCTGCTGCTCTCACCGCTGAGCGGATAAAGACCCCACCAAAACGTCCAGGAGGCCGCAGAAGAGGACGGCTTCCCAATAACAGTAGCAGGCCCAGCACCCCCACCATT
AATGTGCTGGAATCAAAGGATACAGACAGTGATAGGGAAGCAGGGACTGAAACGGGGGGAGAGAACAATGATAAAGAAGAAGAAGAGAAGAAAGATGAAACTTCGAGCTCCTCTGAAGCAAATTCTCGGTGTCAAACACCAATAAAGATGAAGCCA
AATATTGAACCTCCTGAGAATGTGGAGTGGAGTGGTGCTGAAGCCTCAATGTTTAGAGTCCTCATTGGCACTTACTATGACAATTTCTGTGCCATTGCTAGGTTAATTGGGACCAAAACATGTAGACAGGTGTATGAGTTTAGAGTCAAAGAATCT
AGCATCATAGCTCCAGCTCCCGCTGAGGATGTGGATACTCCTCCAAGGAAAAAGAAGAGGAAACACCGGTTGTGGGCTGCACACTGCAGAAAGATACAGCTGAAAAAGGACGGCTCCTCTAACCATGTTTACAACTATCAACCCTGTGATCATCCA
CGGCAGCCTTGTGACAGTTCGTGCCCTTGTGTGATAGCACAAAATTTTTGTGAAAAGTTTTGTCAATGTAGTTCAGAGTGTCAAAACCGCTTTCCGGGATGCCGCTGCAAAGCACAGTGCAACACCAAGCAGTGCCCGTGCTACCTGGCTGTCCGA
GAGTGTGACCCTGACCTCTGTCTTACTTGTGGAGCCGCTGACCATTGGGACAGTAAAAATGTGTCCTGCAAGAACTGCAGTATTCAGCGGGGCTCCAAAAAGCATCTATTGCTGGCACCATCTGACGTGGCAGGCTGGGGGATTTTTATCAAAGAT
CCTGTGCAGAAAAATGAATTCATCTCAGAATACTGTGGAGAGATTATTTCTCAAGATGAAGCTGACAGAAGAGGGAAAGTGTATGATAAATACATGTGCAGCTTTCTGTTCAACTTGAACAATGATTTTGTGGTGGATGCAACCCGCAAGGGTAAC
AAAATTCGTTTTGCAAATCATTCGGTAAATCCAAACTGCTATGCAAAAGTTATGATGGTTAACGGTGATCACAGGATAGGTATTTTTGCCAAGAGAGCCATCCAGACTGGCGAAGAGCTGTTTTTTGATTACAGATACAGCCAGGCTGATGCCCTG
AAGTATGTCGGCATCGAAAGAGAAATGGAAATCCCTTGA'''.replace('\n','')  # EZH2 CDS string

seq = Seq(cds)
wt_prot = seq.translate(to_stop=False)
print(wt_prot[218])   # 'S' — wild-type Serine at position 219

# Apply c.656C>T: position 656 is 0-indexed 655
mut_seq = seq[:655] + "T" + seq[656:]
mut_prot = mut_seq.translate(to_stop=False)
print(mut_prot[218])  # 'F' — mutant Phenylalanine
```

---

## 1.3 Core Unit Tests (`test_cosmic_core.py`)

14 Unit tests are implemented to verify the internal logic of `get_mut_pro_seq()` using minimal synthetic sequences (`Seq("ATGAAATTT")` etc.).

| Class | Coverage |
|---|---|
| `TestGetMutProSeqDnaBranch` | DNA substitution, insertion, single-base deletion, two-position deletion, ambiguous `c.?` routing to protein branch |
| `TestGetMutProSeqProteinBranch` | Missense, nonsense truncation, in-frame insertion, two-position in-frame deletion, single-position in-frame deletion |
| `TestGetMutProSeqEdgeCases` | Mutation beyond sequence length, ambiguous `p.?` returns `""`, non-alpha last char in aa_mut returns `""`, DNA ref-allele mismatch returns `""` |

These tests are fast to run and isolate individual code paths without requiring reference files.

---

# 2. COSMIC Integration Pipeline Tests

## 2.1 COSMIC v103 Schema and Data Model

COSMIC v103 introduced two breaking changes from earlier releases:

**Column name changes** — all columns are now `UPPERCASE_UNDERSCORE`:

| Old (v2/v98) | New (v103) |
|---|---|
| `Gene name` | `GENE_SYMBOL` |
| `Accession Number` | `TRANSCRIPT_ACCESSION` |
| `CDS mutation` | `MUTATION_CDS` |
| `AA Mutation` | `MUTATION_AA` |
| `Type` | `MUTATION_DESCRIPTION` |
| `Primary site` | *(moved to Classification file)* |

**Tissue-type data architecture** — `PRIMARY_SITE` (tissue of origin) is no longer a column in the mutation export. It lives in a separate **Classification file** (`Cosmic_Classification_Tsv_v103_GRCh38.tar`) and is joined via `COSMIC_PHENOTYPE_ID`:

```
mutation_file.COSMIC_PHENOTYPE_ID (e.g. COSO1001)
    → Cosmic_Classification_v103_GRCh38.tsv.COSMIC_PHENOTYPE_ID
    → PRIMARY_SITE (e.g. upper_aerodigestive_tract)
```

`COSMIC_PHENOTYPE_ID` is present at column 6 in both `Cosmic_CompleteTargetedScreensMutant` and `Cosmic_GenomeScreensMutant` files. No hop through the Sample file is needed.

> **Note**: The Sample file (`Cosmic_Sample_Tsv_v103_GRCh38.tsv`) contains 34 columns covering sequencing metadata (`SAMPLE_NAME`, `COSMIC_SAMPLE_ID`, `TUMOUR_ID`, etc.) but does **not** contain `PRIMARY_SITE`. The Classification file is the authoritative source for tissue and histology.

---

## 2.2 Integration Test (`test_cosmic_to_proteindb`)

**Test file**: `pgatk/tests/pgatk_tests.py::PgatkRunnerTests::test_cosmic_to_proteindb`

The test runs the full `cosmic-to-proteindb` pipeline end-to-end using small testdata files that mirror the real v103 schema:

| File | Purpose |
|---|---|
| `testdata/test_cosmic_mutations.tsv` | 12 HRAS mutation rows in v103 27-column format; covers `missense_variant`, `stop_gained`, `inframe_deletion`, `synonymous_variant` (filtered), `frameshift_variant` (`p.?`, skipped) |
| `testdata/test_cosmic_genes.fa` | HRAS CDS FASTA from `Cosmic_Genes_v103_GRCh38.fasta` |
| `testdata/test_cosmic_classification.tsv` | 8 rows mapping `COSMIC_PHENOTYPE_ID` (COSO1001–COSO1008) to `PRIMARY_SITE`; covers `upper_aerodigestive_tract`, `skin`, `bone`, `thyroid`, `liver`, `haematopoietic_and_lymphoid_tissue` |

CLI invocation used in the test:

```
cosmic-to-proteindb
  --config_file config/cosmic_config.yaml
  --input_mutation testdata/test_cosmic_mutations.tsv
  --input_genes testdata/test_cosmic_genes.fa
  --output_db testdata/test_cosmic_mutations_proteindb.fa
  --clinical_sample_file testdata/test_cosmic_classification.tsv
  --filter_column PRIMARY_SITE
  --split_by_filter_column
  --accepted_values all
```

Expected output files (split per `PRIMARY_SITE`):

```
testdata/test_cosmic_mutations_proteindb_bone.fa
testdata/test_cosmic_mutations_proteindb_liver.fa
testdata/test_cosmic_mutations_proteindb_skin.fa
testdata/test_cosmic_mutations_proteindb_thyroid.fa
testdata/test_cosmic_mutations_proteindb_upperaerodigestivetract.fa
```

The `synonymous_variant` row (COSO1006) is pre-filtered and produces no output. The `frameshift_variant` row carrying `p.?` (COSO1008) is skipped by the `p.?` guard. The consolidated `testdata/test_cosmic_mutations_proteindb.fa` contains all non-filtered mutations regardless of tissue type.

---

## 2.3 `--clinical_sample_file` CLI Option

The `--clinical_sample_file` / `-cl` option accepts any TSV file that maps a sample/phenotype ID column to the `--filter_column` value. For COSMIC, this is the Classification file:

| Argument | COSMIC v103 value |
|---|---|
| `--clinical_sample_file` | `Cosmic_Classification_v103_GRCh38.tsv` (or `.gz`) |
| `--filter_column` | `PRIMARY_SITE` |
| Join key used internally | `COSMIC_PHENOTYPE_ID` (auto-detected for COSMIC) |

For **cBioportal**, the equivalent is the clinical sample file with `SAMPLE_ID` as the join key — the `--clinical_sample_file` option uses `SAMPLE_ID` by default and switches to `COSMIC_PHENOTYPE_ID` only when invoked via `cosmic-to-proteindb`.

---

## 2.4 Independent Validation with Real COSMIC Files

After downloading COSMIC v103 files into `use-cases/use-case2/cosmic_data/`, run the pipeline against a subset of real mutations:

```bash
# Extract 500 actionable mutations (skip intronic/UTR/synonymous)
zcat use-cases/use-case2/cosmic_data/Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.tsv.gz \
  | head -1 > /tmp/cosmic_test.tsv
zcat use-cases/use-case2/cosmic_data/Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.tsv.gz \
  | grep -E "missense_variant|stop_gained|inframe_deletion|inframe_insertion|frameshift" \
  | head -500 >> /tmp/cosmic_test.tsv

pgatk cosmic-to-proteindb \
  --config_file pgatk/config/cosmic_config.yaml \
  --input_mutation /tmp/cosmic_test.tsv \
  --input_genes <(zcat use-cases/use-case2/cosmic_data/Cosmic_Genes_v103_GRCh38.fasta.gz) \
  --output_db /tmp/cosmic_out.fa \
  --clinical_sample_file <(zcat use-cases/use-case2/cosmic_data/Cosmic_Classification_v103_GRCh38.tsv.gz) \
  --filter_column PRIMARY_SITE \
  --split_by_filter_column \
  --accepted_values all
```

Expected: ~27 per-tissue output files, including canonical COSMIC tissue types (`lung`, `breast`, `large_intestine`, `skin`, `haematopoietic_and_lymphoid_tissue`, etc.). FASTA headers follow the format:

```
>COSMIC:GENE_SYMBOL:p.MutAA:SO_term
```

To verify the join is working (i.e. phenotype IDs are resolving to tissue names):

```bash
# Count sequences per tissue file
for f in /tmp/cosmic_out_*.fa; do
  echo "$f: $(grep -c '^>' $f) sequences"
done
```

A file named `cosmic_out_NS.fa` is expected for phenotypes where tissue is not specified in the Classification file — these are retained rather than silently dropped.

---

# 3. Variant Translation Tests - Ensembl vcf-to-proteindb

## 3.1 Testdata Files

Small, self-contained testdata covering 10 representative variants across all consequence types. All variants are real (sourced from Ensembl GRCh38.115 release 115, chr22) and involve two genes:

- **OR11H1** (`ENST00000643195`, mRNA, single-exon olfactory receptor gene, CDS 948 nt / 315 AA + stop) — used for all coding consequence types
- **BID** (`ENST00000550946`, ncRNA/retained_intron on chr22) — used for the non-coding transcript example

| File | Role |
|---|---|
| `pgatk/testdata/test_ensembl_v2p.vcf` | 10 variants covering 9 consequence types (Ensembl v115 VCF format) |
| `pgatk/testdata/test_ensembl_v2p.fa` | Transcript sequences for ENST00000643195 and ENST00000550946 |
| `pgatk/testdata/test_ensembl_v2p.gtf` | Gene models for both transcripts (GRCh38.p14) |
| `pgatk/testdata/test_ensembl_v2p_proteindb.fa` | Expected output — 10 sequences (7 mRNA + 3 ncRNA frames) |

The wild-type OR11H1 protein (ENST00000643195) used as comparison baseline:

```
MNVSEPNSSFAFVNEFILQGFSCEWTIQIFLFSLFTTTYALTITGNGAIAFVLW...KVLGSSNII*   (316 chars incl. stop)
```

---

## 3.2 CSQ Annotation Field Format

The VCF `##INFO` header declares the consequence annotation format:

```
##INFO=<ID=CSQ,...,Description="...Format=Allele|Consequence|Feature_type|Feature|Amino_acids|SIFT">
```

Each variant can carry multiple comma-separated CSQ entries (one per overlapping transcript). The pipeline selects entries where `Feature_type` matches `--include_biotypes` and `Consequence` is not in `--exclude_consequences`.

---

## 3.3 Filtering Logic

| Parameter | Default value | Effect |
|---|---|---|
| `--include_biotypes` | `protein_coding,...` | Only transcripts whose `Feature_type` matches are processed |
| `--exclude_consequences` | `downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant, regulatory_region_variant` | Variants with these consequences produce no output |
| `--af_field` | *(empty)* | If set, variants with MAF below `--af_threshold` (default 0.01) are skipped |
| `--num_orfs` | `3` | For ncRNA transcripts, all 3 reading frames are translated and emitted separately |

---

## 3.4 Validated Variant Examples

All mRNA examples use `ENST00000643195` (OR11H1 gene, chr22). The ncRNA example uses `ENST00000550946` (BID gene, chr22).

### Actionable Types (produce output sequences)

| Consequence | rsID | Position | REF→ALT | VCF `Amino_acids` | Expected assertion |
|---|---|---|---|---|---|
| `start_lost` | rs1211697244 | 22:15528192 | A→G | `M/V` | First AA is `V` (not `M`); length unchanged |
| `missense_variant` | rs1410655344 | 22:15528195 | A→G | `N/D` | Position 1 (0-based) is `D`; length unchanged (316 chars) |
| `frameshift_variant` | rs1402769459 | 22:15528197 | TG→T | `V/X` | Frame shifts after position 1; entire downstream sequence differs from WT |
| `stop_gained` | rs1420478920 | 22:15528234 | G→T | `E/*` | Position 14 (0-based) is `*`; protein truncated there |
| `inframe_deletion` | rs1203023715 | 22:15528914 | CTTC→C | `AFS/AS` | `AFS` → `AS`; length is 315 chars (1 AA shorter than WT) |
| `protein_altering_variant` | rs1986039639 | 22:15528961 | G→GCTG | `SS/SCS` | `SS` → `SCS`; length is 317 chars (1 AA longer than WT) |
| `stop_lost` | rs1986046473 | 22:15529137 | T→A | `*/K` | Ends with `K` (not `*`); length is 316 chars (stop replaced by K) |
| `non_coding_transcript_exon_variant` | rs147461488 | 22:17740102 | GGCCACGCTCAACT→G | — | 3 output sequences (`_1`, `_2`, `_3`) for all reading frames |

### No-Output Types

| Consequence | rsID | Position | REF→ALT | VCF `Amino_acids` | Reason no output sequence |
|---|---|---|---|---|---|
| `synonymous_variant` | rs1394965478 | 22:15528197 | T→C | `N/N` | In `--exclude_consequences` by default |
| `stop_retained_variant` | rs1986046579 | 22:15529138 | A→G | `*/*` | Stop codon TAA→TGA: protein unchanged, no distinct variant sequence written |

---

## 3.5 Output FASTA Header Format

```
>var_<rsid>_<chrom>.<pos>.<REF>.<ALT>_<transcript_id>[_<frame>]
```

- **mRNA transcripts** — no frame suffix; the annotated CDS frame is used:
  ```
  >var_rs1410655344_22.15528195.A.G_ENST00000643195
  ```
- **ncRNA transcripts** — three entries per variant with suffixes `_1`, `_2`, `_3` (one per reading frame, controlled by `--num_orfs 3`):
  ```
  >var_rs147461488_22.17740102.GGCCACGCTCAACT.G_ENST00000550946_1
  >var_rs147461488_22.17740102.GGCCACGCTCAACT.G_ENST00000550946_2
  >var_rs147461488_22.17740102.GGCCACGCTCAACT.G_ENST00000550946_3
  ```

---

## 3.6 Independent Validation

### 1. Run the pipeline on testdata

```bash
cd pgatk/  # package root

pgatk vcf-to-proteindb \
  --config_file config/ensembl_config.yaml \
  --vcf testdata/test_ensembl_v2p.vcf \
  --input_fasta testdata/test_ensembl_v2p.fa \
  --gene_annotations_gtf testdata/test_ensembl_v2p.gtf \
  --output_proteindb /tmp/v2p_out.fa \
  --protein_prefix var \
  --annotation_field_name CSQ \
  --biotype_str feature_type \
  --include_biotypes mRNA,ncRNA
```

Expected `Translation summary` log:
```
# variants with invalid record:0
# variants not passing Filter:0
# variants not passing AF threshold:0
# feature IDs from VCF that are not found in the given FASTA file:0
# variants successfully translated:9
```

Expected output: **10 sequences** (7 mRNA + 3 ncRNA frames). The `synonymous_variant` and `stop_retained_variant` produce no output sequence.

### 2. Verify headers in output

```bash
OUTPUT=/tmp/v2p_out.fa

# All 8 actionable mRNA variants should appear
for rsid in rs1211697244 rs1410655344 rs1402769459 rs1420478920 rs1203023715 rs1986039639 rs1986046473; do
  count=$(grep -c "^>var_${rsid}" $OUTPUT)
  echo "$rsid: $count  (expected 1)"
done

# ncRNA variant: 3 frame sequences
grep -c "^>var_rs147461488" $OUTPUT
# Expected: 3

# Synonymous must NOT appear
grep -c "^>var_rs1394965478" $OUTPUT
# Expected: 0

# Stop-retained must NOT appear (stop codon unchanged → protein identical to reference)
grep -c "^>var_rs1986046579" $OUTPUT
# Expected: 0
```

### 3. Confirm amino acid changes

```python
from Bio import SeqIO

records = SeqIO.to_dict(SeqIO.parse("/tmp/v2p_out.fa", "fasta"))

# start_lost (M/V): position 0 is V, not M
sl = str(records["var_rs1211697244_22.15528192.A.G_ENST00000643195"].seq)
assert sl[0] == "V", f"Expected V at pos 0, got {sl[0]}"

# missense (N/D): position 1 is D; length = 316 (WT reference)
mut = str(records["var_rs1410655344_22.15528195.A.G_ENST00000643195"].seq)
assert mut[1] == "D", f"Expected D at pos 1, got {mut[1]}"
wt_len = len(mut)  # 316 — used as reference length below

# stop_gained (E/*): stop at position 14
trunc = str(records["var_rs1420478920_22.15528234.G.T_ENST00000643195"].seq)
assert trunc[14] == "*", f"Expected * at pos 14, got {trunc[14]}"

# inframe_deletion (AFS/AS): 1 AA shorter
del_seq = str(records["var_rs1203023715_22.15528914.CTTC.C_ENST00000643195"].seq)
assert len(del_seq) == wt_len - 1, f"Expected {wt_len-1}, got {len(del_seq)}"

# protein_altering (SS/SCS): 1 AA longer
ins_seq = str(records["var_rs1986039639_22.15528961.G.GCTG_ENST00000643195"].seq)
assert len(ins_seq) == wt_len + 1, f"Expected {wt_len+1}, got {len(ins_seq)}"

# stop_lost (*/K): ends with K, not *; same length as WT (stop replaced by K)
ext = str(records["var_rs1986046473_22.15529137.T.A_ENST00000643195"].seq)
assert ext[-1] == "K", f"Expected K at end, got {ext[-1]}"
assert len(ext) == wt_len, f"Expected {wt_len}, got {len(ext)}"

# ncRNA: 3 frames present
for frame in ["1", "2", "3"]:
    key = f"var_rs147461488_22.17740102.GGCCACGCTCAACT.G_ENST00000550946_{frame}"
    assert key in records, f"Missing ncRNA frame {frame}"

print("All assertions passed.")
```

### 4. Compare output against committed reference

```bash
# The committed expected output can be diff'd against a fresh pipeline run
diff <(grep "^>" testdata/test_ensembl_v2p_proteindb.fa | sort) \
     <(grep "^>" /tmp/v2p_out.fa | sort)
# Expected: no output (headers match exactly)

diff testdata/test_ensembl_v2p_proteindb.fa /tmp/v2p_out.fa
# Expected: no output (files identical)
```

### 5. Extended validation with AF filtering

To test the `--af_field MAF` filter (the testdata VCF includes a `MAF` field):

```bash
# With AF threshold: only variants with MAF ≥ 0.01 are translated
# rs1410655344 MAF=0.012 → passes; rs1394965478 MAF=0.008 → filtered (also excluded by consequence)
pgatk vcf-to-proteindb \
  --config_file config/ensembl_config.yaml \
  --vcf testdata/test_ensembl_v2p.vcf \
  --input_fasta testdata/test_ensembl_v2p.fa \
  --gene_annotations_gtf testdata/test_ensembl_v2p.gtf \
  --output_proteindb /tmp/v2p_af.fa \
  --protein_prefix var \
  --annotation_field_name CSQ \
  --biotype_str feature_type \
  --include_biotypes mRNA,ncRNA \
  --af_field MAF

grep -c "^>" /tmp/v2p_af.fa
# Expected: 4  (missense MAF=0.012, ncRNA MAF=0.168 → 1 mRNA + 3 ncRNA frames)
```
