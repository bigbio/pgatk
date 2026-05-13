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
