"""Translation-correctness tests for gnomAD vcf-to-proteindb.

Unit tests (TestGetAltseqGnomADOR11H1) call get_altseq() + get_orfs_vcf() directly
on OR11H1 (ENST00000643195, chr22 + strand, single-exon) using 3-frame exon translation
regardless of the FASTA header.

Integration tests (TestGnomADPipelineTranslation) run _vcf_to_proteindb_chunk
end-to-end with the gffread-generated FASTA (CDS=1-948 header), which triggers
1-frame CDS translation: one record per variant, no _1/_2/_3 suffixes in record IDs.

OR11H1 genomic layout (GRCh38, chr22, + strand):
  Exon / CDS: 15528192–15529136  (945 bp = 315 coding codons)
  Stop codon: 15529137–15529139  (3 bp)
  Total FASTA: 948 bp; CDS=1-948 in gffread header
  No 5′ UTR — FASTA sequence begins at ATG start codon.

GENCODE FASTA key (gffread format, space-delimited, CDS= header):
  SeqIO key: ENST00000643195.1
  Header: >ENST00000643195.1 CDS=1-948 gene_type=protein_coding;gene_name=OR11H1;...

Variants in test_gnomad_bgz.vcf.gz (real gnomAD v4.1.1 records, 44-field VEP):
  rs1410655344  22:15528195 A>G   missense N2D     AF_afr=6.0507e-05  → FILTERED (AF < 0.01)
  rs751806421   22:15528246 C>T   stop_gained Q19* AF_afr=0.02        → FILTERED (consequence)
  rs199856986   22:15528424 T>C   missense V78A    AF_afr=0.1239      → PASS
  rs1203023715  22:15528913 CCTT>C inframe_del     AF_afr=0.015       → PASS
"""
from __future__ import annotations

import warnings
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from pgatk.toolbox.vcf_utils import get_altseq, get_orfs_vcf

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Absolute testdata directory — needed for module-scoped fixtures that run
# before the autouse _chdir_to_package_root function-scoped fixture fires.
_TESTDATA = Path(__file__).resolve().parent.parent / "testdata"

# Exon features for OR11H1 ENST00000643195 (1-based, inclusive, + strand)
OR11H1_EXON_FEAT = [[15528192, 15529139, "exon"]]

# gffread-format FASTA key — SeqIO uses only the text before the first space
_GENCODE_OR11H1_KEY = "ENST00000643195.1"


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def or11h1_gnomad_seq():
    """OR11H1 sequence from test_gnomad_gencode.fa (gffread, CDS=1-948 header)."""
    idx = SeqIO.index(str(_TESTDATA / "test_gnomad_gencode.fa"), "fasta")
    return idx[_GENCODE_OR11H1_KEY].seq


def _altseq_protein(ref_seq, ref_allele, alt_allele, pos, features=None):
    """Return (ref_protein_frame1, alt_protein_frame1) for a + strand variant.

    Always uses 3-frame exon translation (explicit exon features, no cds_info);
    returns the frame-1 ORF (index 0).  Since OR11H1 has no 5′ UTR, frame 1
    equals the canonical CDS protein in all cases.
    """
    features = features if features is not None else OR11H1_EXON_FEAT
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coding_ref, coding_alt = get_altseq(
            ref_seq, Seq(ref_allele), Seq(alt_allele), pos, "+", features, None
        )
        ref_orfs, alt_orfs = get_orfs_vcf(coding_ref, coding_alt, 1, 3)
    return str(ref_orfs[0]), str(alt_orfs[0])


# ===========================================================================
# TestGetAltseqGnomADOR11H1 — unit tests per mutation type
# ===========================================================================

class TestGetAltseqGnomADOR11H1:
    """Direct get_altseq / get_orfs_vcf tests on OR11H1 using 3-frame exon mode.

    3-frame exon translation is used (cds_info=None, features=exon).
    Since OR11H1 has no 5′ UTR, frame 1 equals the canonical CDS protein.
    """

    def test_reference_protein_starts_with_met(self, or11h1_gnomad_seq):
        """Frame-1 translation begins with Met (ATG at position 0 of FASTA)."""
        ref_p, _ = _altseq_protein(or11h1_gnomad_seq, "A", "A", 15528192)
        assert ref_p[0] == "M", f"expected M, got {ref_p[0]}"

    def test_reference_protein_length(self, or11h1_gnomad_seq):
        """FASTA is 948 nt → 948/3 = 316 codons (315 coding + 1 stop)."""
        ref_p, _ = _altseq_protein(or11h1_gnomad_seq, "A", "A", 15528192)
        assert len(ref_p) == 316

    def test_reference_protein_ends_with_stop(self, or11h1_gnomad_seq):
        """Canonical frame-1 translation ends with '*'."""
        ref_p, _ = _altseq_protein(or11h1_gnomad_seq, "A", "A", 15528192)
        assert ref_p[-1] == "*"

    def test_three_frame_translation_returns_three_orfs(self, or11h1_gnomad_seq):
        """get_orfs_vcf with num_orfs=3 returns exactly 3 ORF sequences."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coding_ref, coding_alt = get_altseq(
                or11h1_gnomad_seq, Seq("A"), Seq("A"), 15528192,
                "+", OR11H1_EXON_FEAT, None,
            )
            ref_orfs, alt_orfs = get_orfs_vcf(coding_ref, coding_alt, 1, 3)
        assert len(ref_orfs) == 3, f"expected 3 ORFs, got {len(ref_orfs)}"
        assert len(alt_orfs) == 3

    def test_frame1_is_canonical_protein(self, or11h1_gnomad_seq):
        """Frame-1 ORF (index 0) gives canonical MNVSE... protein."""
        ref_p, _ = _altseq_protein(or11h1_gnomad_seq, "A", "A", 15528192)
        assert ref_p[:5] == "MNVSE", f"expected MNVSE, got {ref_p[:5]}"

    # --- missense (rs199856986, 15528424 T>C, Val78Ala) -----------------------

    def test_missense_changes_single_aa(self, or11h1_gnomad_seq):
        """22:15528424 T>C (rs199856986): codon 78 gTc→gCc (Val→Ala)."""
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "T", "C", 15528424)
        assert ref_p[77] == "V", f"expected V at AA 77, got {ref_p[77]}"
        assert alt_p[77] == "A", f"expected A at AA 77, got {alt_p[77]}"

    def test_missense_length_unchanged(self, or11h1_gnomad_seq):
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "T", "C", 15528424)
        assert len(ref_p) == len(alt_p)

    def test_missense_only_one_aa_differs(self, or11h1_gnomad_seq):
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "T", "C", 15528424)
        diffs = sum(1 for r, a in zip(ref_p, alt_p) if r != a)
        assert diffs == 1

    # --- inframe_deletion (rs1203023715, 15528913 CCTT>C, Phe242del) ----------

    def test_inframe_deletion_shortens_protein_by_one_aa(self, or11h1_gnomad_seq):
        """22:15528913 CCTT>C (rs1203023715): 3-bp deletion removes 1 codon (Phe242)."""
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "CCTT", "C", 15528913)
        assert len(alt_p) == len(ref_p) - 1

    def test_inframe_deletion_preserves_upstream_aas(self, or11h1_gnomad_seq):
        """Residues before the deletion site (codons 1–240) are identical in ref and alt."""
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "CCTT", "C", 15528913)
        del_aa = (15528913 - 15528192) // 3  # 721 // 3 = 240
        assert ref_p[:del_aa] == alt_p[:del_aa]

    # --- AF-filtered missense (rs1410655344, 15528195 A>G, Asn2Asp) -----------

    def test_af_filtered_missense_is_real_missense(self, or11h1_gnomad_seq):
        """22:15528195 A>G (rs1410655344): real missense N→D at AA 1 (AF_afr=6e-05)."""
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "A", "G", 15528195)
        assert ref_p[1] == "N", f"expected N at AA 1, got {ref_p[1]}"
        assert alt_p[1] == "D", f"expected D at AA 1, got {alt_p[1]}"

    # --- stop_gained (rs751806421, 15528246 C>T) — filtered by consequence ---

    def test_stop_gained_inserts_stop_at_aa18(self, or11h1_gnomad_seq):
        """22:15528246 C>T (rs751806421): codon 19 Caa→Taa (Gln→*) at index 18."""
        ref_p, alt_p = _altseq_protein(or11h1_gnomad_seq, "C", "T", 15528246)
        assert ref_p[18] == "Q", f"expected Q at AA 18, got {ref_p[18]}"
        assert alt_p[18] == "*", f"expected * at AA 18, got {alt_p[18]}"


# ===========================================================================
# TestGnomADPipelineTranslation — full pipeline integration
# ===========================================================================

_GNOMAD_CHUNK_ARGS = {
    "annotation_field_name": "vep",
    "transcript_str": "Feature",
    "consequence_str": "Consequence",
    "biotype_str": "BIOTYPE",
    "include_biotypes": "protein_coding",
    "include_consequences": "missense_variant,inframe_deletion",
    "af_field": "AF_afr",
    "af_threshold": 0.01,
    "protein_prefix": "gnomadvar",
}


@pytest.fixture(scope="module")
def gnomad_pipeline_records(tmp_path_factory):
    """Run _vcf_to_proteindb_chunk on gnomAD testdata; return parsed records."""
    from pgatk.ensembl.ensembl import EnsemblDataService
    out = str(tmp_path_factory.mktemp("gnomad_translation") / "out.fa")
    svc = EnsemblDataService({}, _GNOMAD_CHUNK_ARGS)
    svc._vcf_to_proteindb_chunk(
        str(_TESTDATA / "test_gnomad_bgz.vcf.gz"),
        str(_TESTDATA / "test_gnomad_gencode.fa"),
        str(_TESTDATA / "test_gnomad_gencode.gtf"),
        out,
    )
    records = {}
    if Path(out).exists():
        for rec in SeqIO.parse(out, "fasta"):
            records[rec.id] = str(rec.seq)
    return records


@pytest.fixture(scope="module")
def gnomad_ref_protein():
    """Reference OR11H1 protein (frame 1, 3-frame exon mode) for diff comparison.

    Since OR11H1 has no 5′ UTR, this matches the pipeline's CDS-mode output.
    """
    idx = SeqIO.index(str(_TESTDATA / "test_gnomad_gencode.fa"), "fasta")
    seq = idx[_GENCODE_OR11H1_KEY].seq
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coding_ref, _ = get_altseq(
            seq, Seq("A"), Seq("A"), 15528192, "+", OR11H1_EXON_FEAT, None
        )
        ref_orfs, _ = get_orfs_vcf(coding_ref, coding_ref, 1, 3)
    return str(ref_orfs[0])  # frame-1 canonical protein


class TestGnomADPipelineTranslation:
    """End-to-end pipeline correctness for gnomAD variant filtering and translation.

    The FASTA has CDS=1-948 → 1-frame CDS translation → one record per variant,
    no _1/_2/_3 frame suffixes in record IDs.
    """

    def test_pipeline_produces_output(self, gnomad_pipeline_records):
        assert gnomad_pipeline_records, "pipeline output must not be empty"

    def test_af_filtered_variant_absent(self, gnomad_pipeline_records):
        """rs1410655344 (missense, AF_afr=6e-05 < 0.01) must be absent."""
        assert not any("rs1410655344" in k for k in gnomad_pipeline_records)

    def test_consequence_filtered_variant_absent(self, gnomad_pipeline_records):
        """rs751806421 (stop_gained, AF_afr=0.02) must be absent — not in include_consequences."""
        assert not any("rs751806421" in k for k in gnomad_pipeline_records)

    def test_missense_produces_output(self, gnomad_pipeline_records):
        """rs199856986 (missense V78A, AF_afr=0.124) passes all filters."""
        assert any("rs199856986" in k for k in gnomad_pipeline_records)

    def test_inframe_deletion_produces_output(self, gnomad_pipeline_records):
        """rs1203023715 (inframe_del Phe242del, AF_afr=0.015) passes all filters."""
        assert any("rs1203023715" in k for k in gnomad_pipeline_records)

    def test_one_orf_emitted_per_variant(self, gnomad_pipeline_records):
        """CDS= header triggers 1-frame translation → exactly 1 record per passing variant."""
        missense_keys = [k for k in gnomad_pipeline_records if "rs199856986" in k]
        assert len(missense_keys) == 1, (
            f"expected 1 CDS-mode record for missense, got {len(missense_keys)}"
        )

    def test_missense_differs_at_exactly_one_position(
        self, gnomad_pipeline_records, gnomad_ref_protein
    ):
        """rs199856986 (V78A): alt protein differs from reference at exactly 1 position."""
        seqs = [s for k, s in gnomad_pipeline_records.items() if "rs199856986" in k]
        assert seqs, "no output for missense variant"
        alt_p = seqs[0]
        ref_stripped = gnomad_ref_protein.rstrip("*")
        alt_stripped = alt_p.rstrip("*")
        diffs = sum(1 for r, a in zip(ref_stripped, alt_stripped) if r != a)
        assert diffs == 1, f"expected 1 AA diff, got {diffs}"

    def test_missense_has_ala_at_position77(
        self, gnomad_pipeline_records, gnomad_ref_protein
    ):
        """rs199856986 (V78A): alt protein has A at position 77 (0-indexed); reference has V."""
        seqs = [s for k, s in gnomad_pipeline_records.items() if "rs199856986" in k]
        assert seqs, "no output for missense variant"
        assert gnomad_ref_protein[77] == "V", "reference must have V at AA 77"
        assert seqs[0][77] == "A", f"expected A at AA 77, got {seqs[0][77]}"

    def test_inframe_deletion_shorter(
        self, gnomad_pipeline_records, gnomad_ref_protein
    ):
        """rs1203023715 (Phe242del): alt protein is 1 AA shorter than reference."""
        seqs = [s for k, s in gnomad_pipeline_records.items() if "rs1203023715" in k]
        assert seqs, "no output for inframe_deletion variant"
        alt_p = seqs[0]
        ref_len = len(gnomad_ref_protein.rstrip("*"))
        alt_len = len(alt_p.rstrip("*"))
        assert alt_len == ref_len - 1, (
            f"expected ref_len-1={ref_len - 1}, got {alt_len}"
        )
