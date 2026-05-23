"""Translation-correctness tests for Ensembl vcf-to-proteindb.

Unit tests (TestGetAltseqOR11H1) call get_altseq() + get_orfs_vcf() directly
on OR11H1 (ENST00000643195, chr22 + strand, single-exon CDS=1-948), covering
all mutation types present in test_ensembl_v2p.vcf.

Integration tests (TestEnsemblPipelineTranslation) run _vcf_to_proteindb_chunk
end-to-end on those same files and assert that each output protein has the
correct biological properties.

OR11H1 genomic layout (GRCh38, chr22, + strand):
  CDS exon:   15528192–15529136  (945 bp → 315 coding codons)
  stop codon: 15529137–15529139  (3 bp)
  CDS=1-948 in gffread FASTA header (CDS + stop codon included)

Variants from test_ensembl_v2p.vcf used for unit tests:
  rs1211697244  22:15528192 A>G   start_lost       M→V   at AA 0
  rs1410655344  22:15528195 A>G   missense         N→D   at AA 1
  rs1394965478  22:15528197 T>C   synonymous       no protein change
  rs1402769459  22:15528197 TG>T  frameshift       diverges at AA 2
  rs1420478920  22:15528234 G>T   stop_gained      E→*   at AA 14
  rs1203023715  22:15528914 CTTC>C inframe_del     removes 1 AA (3 bp del)
  rs1986039639  22:15528961 G>GCTG inframe_ins     adds 1 AA (3 bp ins)
  rs1986046473  22:15529137 T>A   stop_lost        * removed, protein extends
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

# Genomic CDS features for OR11H1 ENST00000643195 (1-based, inclusive)
OR11H1_CDS_FEAT = [[15528192, 15529136, "CDS"]]
OR11H1_STOP_FEAT = [[15529137, 15529139, "stop_codon"]]
OR11H1_CDS_AND_STOP = OR11H1_CDS_FEAT + OR11H1_STOP_FEAT

# CDS positions in the gffread-generated FASTA (CDS=1-948 includes stop codon)
OR11H1_CDS_INFO = [1, 948]

# Absolute testdata directory — needed for module-scoped fixtures that run
# before the autouse _chdir_to_package_root function-scoped fixture fires.
_TESTDATA = Path(__file__).resolve().parent.parent / "testdata"


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def or11h1_seq():
    """OR11H1 CDS sequence from test_ensembl_v2p.fa (gffread, CDS=1-948)."""
    idx = SeqIO.index(str(_TESTDATA / "test_ensembl_v2p.fa"), "fasta")
    return idx["ENST00000643195"].seq


def _translate(seq: Seq) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return str(seq.translate(to_stop=False))


def _altseq_protein(ref_seq, ref_allele, alt_allele, pos, features=None, cds_info=None):
    """Return (ref_protein, alt_protein) for a variant on the + strand."""
    features = features if features is not None else OR11H1_CDS_FEAT
    cds_info = cds_info if cds_info is not None else OR11H1_CDS_INFO
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coding_ref, coding_alt = get_altseq(
            ref_seq, Seq(ref_allele), Seq(alt_allele), pos, "+", features, cds_info
        )
        ref_orfs, alt_orfs = get_orfs_vcf(coding_ref, coding_alt, 1, 1)
    return str(ref_orfs[0]), str(alt_orfs[0])


# ===========================================================================
# TestGetAltseqOR11H1 — unit tests per mutation type
# ===========================================================================

class TestGetAltseqOR11H1:
    """Direct get_altseq / get_orfs_vcf tests on real OR11H1 CDS coordinates."""

    def test_reference_protein_starts_with_met(self, or11h1_seq):
        """The reference translation begins with Met."""
        ref_p, _ = _altseq_protein(or11h1_seq, "A", "A", 15528192)
        assert ref_p[0] == "M", f"expected M, got {ref_p[0]}"

    def test_reference_protein_length(self, or11h1_seq):
        """CDS=1-948 → 948/3 = 316 codons (315 coding + 1 stop)."""
        ref_p, _ = _altseq_protein(or11h1_seq, "A", "A", 15528192)
        assert len(ref_p) == 316

    def test_reference_protein_ends_with_stop(self, or11h1_seq):
        """Canonical translation ends with '*' (stop codon)."""
        ref_p, _ = _altseq_protein(or11h1_seq, "A", "A", 15528192)
        assert ref_p[-1] == "*"

    # --- start_lost ---------------------------------------------------------

    def test_start_lost_changes_first_aa(self, or11h1_seq):
        """22:15528192 A>G (rs1211697244): Met codon ATG → GTG (Met → Val)."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528192)
        assert ref_p[0] == "M"
        assert alt_p[0] == "V"

    def test_start_lost_does_not_alter_length(self, or11h1_seq):
        """start_lost is a missense at AA 0 — protein length is unchanged."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528192)
        assert len(ref_p) == len(alt_p)

    def test_start_lost_downstream_aas_unchanged(self, or11h1_seq):
        """All residues after position 0 must be identical to reference."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528192)
        assert ref_p[1:] == alt_p[1:]

    # --- missense -----------------------------------------------------------

    def test_missense_changes_single_aa(self, or11h1_seq):
        """22:15528195 A>G (rs1410655344): codon 2 AAT→GAT (Asn→Asp)."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528195)
        assert ref_p[1] == "N"
        assert alt_p[1] == "D"

    def test_missense_length_unchanged(self, or11h1_seq):
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528195)
        assert len(ref_p) == len(alt_p)

    def test_missense_only_one_aa_differs(self, or11h1_seq):
        ref_p, alt_p = _altseq_protein(or11h1_seq, "A", "G", 15528195)
        diffs = sum(1 for r, a in zip(ref_p, alt_p) if r != a)
        assert diffs == 1

    # --- synonymous ---------------------------------------------------------

    def test_synonymous_leaves_protein_unchanged(self, or11h1_seq):
        """22:15528197 T>C (rs1394965478): synonymous — protein identical."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "T", "C", 15528197)
        assert ref_p == alt_p

    # --- stop_gained --------------------------------------------------------

    def test_stop_gained_inserts_stop_at_aa14(self, or11h1_seq):
        """22:15528234 G>T (rs1420478920): codon 15 GAG→TAG (Glu→*)."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "G", "T", 15528234)
        assert ref_p[14] == "E"
        assert alt_p[14] == "*"

    def test_stop_gained_aas_before_stop_unchanged(self, or11h1_seq):
        ref_p, alt_p = _altseq_protein(or11h1_seq, "G", "T", 15528234)
        assert ref_p[:14] == alt_p[:14]

    # --- frameshift ---------------------------------------------------------

    def test_frameshift_diverges_from_ref(self, or11h1_seq):
        """22:15528197 TG>T (rs1402769459): -1 deletion; frame shifts at AA 2."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "TG", "T", 15528197)
        assert ref_p[:2] == alt_p[:2], "first 2 AAs precede the deletion"
        assert ref_p[2] != alt_p[2], "AA 2 must change after frameshift"

    def test_frameshift_proteins_are_different(self, or11h1_seq):
        ref_p, alt_p = _altseq_protein(or11h1_seq, "TG", "T", 15528197)
        assert ref_p != alt_p

    # --- inframe_deletion ---------------------------------------------------

    def test_inframe_deletion_shortens_protein_by_one_aa(self, or11h1_seq):
        """22:15528914 CTTC>C (rs1203023715): 3-bp deletion removes 1 codon."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "CTTC", "C", 15528914)
        assert len(alt_p) == len(ref_p) - 1

    def test_inframe_deletion_preserves_upstream_aas(self, or11h1_seq):
        """Residues before the deletion are identical in ref and alt."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "CTTC", "C", 15528914)
        # deletion starts at CDS index 722 → AA 240
        del_aa = (15528914 - 15528192) // 3
        assert ref_p[:del_aa] == alt_p[:del_aa]

    # --- inframe_insertion --------------------------------------------------

    def test_inframe_insertion_lengthens_protein_by_one_aa(self, or11h1_seq):
        """22:15528961 G>GCTG (rs1986039639): 3-bp insertion adds 1 codon."""
        ref_p, alt_p = _altseq_protein(or11h1_seq, "G", "GCTG", 15528961)
        assert len(alt_p) == len(ref_p) + 1

    def test_inframe_insertion_preserves_upstream_aas(self, or11h1_seq):
        ref_p, alt_p = _altseq_protein(or11h1_seq, "G", "GCTG", 15528961)
        ins_aa = (15528961 - 15528192) // 3
        assert ref_p[:ins_aa] == alt_p[:ins_aa]

    # --- stop_lost ----------------------------------------------------------

    def test_stop_lost_removes_stop_codon(self, or11h1_seq):
        """22:15529137 T>A (rs1986046473): stop codon TAA→AAA; '*' disappears."""
        ref_p, alt_p = _altseq_protein(
            or11h1_seq, "T", "A", 15529137,
            features=OR11H1_CDS_AND_STOP,
        )
        assert ref_p[-1] == "*", "reference must end with stop"
        assert "*" not in alt_p, "stop_lost alt protein must have no stop codon"

    def test_stop_lost_alt_protein_extends_past_normal_terminus(self, or11h1_seq):
        """Alt protein is longer than reference (stop removed, reads through)."""
        ref_p, alt_p = _altseq_protein(
            or11h1_seq, "T", "A", 15529137,
            features=OR11H1_CDS_AND_STOP,
        )
        assert len(alt_p) > len(ref_p) - 1, "alt must extend past normal stop"


# ===========================================================================
# TestEnsemblPipelineTranslation — full pipeline integration
# ===========================================================================

_V2P_PIPELINE_ARGS = {
    "annotation_field_name": "CSQ",
    "transcript_str": "Feature",
    "consequence_str": "Consequence",
    "biotype_str": "feature_type",
    "include_biotypes": "mRNA,ncRNA",
    "af_field": "MAF",
    "af_threshold": 0.001,
    "protein_prefix": "var",
}


@pytest.fixture(scope="module")
def ensembl_pipeline_records(tmp_path_factory):
    """Run _vcf_to_proteindb_chunk on OR11H1 testdata; return parsed records."""
    from pgatk.ensembl.ensembl import EnsemblDataService
    out = str(tmp_path_factory.mktemp("ensembl_translation") / "out.fa")
    svc = EnsemblDataService({}, _V2P_PIPELINE_ARGS)
    svc._vcf_to_proteindb_chunk(
        str(_TESTDATA / "test_ensembl_v2p.vcf"),
        str(_TESTDATA / "test_ensembl_v2p.fa"),
        str(_TESTDATA / "test_ensembl_v2p.gtf"),
        out,
    )
    records = {}
    if Path(out).exists():
        for rec in SeqIO.parse(out, "fasta"):
            records[rec.id] = str(rec.seq)
    return records


@pytest.fixture(scope="module")
def or11h1_ref_protein():
    """Reference OR11H1 protein computed with identity variant."""
    idx = SeqIO.index(str(_TESTDATA / "test_ensembl_v2p.fa"), "fasta")
    seq = idx["ENST00000643195"].seq
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coding_ref, _ = get_altseq(seq, Seq("A"), Seq("A"), 15528192, "+",
                                   OR11H1_CDS_FEAT, OR11H1_CDS_INFO)
        ref_orfs, _ = get_orfs_vcf(coding_ref, coding_ref, 1, 1)
    return str(ref_orfs[0])


class TestEnsemblPipelineTranslation:
    """End-to-end pipeline correctness for each OR11H1 mutation type."""

    def test_pipeline_produces_output(self, ensembl_pipeline_records):
        assert ensembl_pipeline_records, "pipeline output must not be empty"

    def test_synonymous_variant_absent_from_output(self, ensembl_pipeline_records):
        """rs1394965478 (synonymous) should not appear — identical to reference."""
        assert not any("rs1394965478" in k for k in ensembl_pipeline_records)

    def test_missense_produces_output(self, ensembl_pipeline_records):
        assert any("rs1410655344" in k for k in ensembl_pipeline_records)

    def test_stop_gained_produces_output(self, ensembl_pipeline_records):
        assert any("rs1420478920" in k for k in ensembl_pipeline_records)

    def test_frameshift_produces_output(self, ensembl_pipeline_records):
        assert any("rs1402769459" in k for k in ensembl_pipeline_records)

    def test_inframe_deletion_produces_output(self, ensembl_pipeline_records):
        assert any("rs1203023715" in k for k in ensembl_pipeline_records)

    def test_inframe_insertion_produces_output(self, ensembl_pipeline_records):
        assert any("rs1986039639" in k for k in ensembl_pipeline_records)

    def test_stop_lost_produces_output(self, ensembl_pipeline_records):
        assert any("rs1986046473" in k for k in ensembl_pipeline_records)

    def test_missense_differs_at_exactly_one_position(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1410655344 (N→D) must change exactly 1 AA."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1410655344" in k]
        assert seqs, "no output for missense variant"
        alt_p = seqs[0]
        ref_stripped = or11h1_ref_protein.rstrip("*")
        alt_stripped = alt_p.rstrip("*")
        diffs = sum(1 for r, a in zip(ref_stripped, alt_stripped) if r != a)
        assert diffs == 1, f"expected 1 AA diff, got {diffs}"

    def test_stop_gained_introduces_early_stop(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1420478920 (E→* at AA 14): alt protein has '*' at index 14."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1420478920" in k]
        assert seqs, "no output for stop_gained variant"
        alt_p = seqs[0]
        assert alt_p[14] == "*", f"expected * at AA 14, got {alt_p[14]}"
        assert or11h1_ref_protein[14] == "E"

    def test_frameshift_diverges_from_ref(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1402769459 (−1 del): alt and ref diverge from AA 2 onward."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1402769459" in k]
        assert seqs, "no output for frameshift variant"
        alt_p = seqs[0]
        assert alt_p[:2] == or11h1_ref_protein[:2], "first 2 AAs must match"
        assert alt_p[2] != or11h1_ref_protein[2], "AA 2 must diverge"

    def test_inframe_deletion_alt_protein_shorter(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1203023715 (3-bp del): alt protein is 1 AA shorter than reference."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1203023715" in k]
        assert seqs, "no output for inframe_deletion variant"
        alt_p = seqs[0]
        ref_len = len(or11h1_ref_protein.rstrip("*"))
        alt_len = len(alt_p.rstrip("*"))
        assert alt_len == ref_len - 1, (
            f"expected ref_len-1={ref_len - 1}, got {alt_len}"
        )

    def test_inframe_insertion_alt_protein_longer(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1986039639 (3-bp ins): alt protein is 1 AA longer than reference."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1986039639" in k]
        assert seqs, "no output for inframe_insertion variant"
        alt_p = seqs[0]
        ref_len = len(or11h1_ref_protein.rstrip("*"))
        alt_len = len(alt_p.rstrip("*"))
        assert alt_len == ref_len + 1, (
            f"expected ref_len+1={ref_len + 1}, got {alt_len}"
        )

    def test_stop_lost_alt_protein_extends_past_stop(
        self, ensembl_pipeline_records, or11h1_ref_protein
    ):
        """rs1986046473 (stop_lost): alt protein extends past normal terminus."""
        seqs = [s for k, s in ensembl_pipeline_records.items() if "rs1986046473" in k]
        assert seqs, "no output for stop_lost variant"
        alt_p = seqs[0]
        assert len(alt_p) > len(or11h1_ref_protein) - 1, (
            "stop_lost alt must be longer than canonical coding sequence"
        )
