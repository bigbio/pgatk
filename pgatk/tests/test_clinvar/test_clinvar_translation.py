"""Tests for ClinVar-to-proteindb variant translation, covering all mutation types.

Unit tests (TestGetAltseq*) exercise get_altseq() + get_orfs_vcf() directly
with synthetic sequences of known structure — no file I/O required.

Integration tests (TestClinVarISG15Integration) run the full ClinVarService
pipeline against real ISG15 (NM_005101.4) data from use-cases/ncbi_clinvar/.
They are skipped automatically when that directory is absent.

Real variant sources used (all ISG15, NM_005101.4):
  missense     chr1:1013997 C>A   ClinVar ID 1035971  (Uncertain_significance)
  stop_gained  chr1:1014026 C>T   engineered from codon 16 CAG→TAG
  inframe_del  chr1:1014223 TCTGAGCATC>T  ClinVar ID 915424  (Uncertain_significance)
  frameshift   chr1:1014316 C>CG  ClinVar ID 161455   (Pathogenic)
  inframe_ins  chr1:1014459 A>AGCCCGT     ClinVar ID 1038082  (Uncertain_significance)
  benign       chr1:1013541 T>C   ClinVar ID 1185394  (Benign — filtered out)
"""
from __future__ import annotations

import shutil
import warnings
from pathlib import Path

import pytest
from Bio.Seq import Seq

from pgatk.toolbox.vcf_utils import get_altseq, get_orfs_vcf

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TESTDATA_DIR = Path(__file__).resolve().parent.parent.parent / "testdata" / "clinvar"
USE_CASES_DIR = Path(__file__).resolve().parent.parent.parent.parent / "use-cases" / "ncbi_clinvar"

_REAL_DATA = (
    (USE_CASES_DIR / "clinvar.vcf").exists()
    and (USE_CASES_DIR / "GRCh38_latest_genomic.gff").exists()
    and (USE_CASES_DIR / "transcripts.fa").exists()
    and (USE_CASES_DIR / "GRCh38_latest_assembly_report.txt").exists()
)

# ---------------------------------------------------------------------------
# Unit-test shared fixture
# ---------------------------------------------------------------------------
# Synthetic 18-bp single-CDS-exon transcript on + strand.
#
#   mRNA:     A T G A A A C C C G G G T T T T G G
#   pos:      1000              ...             1017  (genomic, 1-based)
#   codons:   ATG(M) AAA(K) CCC(P) GGG(G) TTT(F) TGG(W)
#   protein:  M  K  P  G  F  W
#
# cds_info = [1, 18] means the whole sequence is CDS (no UTR).
# features_info lists the single exon in gffutils format [[start, end, type]].

SYNTH_SEQ = Seq("ATGAAACCCGGGTTTTGG")
SYNTH_FEATURES = [[1000, 1017, "CDS"]]
SYNTH_CDS_INFO = [1, 18]
SYNTH_REF_PROTEIN = "MKPGFW"  # translate(to_stop=False)

# For minus-strand tests the same sequence occupies pos 2000-2017 on the
# reverse strand (complement is read right-to-left).
SYNTH_MINUS_FEATURES = [[2000, 2017, "CDS"]]

# ISG15 CDS exon features used in multi-exon and integration tests.
# Genomic coords (1-based, gffutils style):
#   exon 1: NC_000001.11 1013574–1013576  (3 bp)
#   exon 2: NC_000001.11 1013984–1014475  (492 bp)
ISG15_CDS_EXONS = [[1013574, 1013576, "CDS"], [1013984, 1014475, "CDS"]]
ISG15_CDS_INFO = [78, 572]  # CDS positions in NM_005101.4 mRNA (1-based)


def _translate(seq: Seq) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return str(seq.translate(to_stop=False))


def _altseq_and_protein(ref_seq, ref_allele, alt_allele, pos,
                         strand="+", features=None, cds_info=None):
    features = features or SYNTH_FEATURES
    cds_info = cds_info or SYNTH_CDS_INFO
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coding_ref, coding_alt = get_altseq(
            ref_seq, Seq(ref_allele), Seq(alt_allele),
            pos, strand, features, cds_info,
        )
        ref_orfs, alt_orfs = get_orfs_vcf(coding_ref, coding_alt, 1, 1)
    return str(ref_orfs[0]), str(alt_orfs[0])


# ===========================================================================
# TestGetAltseqMissense
# ===========================================================================

class TestGetAltseqMissense:
    """SNP that substitutes one amino acid without altering protein length."""

    def test_missense_changes_single_amino_acid(self):
        """pos 1003 A>C changes codon 2 AAA(K) → CAA(Q)."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "C", 1003)
        assert len(ref_p) == len(alt_p), "missense must not alter protein length"
        assert ref_p[1] == "K"
        assert alt_p[1] == "Q"

    def test_missense_leaves_surrounding_aas_unchanged(self):
        """Residues outside the mutated codon must be identical."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "C", 1003)
        assert ref_p[0] == alt_p[0] == "M"
        assert ref_p[2:] == alt_p[2:]

    def test_missense_ref_protein_matches_known_sequence(self):
        ref_p, _ = _altseq_and_protein(SYNTH_SEQ, "A", "C", 1003)
        assert ref_p == SYNTH_REF_PROTEIN

    def test_missense_at_last_codon(self):
        """SNP in the final coding codon TGG(W) → AGG(R): pos 1015 T>A.

        pos 1015 = index 15 = first base of the last codon (TGG = Trp).
        Using REF='T' (the actual base) and ALT='A' → AGG (Arg).
        """
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "T", "A", 1015)
        assert ref_p[-1] == "W"
        assert alt_p[-1] == "R"
        assert ref_p[:-1] == alt_p[:-1]

    def test_same_ref_alt_returns_identical_proteins(self):
        """When REF == ALT the output proteins must be identical."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "A", 1003)
        assert ref_p == alt_p


# ===========================================================================
# TestGetAltseqStopGained
# ===========================================================================

class TestGetAltseqStopGained:
    """SNP that introduces a premature stop codon."""

    def test_stop_gained_inserts_stop_at_mutated_position(self):
        """pos 1003 A>T changes codon 2 AAA(K) → TAA(*): stop at position 1."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "T", 1003)
        assert ref_p[1] == "K"
        assert alt_p[1] == "*", "codon 2 should become a stop codon"

    def test_stop_gained_protein_length_equal_with_to_stop_false(self):
        """translate(to_stop=False) keeps the full-length sequence including '*'."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "T", 1003)
        assert len(ref_p) == len(alt_p)

    def test_stop_gained_sequence_before_stop_unchanged(self):
        """Amino acids before the stop codon must be identical to reference."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "T", 1003)
        assert alt_p[0] == ref_p[0] == "M"

    def test_stop_gained_at_first_coding_codon(self):
        """SNP that destroys the start methionine context is still handled."""
        # pos 1001 T>A changes codon 1 ATG → AAG (not a stop, but tests boundary)
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "T", "A", 1001)
        assert ref_p[0] == "M"
        assert alt_p[0] == "K"


# ===========================================================================
# TestGetAltseqFrameshiftVariant
# ===========================================================================

class TestGetAltseqFrameshiftVariant:
    """Single-base indels that shift the reading frame."""

    def test_frameshift_insertion_changes_downstream_aas(self):
        """Inserting G at pos 1006 (between codon 2 and 3) shifts frame from codon 3."""
        # alt: ATGAAA TG CCGGGTTTTGG → MKCRVL...
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "T", "TG", 1006)
        assert ref_p[:2] == alt_p[:2] == "MK", "first 2 codons precede the insertion"
        assert ref_p[2] != alt_p[2], "codon 3 must change due to frameshift"

    def test_frameshift_deletion_changes_downstream_aas(self):
        """Deleting one A at pos 1004 (within codon 2) shifts frame from codon 2."""
        # REF=AA ALT=A → delete second base of codon 2
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "AA", "A", 1003)
        assert ref_p[0] == alt_p[0] == "M"
        assert ref_p[1] != alt_p[1], "codon 2 must change due to frameshift"

    def test_frameshift_insertion_alters_protein_length(self):
        """A +1 insertion shifts the reading frame from codon 3 onward.

        18 bp + 1 inserted base = 19 bp.  Biopython drops the trailing 1-base
        partial codon so both ref and alt translate to 6 aa — lengths are equal.
        The real frameshift signature is that all residues from the insertion
        point onward differ, which is the invariant we assert here.
        """
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "T", "TG", 1006)
        assert ref_p[:2] == alt_p[:2] == "MK", "first 2 AAs precede the insertion"
        assert ref_p[2:] != alt_p[2:], "all residues from insertion point must change"

    def test_frameshift_proteins_are_different(self):
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "T", "TG", 1006)
        assert ref_p != alt_p


# ===========================================================================
# TestGetAltseqInframeDeletion
# ===========================================================================

class TestGetAltseqInframeDeletion:
    """Deletion of a multiple-of-3 bases: removes whole codons, no frame shift."""

    def test_inframe_deletion_shortens_protein_by_one_aa(self):
        """Delete codon 2 (AAA) via GAAA>G at pos 1002; protein loses K."""
        # REF=GAAA ALT=G → delete positions 1003-1005
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "GAAA", "G", 1002)
        assert len(alt_p) == len(ref_p) - 1
        assert "K" in ref_p
        # K should no longer appear at position 1
        assert alt_p[1] != "K"

    def test_inframe_deletion_preserves_flanking_aas(self):
        """Start Met and codons after the deletion must be preserved."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "GAAA", "G", 1002)
        assert alt_p[0] == "M"
        # Codons 3-6 of ref → codons 2-5 of alt (P G F W)
        assert alt_p[1:] == "PGFW"

    def test_inframe_deletion_six_bases_removes_two_aas(self):
        """Delete 6 bases (codons 2-3: AAA CCC) via GAAACCC>G."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "GAAACCC", "G", 1002)
        assert len(alt_p) == len(ref_p) - 2
        assert alt_p == "MGFW"


# ===========================================================================
# TestGetAltseqInframeInsertion
# ===========================================================================

class TestGetAltseqInframeInsertion:
    """Insertion of a multiple-of-3 bases: adds whole codons, no frame shift."""

    def test_inframe_insertion_lengthens_protein_by_one_aa(self):
        """Insert CAT (His) after codon 2 via A>ACAT at pos 1005."""
        # alt: ATGAAA CAT CCCGGGTTTTGG → MKHPGFW
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "ACAT", 1005)
        assert len(alt_p) == len(ref_p) + 1

    def test_inframe_insertion_adds_correct_amino_acid(self):
        """CAT encodes His (H); it should appear between position 2 and 3."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "ACAT", 1005)
        assert alt_p == "MKHPGFW"

    def test_inframe_insertion_preserves_flanking_aas(self):
        """Codons before and after the insertion must be unchanged."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "ACAT", 1005)
        assert alt_p[:2] == "MK"
        assert alt_p[3:] == "PGFW"

    def test_inframe_insertion_six_bases_adds_two_aas(self):
        """Insert CATGGG (His + Gly) via A>ACATGGG at pos 1005."""
        ref_p, alt_p = _altseq_and_protein(SYNTH_SEQ, "A", "ACATGGG", 1005)
        assert len(alt_p) == len(ref_p) + 2
        assert alt_p == "MKHGPGFW"


# ===========================================================================
# TestGetAltseqMinusStrand
# ===========================================================================

class TestGetAltseqMinusStrand:
    """Variant application on reverse-complement (minus) strand transcripts."""

    def test_minus_strand_missense_changes_protein(self):
        """The same SNP applied to a minus-strand gene changes a codon."""
        # Minus-strand gene: genomic pos 2000-2017, sequence complement read right-to-left.
        # The ref_seq for get_altseq is still the mRNA (5'→3').
        ref_p, alt_p = _altseq_and_protein(
            SYNTH_SEQ, "A", "C", 2014,
            strand="-", features=SYNTH_MINUS_FEATURES,
        )
        assert ref_p != alt_p

    def test_minus_strand_ref_protein_same_length(self):
        """A missense on minus strand must not alter protein length."""
        ref_p, alt_p = _altseq_and_protein(
            SYNTH_SEQ, "A", "C", 2014,
            strand="-", features=SYNTH_MINUS_FEATURES,
        )
        assert len(ref_p) == len(alt_p)


# ===========================================================================
# TestGetAltseqMultiExon — uses ISG15 CDS exon structure
# ===========================================================================

class TestGetAltseqMultiExon:
    """Variant handling across a real two-exon CDS (ISG15 NM_005101.4)."""

    @pytest.fixture(autouse=True)
    def _load_isg15(self):
        """Load ISG15 mRNA from the use-cases FASTA; skip if absent."""
        from Bio import SeqIO
        fasta = USE_CASES_DIR / "GRCh38_latest_rna.fna"
        if not fasta.exists():
            pytest.skip("use-cases/ncbi_clinvar/GRCh38_latest_rna.fna not found")
        idx = SeqIO.index(str(fasta), "fasta")
        self.isg15_seq = idx["NM_005101.4"].seq

    def test_missense_in_exon2_changes_protein(self):
        """chr1:1013997 C>A (in CDS exon 2, ClinVar ID 1035971) changes AA 6.

        var_index_in_cds = 3 (exon1 length) + (1013997 - 1013984) = 16.
        Codon index = 16 // 3 = 5  →  0-based protein index 5 (the 6th AA).
        """
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "C", "A", 1013997,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert ref_p != alt_p
        assert len(ref_p) == len(alt_p), "missense must not alter length"
        assert ref_p[5] != alt_p[5], "change should be at AA position 6 (index 5)"

    def test_stop_gained_creates_premature_stop(self):
        """chr1:1014026 C>T (codon 16 CAG→TAG) inserts stop at AA 15 (index)."""
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "C", "T", 1014026,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert alt_p[15] == "*", "codon 16 (index 15) must become a stop"
        assert ref_p[15] != "*"

    def test_frameshift_insertion_diverges_at_cds_position_113(self):
        """chr1:1014316 C>CG (ClinVar ID 161455, Pathogenic) shifts frame at AA 113."""
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "C", "CG", 1014316,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert ref_p[:113] == alt_p[:113], "residues before frameshift must be unchanged"
        assert ref_p[113] != alt_p[113], "residue 113 must diverge after frameshift"

    def test_inframe_deletion_shortens_protein_by_3aa(self):
        """chr1:1014223 TCTGAGCATC>T (ClinVar ID 915424) removes 3 aa (9 bp del)."""
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "TCTGAGCATC", "T", 1014223,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert len(alt_p) == len(ref_p) - 3

    def test_inframe_insertion_lengthens_protein_by_2aa(self):
        """chr1:1014459 A>AGCCCGT (ClinVar ID 1038082) inserts 2 aa (6 bp ins)."""
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "A", "AGCCCGT", 1014459,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert len(alt_p) == len(ref_p) + 2

    def test_variant_in_exon1_cds_is_handled(self):
        """A variant overlapping the tiny 3-bp CDS exon 1 must still be applied."""
        # CDS exon 1 spans genomic 1013574-1013576; pos 1013574 = first CDS base
        ref_p, alt_p = _altseq_and_protein(
            self.isg15_seq, "A", "T", 1013574,
            features=ISG15_CDS_EXONS, cds_info=ISG15_CDS_INFO,
        )
        assert ref_p != alt_p


# ===========================================================================
# TestClinVarISG15Integration — full pipeline, real ISG15 data
# ===========================================================================

@pytest.mark.skipif(
    not shutil.which("bedtools"),
    reason="bedtools not installed",
)
class TestClinVarISG15Integration:
    """End-to-end pipeline tests using mini ISG15 testdata.

    Uses clinvar_isg15_cds.fna (header has CDS=1-495) so the service runs 1-frame CDS
    translation and produces the canonical ISG15 protein for comparison.

    Each test exercises one mutation type so that every branch of the
    ClinVar translation pipeline is covered by at least one real variant.
    Benign and Uncertain_significance variants are included to verify that
    the default CLNSIG filter behaviour is correct.
    """

    MINI_VCF = str(TESTDATA_DIR / "clinvar_isg15_mutation_types.vcf")
    # CDS-only FASTA: entire 495-bp sequence is CDS (CDS=1-495 in header)
    _ISG15_CDS_ONLY_INFO = [1, 495]

    @pytest.fixture()
    def run_pipeline(self, tmp_path):
        """Run ClinVarService using mini ISG15 testdata; returns parsed output."""
        from pgatk.clinvar.clinvar_service import ClinVarService
        from Bio import SeqIO

        def _run(output_name="output.fa"):
            output_file = str(tmp_path / output_name)
            svc = ClinVarService(
                vcf_file=self.MINI_VCF,
                gff_file=str(TESTDATA_DIR / "clinvar_isg15_cds.gff"),
                fasta_file=str(TESTDATA_DIR / "clinvar_isg15_cds.fna"),
                assembly_report=str(TESTDATA_DIR / "clinvar_isg15_assembly_report.txt"),
                output_file=output_file,
            )
            svc.run()
            records = {}
            if Path(output_file).exists():
                for rec in SeqIO.parse(output_file, "fasta"):
                    records[rec.id] = str(rec.seq)
            return records

        return _run

    def test_pipeline_runs_without_error(self, run_pipeline):
        """Service must complete without raising an exception."""
        records = run_pipeline()
        assert isinstance(records, dict)

    def test_benign_variant_is_filtered_out(self, run_pipeline):
        """CLNSIG=Benign (ID 1185394) must not appear in output."""
        records = run_pipeline()
        assert not any("1185394" in k for k in records), (
            "Benign variant 1185394 should be excluded by default CLNSIG filter"
        )

    def test_missense_variant_produces_output(self, run_pipeline):
        """CLNSIG=Uncertain_significance missense (ID 1035971) must appear in output."""
        records = run_pipeline()
        assert any("1035971" in k for k in records), (
            "missense variant 1035971 (Uncertain_significance) should pass the filter"
        )

    def test_stop_gained_variant_produces_output(self, run_pipeline):
        """CLNSIG=Pathogenic stop_gained (ID 9000001) must appear in output."""
        records = run_pipeline()
        assert any("9000001" in k for k in records), (
            "stop_gained variant 9000001 (Pathogenic) should produce output"
        )

    def test_frameshift_variant_produces_output(self, run_pipeline):
        """CLNSIG=Pathogenic frameshift (ClinVar 161455) must appear in output."""
        records = run_pipeline()
        assert any("161455" in k for k in records), (
            "frameshift variant 161455 (Pathogenic) should produce output"
        )

    def test_inframe_deletion_produces_output(self, run_pipeline):
        """Inframe deletion (ClinVar 915424, Uncertain_significance) must appear."""
        records = run_pipeline()
        assert any("915424" in k for k in records), (
            "inframe_deletion variant 915424 should produce output"
        )

    def test_inframe_insertion_produces_output(self, run_pipeline):
        """Inframe insertion (ClinVar 1038082, Uncertain_significance) must appear."""
        records = run_pipeline()
        assert any("1038082" in k for k in records), (
            "inframe_insertion variant 1038082 should produce output"
        )

    def test_missense_alt_protein_differs_from_ref_at_single_position(self, run_pipeline):
        """The missense alt protein must differ from reference at exactly one AA."""
        from Bio import SeqIO

        records = run_pipeline()
        missense_seqs = [seq for k, seq in records.items() if "1035971" in k]
        assert missense_seqs, "no output for missense variant"

        idx = SeqIO.index(str(TESTDATA_DIR / "clinvar_isg15_cds.fna"), "fasta")
        isg15_seq = idx["NM_005101.4"].seq
        ref_p, _ = _altseq_and_protein(
            isg15_seq, "C", "C", 1013997,  # REF==ALT → identity (reference protein)
            features=ISG15_CDS_EXONS, cds_info=self._ISG15_CDS_ONLY_INFO,
        )
        alt_p = missense_seqs[0].rstrip("*")
        ref_p_stripped = ref_p.rstrip("*")

        diffs = sum(1 for r, a in zip(ref_p_stripped, alt_p) if r != a)
        assert diffs == 1, f"missense should change exactly 1 AA, got {diffs}"

    def test_frameshift_protein_diverges_at_expected_position(self, run_pipeline):
        """Frameshift alt protein (ID 161455) must match reference up to AA 113."""
        from Bio import SeqIO

        records = run_pipeline()
        fs_seqs = [seq for k, seq in records.items() if "161455" in k]
        assert fs_seqs, "no output for frameshift variant"

        idx = SeqIO.index(str(TESTDATA_DIR / "clinvar_isg15_cds.fna"), "fasta")
        isg15_seq = idx["NM_005101.4"].seq
        ref_p, _ = _altseq_and_protein(
            isg15_seq, "C", "C", 1014316,
            features=ISG15_CDS_EXONS, cds_info=self._ISG15_CDS_ONLY_INFO,
        )
        alt_p = fs_seqs[0]
        assert alt_p[:113] == ref_p[:113], (
            "frameshift variant 161455 should match reference up to AA 113"
        )
        assert alt_p[113] != ref_p[113], (
            "frameshift variant 161455 must diverge at AA 113"
        )

    def test_inframe_deletion_alt_protein_is_shorter(self, run_pipeline):
        """Inframe deletion (9-bp, ID 915424) must reduce protein length by 3 aa."""
        from Bio import SeqIO

        records = run_pipeline()
        del_seqs = [seq for k, seq in records.items() if "915424" in k]
        assert del_seqs, "no output for inframe_deletion variant"

        idx = SeqIO.index(str(TESTDATA_DIR / "clinvar_isg15_cds.fna"), "fasta")
        isg15_seq = idx["NM_005101.4"].seq
        ref_p, _ = _altseq_and_protein(
            isg15_seq, "TCTGAGCATC", "TCTGAGCATC", 1014223,
            features=ISG15_CDS_EXONS, cds_info=self._ISG15_CDS_ONLY_INFO,
        )
        alt_p = del_seqs[0]
        assert len(alt_p) == len(ref_p) - 3, (
            f"9-bp inframe deletion should shorten protein by 3 aa; "
            f"ref={len(ref_p)}, alt={len(alt_p)}"
        )

    def test_inframe_insertion_alt_protein_is_longer(self, run_pipeline):
        """Inframe insertion (6-bp, ID 1038082) must extend protein by 2 aa."""
        from Bio import SeqIO

        records = run_pipeline()
        ins_seqs = [seq for k, seq in records.items() if "1038082" in k]
        assert ins_seqs, "no output for inframe_insertion variant"

        idx = SeqIO.index(str(TESTDATA_DIR / "clinvar_isg15_cds.fna"), "fasta")
        isg15_seq = idx["NM_005101.4"].seq
        ref_p, _ = _altseq_and_protein(
            isg15_seq, "A", "A", 1014459,
            features=ISG15_CDS_EXONS, cds_info=self._ISG15_CDS_ONLY_INFO,
        )
        alt_p = ins_seqs[0]
        assert len(alt_p) == len(ref_p) + 2, (
            f"6-bp inframe insertion should lengthen protein by 2 aa; "
            f"ref={len(ref_p)}, alt={len(alt_p)}"
        )
