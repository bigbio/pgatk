"""Tests for cbioportal_to_proteindb biotype filtering and improved translation.

Covers:
- _parse_cbio_fasta_header: CDS= and gene_biotype= extraction
- Variant_Classification filtering (Nonsense_Mutation excluded by default)
- Biotype filtering (protein_coding default, lncRNA excluded)
- 3-frame translation when include_biotypes=all for non-CDS transcripts
- CDS slice extraction: HGVSc positions are CDS-relative after slicing
"""

from pathlib import Path

import pytest
from Bio import SeqIO

from pgatk.cgenomes.cgenomes_proteindb import CancerGenomesService

_TESTDATA = Path(__file__).resolve().parent.parent / "testdata"


def _translate_ref_cds(fasta_path, nm_id: str) -> str:
    """Translate the CDS of an NCBI gffread transcript record.

    Strips the 'rna-' prefix and version suffix to match bare NM_ IDs.
    Returns an empty string when the transcript is not found.
    """
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        if rec.id.removeprefix("rna-").split(".")[0] != nm_id:
            continue
        for tok in rec.description.split():
            stripped = tok.strip("[]")
            if stripped.startswith("CDS="):
                s, e = (int(x) for x in stripped.split("=")[1].split("-"))
                return str(rec.seq[s - 1:e].translate())
    return ""

# Minimal transcript FASTA (gffread -F -w format):
#   - NM_PROT01: protein_coding, CDS=4-27 in a 30-base transcript
#     full seq:  AAA  ATG AAA CGT TTT AAG GGG TGG TAA  TTT
#     positions: 1-3  4----------------------------27  28-30
#     CDS protein reference: M-K-R-F-K-G-W-* (MKRFKGW*, length 8)
#   - NR_LNCRNA1: lncRNA, no CDS= header
#     seq: ATGAAATTTCCCGGGCAGTTTTAA (24 bases)
#     position 3 = G; c.3G>T produces ATTAAATTTCCCGGGCAGTTTTAA
_FASTA = """\
>rna-NM_PROT01.1 CDS=4-27 Dbxref=X;gene_biotype=protein_coding;transcript_id=NM_PROT01.1
AAAATGAAACGTTTTAAGGGGTGGTAATTT
>rna-NR_LNCRNA1.1 Dbxref=X;gene_biotype=lncRNA;transcript_id=NR_LNCRNA1.1
ATGAAATTTCCCGGGCAGTTTTAA
"""

# MAF-style mutation table (only columns required by cbioportal_to_proteindb)
# c.7C>T on NM_PROT01 (CDS pos 7, C→T): CGT→TGT, R→C → protein MKCFKGW*
# Same position but Nonsense_Mutation: tests default exclusion
# lncRNA mutation c.3G>T: tests biotype filter (excluded by default)
# lncRNA mutation c.3G>T again: used in include_biotypes=all test
_MAF_HEADER = "Hugo_Symbol\tTranscript_ID\tHGVSc\tHGVSp_Short\tVariant_Type\tVariant_Classification\tTumor_Sample_Barcode"
_MISSENSE_ROW = "GENE1\tNM_PROT01\tNM_PROT01.1:c.7C>T\tp.R3C\tSNP\tMissense_Mutation\tSAMPLE1"
_NONSENSE_ROW = "GENE1\tNM_PROT01\tNM_PROT01.1:c.7C>T\tp.R3*\tSNP\tNonsense_Mutation\tSAMPLE1"
_LNCRNA_ROW = "GENE2\tNR_LNCRNA1\tNR_LNCRNA1.1:c.3G>T\tp.?\tSNP\tMissense_Mutation\tSAMPLE1"


def _make_service(tmp_path, fasta_text, maf_rows, extra_args=None):
    """Write temp FASTA + MAF, return a configured CancerGenomesService."""
    fasta_path = tmp_path / "transcripts.fa"
    fasta_path.write_text(fasta_text)

    maf_path = tmp_path / "mutations.txt"
    maf_path.write_text("\n".join(maf_rows) + "\n")

    output_path = tmp_path / "output.fa"

    pipeline_args = {
        CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: str(maf_path),
        CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: str(fasta_path),
        CancerGenomesService.CONFIG_OUTPUT_FILE: str(output_path),
    }
    if extra_args:
        pipeline_args.update(extra_args)

    return CancerGenomesService({}, pipeline_args), output_path


class TestParseCbioFastaHeader:
    """Unit tests for the CDS= / biotype header parser."""

    def test_protein_coding_cds_present(self):
        desc = "rna-NM_PROT01.1 CDS=4-27 Dbxref=X;gene_biotype=protein_coding;transcript_id=NM_PROT01.1"
        cds_info, biotype, _ = CancerGenomesService._parse_cbio_fasta_header(desc, "gene_biotype")
        assert cds_info == [4, 27]
        assert biotype == "protein_coding"

    def test_lncrna_no_cds(self):
        desc = "rna-NR_LNCRNA1.1 Dbxref=X;gene_biotype=lncRNA;transcript_id=NR_LNCRNA1.1"
        cds_info, biotype, _ = CancerGenomesService._parse_cbio_fasta_header(desc, "gene_biotype")
        assert cds_info == []
        assert biotype == "lncRNA"

    def test_missing_biotype(self):
        desc = "rna-NM_NOBIOTYPE.1 CDS=1-9"
        cds_info, biotype, _ = CancerGenomesService._parse_cbio_fasta_header(desc, "gene_biotype")
        assert cds_info == [1, 9]
        assert biotype == ""

    def test_no_cds_no_biotype(self):
        desc = "some_transcript_without_annotations"
        cds_info, biotype, _ = CancerGenomesService._parse_cbio_fasta_header(desc, "gene_biotype")
        assert cds_info == []
        assert biotype == ""


class TestCbioportalToProteindb:
    """Integration tests for cbioportal_to_proteindb with biotype + classification filters."""

    def test_missense_protein_coding_included(self, tmp_path):
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _MISSENSE_ROW],
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        # Missense on protein_coding transcript must be translated
        assert "cbiomut:NM_PROT01:GENE1:p.R3C:Missense_Mutation" in output
        # No reading-frame suffix for CDS transcripts (1-frame)
        assert "_RF1" not in output

    def test_missense_protein_correct_substitution(self, tmp_path):
        """c.7C>T on CDS=4-27: codon 3 CGT→TGT gives R→C; expect MKCFKGW*."""
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _MISSENSE_ROW],
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        lines = output.strip().split("\n")
        # Find the protein sequence for the missense entry
        for i, line in enumerate(lines):
            if "cbiomut:NM_PROT01" in line:
                protein = lines[i + 1]
                # Position 3 should be C (cysteine), not R (arginine)
                assert protein[2] == "C", f"Expected C at pos 3, got: {protein}"
                assert protein.startswith("MK")
                break
        else:
            pytest.fail("Missense entry not found in output")

    def test_nonsense_excluded_by_default(self, tmp_path):
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _NONSENSE_ROW],
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        assert "Nonsense_Mutation" not in output

    def test_nonsense_included_when_configured(self, tmp_path):
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _NONSENSE_ROW],
            extra_args={CancerGenomesService.EXCLUDE_VARIANT_CLASSIFICATIONS: ""},
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        assert "Nonsense_Mutation" in output

    def test_lncrna_excluded_by_default(self, tmp_path):
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _LNCRNA_ROW],
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        # lncRNA excluded because default include_biotypes=protein_coding
        assert "NR_LNCRNA1" not in output

    def test_lncrna_3frame_when_include_biotypes_all(self, tmp_path):
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _LNCRNA_ROW],
            extra_args={CancerGenomesService.INCLUDE_BIOTYPES: "all"},
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        # lncRNA has no CDS=, so 3-frame translation produces _RF1/_RF2/_RF3 headers
        assert "_RF1" in output
        assert "_RF2" in output
        assert "_RF3" in output
        # All entries belong to the lncRNA transcript
        assert "NR_LNCRNA1" in output

    def test_include_biotypes_protein_coding_excludes_lncrna(self, tmp_path):
        """Explicit protein_coding filter: same as default but set explicitly."""
        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _MISSENSE_ROW, _LNCRNA_ROW],
            extra_args={CancerGenomesService.INCLUDE_BIOTYPES: "protein_coding"},
        )
        service.cbioportal_to_proteindb()

        output = output_path.read_text()
        assert "NM_PROT01" in output
        assert "NR_LNCRNA1" not in output


class TestParallelEquivalence:
    """workers=2 must produce the same protein set as workers=1 (sequential)."""

    def _run(self, tmp_path, workers):
        output_path = tmp_path / f"proteindb_w{workers}.fa"
        service = CancerGenomesService({}, {
            CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: str(
                _TESTDATA / "test_cbioportal_grch37_mutations.txt"),
            CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: str(
                _TESTDATA / "test_cbioportal_ncbi_grch37_transcripts.fa"),
            CancerGenomesService.CONFIG_OUTPUT_FILE: str(output_path),
            CancerGenomesService.CONFIG_GFF_FILE: str(
                _TESTDATA / "test_cbioportal_ncbi_grch37.gff"),
            CancerGenomesService.WORKERS: workers,
            CancerGenomesService.CBIO_BATCH_SIZE: 3,  # small batch to exercise multi-batch logic
        })
        service.cbioportal_to_proteindb()
        return {r.id: str(r.seq) for r in SeqIO.parse(str(output_path), "fasta")}

    def test_parallel_matches_sequential(self, tmp_path):
        seq = self._run(tmp_path, workers=1)
        par = self._run(tmp_path, workers=2)
        assert set(par.keys()) == set(seq.keys()), (
            f"Header mismatch: only_sequential={set(seq.keys()) - set(par.keys())} "
            f"only_parallel={set(par.keys()) - set(seq.keys())}"
        )
        for key in seq:
            assert par[key] == seq[key], f"Sequence mismatch for {key}"


class TestSplitByFilterColumn:
    """Tests for split_by_filter_column: per-group proteinDB file generation.

    Covers:
    - A per-group FASTA is created alongside the main output (group name
      is derived by stripping non-alpha characters from the clinical
      CANCER_TYPE value, e.g. "Breast Cancer" -> "BreastCancer").
    - The group file IDs match the reference snapshot.
    - The group file is a subset of the main output.
    - Samples absent from the clinical file are skipped from both outputs.
    """

    def _run_split(self, tmp_path):
        output_path = tmp_path / "proteindb.fa"
        service = CancerGenomesService({}, {
            CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: str(
                _TESTDATA / "test_cbioportal_grch37_mutations.txt"),
            CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: str(
                _TESTDATA / "test_cbioportal_ncbi_grch37_transcripts.fa"),
            CancerGenomesService.CONFIG_OUTPUT_FILE: str(output_path),
            CancerGenomesService.CONFIG_GFF_FILE: str(
                _TESTDATA / "test_cbioportal_ncbi_grch37.gff"),
            CancerGenomesService.CLINICAL_SAMPLE_FILE: str(
                _TESTDATA / "test_cbioportal_grch37_clinical_sample.txt"),
            CancerGenomesService.SPLIT_BY_FILTER_COLUMN: True,
        })
        service.cbioportal_to_proteindb()
        return output_path

    def test_split_creates_group_file(self, tmp_path):
        output_path = self._run_split(tmp_path)
        group_file = output_path.with_name("proteindb_BreastCancer.fa")
        assert group_file.exists(), "Expected per-group output file for BreastCancer"
        assert group_file.stat().st_size > 0

    def test_split_group_file_matches_reference(self, tmp_path):
        output_path = self._run_split(tmp_path)
        group_file = output_path.with_name("proteindb_BreastCancer.fa")

        expected = {r.id for r in SeqIO.parse(
            str(_TESTDATA / "test_cbioportal_grch37_proteindb_BreastCancer.fa"), "fasta")}
        actual = {r.id for r in SeqIO.parse(str(group_file), "fasta")}
        assert actual == expected

    def test_split_group_is_subset_of_main_output(self, tmp_path):
        output_path = self._run_split(tmp_path)
        group_file = output_path.with_name("proteindb_BreastCancer.fa")

        main_ids = {r.id for r in SeqIO.parse(str(output_path), "fasta")}
        group_ids = {r.id for r in SeqIO.parse(str(group_file), "fasta")}
        assert group_ids.issubset(main_ids)

    def test_split_unknown_sample_produces_empty_output(self, tmp_path):
        """Samples absent from the clinical file are skipped when split_by_filter_column=True."""
        clinical_path = tmp_path / "clinical.txt"
        clinical_path.write_text("SAMPLE_ID\tCANCER_TYPE\nSAMPLE_OTHER\tBreast Cancer\n")

        service, output_path = _make_service(
            tmp_path, _FASTA,
            [_MAF_HEADER, _MISSENSE_ROW],
            extra_args={
                CancerGenomesService.SPLIT_BY_FILTER_COLUMN: True,
                CancerGenomesService.CLINICAL_SAMPLE_FILE: str(clinical_path),
            },
        )
        service.cbioportal_to_proteindb()

        records = list(SeqIO.parse(str(output_path), "fasta"))
        assert records == [], "Unknown sample must produce no output when split_by_filter_column=True"


class TestRealDataHgvsVerification:
    """Verify pipeline output against HGVSp annotations using real TCGA BRCA data.

    Each test runs cbioportal_to_proteindb on real GRCh37/GRCh38 NCBI data and
    checks that every translated protein matches its annotated HGVSp:
    - missense: exact amino acid at the mutated position
    - frameshift: prefix up to mutation site is identical to reference translation
    """

    def _run(self, tmp_path, mutations, fasta, gff):
        output_path = tmp_path / "proteindb.fa"
        service = CancerGenomesService({}, {
            CancerGenomesService.CONFIG_CANCER_GENOMES_MUTATION_FILE: str(mutations),
            CancerGenomesService.CONFIG_COMPLETE_GENES_FILE: str(fasta),
            CancerGenomesService.CONFIG_OUTPUT_FILE: str(output_path),
            CancerGenomesService.CONFIG_GFF_FILE: str(gff),
        })
        service.cbioportal_to_proteindb()
        return {r.id: str(r.seq)
                for r in SeqIO.parse(str(output_path), "fasta")}

    def test_grch37_hgvsp_verification(self, tmp_path):
        prots = self._run(
            tmp_path,
            _TESTDATA / "test_cbioportal_grch37_mutations.txt",
            _TESTDATA / "test_cbioportal_ncbi_grch37_transcripts.fa",
            _TESTDATA / "test_cbioportal_ncbi_grch37.gff",
        )
        fa = _TESTDATA / "test_cbioportal_ncbi_grch37_transcripts.fa"

        # CDC20 p.H410Y — missense SNP, position 410 His→Tyr
        key = "cbiomut:NM_001255:CDC20:p.H410Y:Missense_Mutation"
        assert key in prots, "CDC20 p.H410Y entry missing"
        ref = _translate_ref_cds(fa, "NM_001255")
        mut = prots[key]
        assert ref[409] == "H", "CDC20 ref pos 410 should be His"
        assert mut[409] == "Y", "CDC20 mut pos 410 should be Tyr (p.H410Y)"
        assert ref[400:409] == mut[400:409], "CDC20 flanking residues 401-409 unchanged"

        # GOLGA2 p.A167V — missense SNP, position 167 Ala→Val
        key = "cbiomut:NM_004486:GOLGA2:p.A167V:Missense_Mutation"
        assert key in prots, "GOLGA2 p.A167V entry missing"
        ref = _translate_ref_cds(fa, "NM_004486")
        mut = prots[key]
        assert ref[166] == "A", "GOLGA2 ref pos 167 should be Ala"
        assert mut[166] == "V", "GOLGA2 mut pos 167 should be Val (p.A167V)"
        assert ref[157:166] == mut[157:166], "GOLGA2 flanking residues 158-166 unchanged"

        # EP300 p.M1339Ifs*2 — frameshift insertion, prefix 1-1338 identical
        key = "cbiomut:NM_001429:EP300:p.M1339Ifs*2:Frame_Shift_Ins"
        assert key in prots, "EP300 p.M1339Ifs*2 entry missing"
        ref = _translate_ref_cds(fa, "NM_001429")
        mut = prots[key]
        assert mut[:1338] == ref[:1338], "EP300 prefix 1-1338 should match reference"
        assert mut[1338:1342] != ref[1338:1342], "EP300 residues after frameshift should differ"

        # TBX3 p.Y265Lfs*12 — frameshift deletion, prefix 1-264 identical
        key = "cbiomut:NM_016569:TBX3:p.Y265Lfs*12:Frame_Shift_Del"
        assert key in prots, "TBX3 p.Y265Lfs*12 entry missing"
        ref = _translate_ref_cds(fa, "NM_016569")
        mut = prots[key]
        assert mut[:264] == ref[:264], "TBX3 prefix 1-264 should match reference"
        assert mut[264] != ref[264], "TBX3 position 265 should differ (Y265L frameshift)"

        # PLK2 p.R387Sfs*3 — frameshift deletion, prefix 1-386 identical
        key = "cbiomut:NM_006622:PLK2:p.R387Sfs*3:Frame_Shift_Del"
        assert key in prots, "PLK2 p.R387Sfs*3 entry missing"
        ref = _translate_ref_cds(fa, "NM_006622")
        mut = prots[key]
        assert mut[:386] == ref[:386], "PLK2 prefix 1-386 should match reference"
        assert mut[386] != ref[386], "PLK2 position 387 should differ (R387S frameshift)"

        # KHDRBS1 p.E143del — in-frame deletion: one AA shorter, prefix 1-142 intact, suffix shifts left by 1
        key = "cbiomut:NM_006559:KHDRBS1:p.E143del:In_Frame_Del"
        assert key in prots, "KHDRBS1 p.E143del entry missing"
        ref = _translate_ref_cds(fa, "NM_006559")
        mut = prots[key]
        assert len(mut) == len(ref) - 1, "KHDRBS1 p.E143del: protein should be 1 AA shorter"
        assert ref[142] == "E", "KHDRBS1 ref position 143 should be Glu"
        assert mut[:142] == ref[:142], "KHDRBS1 prefix 1-142 should match reference"
        assert mut[142:150] == ref[143:151], "KHDRBS1 suffix after deletion should shift left by 1"

        # FLNB p.S2107dup — in-frame insertion: one AA longer, C-terminal suffix preserved
        key = "cbiomut:NM_001457:FLNB:p.S2107dup:In_Frame_Ins"
        assert key in prots, "FLNB p.S2107dup entry missing"
        ref = _translate_ref_cds(fa, "NM_001457")
        mut = prots[key]
        assert len(mut) == len(ref) + 1, "FLNB p.S2107dup: protein should be 1 AA longer"
        assert ref[2106] == "S", "FLNB ref position 2107 should be Ser"
        assert mut[-50:] == ref[-50:], "FLNB C-terminal 50 AA should be unchanged after insertion"

        # HNRNPH3 p.*347Sext*32 — nonstop: stop codon at 347 replaced by AA; extension into 3'UTR not modelled
        key = "cbiomut:NM_012207:HNRNPH3:p.*347Sext*32:Nonstop_Mutation"
        assert key in prots, "HNRNPH3 p.*347Sext*32 entry missing"
        ref = _translate_ref_cds(fa, "NM_012207")
        mut = prots[key]
        assert ref[346] == "*", "HNRNPH3 ref position 347 should be stop (*)"
        assert mut[346] != "*", "HNRNPH3 mut position 347 should be a non-stop AA"
        assert mut[:346] == ref[:346], "HNRNPH3 prefix 1-346 should match reference"

    def test_grch38_hgvsp_verification(self, tmp_path):
        prots = self._run(
            tmp_path,
            _TESTDATA / "test_cbioportal_grch38_mutations.txt",
            _TESTDATA / "test_cbioportal_ncbi_grch38_transcripts.fa",
            _TESTDATA / "test_cbioportal_ncbi_grch38.gff",
        )
        fa = _TESTDATA / "test_cbioportal_ncbi_grch38_transcripts.fa"

        # EP300 p.M1339Ifs*2 — frameshift insertion, prefix 1-1338 identical
        key = "cbiomut:NM_001429:EP300:p.M1339Ifs*2:Frame_Shift_Ins"
        assert key in prots, "EP300 p.M1339Ifs*2 entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_001429")
        mut = prots[key]
        assert mut[:1338] == ref[:1338], "EP300 prefix 1-1338 should match reference (GRCh38)"
        assert mut[1338:1342] != ref[1338:1342], "EP300 residues after frameshift should differ (GRCh38)"

        # TBX3 p.Y265Lfs*17 — frameshift deletion, prefix 1-264 identical
        key = "cbiomut:NM_016569:TBX3:p.Y265Lfs*17:Frame_Shift_Del"
        assert key in prots, "TBX3 p.Y265Lfs*17 entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_016569")
        mut = prots[key]
        assert mut[:264] == ref[:264], "TBX3 prefix 1-264 should match reference (GRCh38)"
        assert mut[264] != ref[264], "TBX3 position 265 should differ (GRCh38)"

        # PLK2 p.R387Sfs*3 — frameshift deletion, prefix 1-386 identical
        key = "cbiomut:NM_006622:PLK2:p.R387Sfs*3:Frame_Shift_Del"
        assert key in prots, "PLK2 p.R387Sfs*3 entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_006622")
        mut = prots[key]
        assert mut[:386] == ref[:386], "PLK2 prefix 1-386 should match reference (GRCh38)"
        assert mut[386] != ref[386], "PLK2 position 387 should differ (GRCh38)"

        # KHDRBS1 p.E143del — in-frame deletion, prefix 1-142 intact, 1 AA shorter
        key = "cbiomut:NM_006559:KHDRBS1:p.E143del:In_Frame_Del"
        assert key in prots, "KHDRBS1 p.E143del entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_006559")
        mut = prots[key]
        assert len(mut) == len(ref) - 1, "KHDRBS1 p.E143del: protein should be 1 AA shorter (GRCh38)"
        assert ref[142] == "E", "KHDRBS1 ref position 143 should be Glu (GRCh38)"
        assert mut[:142] == ref[:142], "KHDRBS1 prefix 1-142 should match reference (GRCh38)"
        assert mut[142:150] == ref[143:151], "KHDRBS1 suffix after deletion should shift left by 1 (GRCh38)"

        # PIK3R1 p.K575_T576dup — in-frame insertion: 2 AAs longer (6-nt insertion), C-terminal preserved
        key = "cbiomut:NM_181523:PIK3R1:p.K575_T576dup:In_Frame_Ins"
        assert key in prots, "PIK3R1 p.K575_T576dup entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_181523")
        mut = prots[key]
        assert len(mut) == len(ref) + 2, "PIK3R1 p.K575_T576dup: protein should be 2 AAs longer (GRCh38)"
        assert mut[-50:] == ref[-50:], "PIK3R1 C-terminal 50 AA should be unchanged after insertion (GRCh38)"

        # TNK2 p.*1039Sext*15 — nonstop: stop at 1039 becomes non-stop AA; 3'UTR extension not modelled
        key = "cbiomut:NM_005781:TNK2:p.*1039Sext*15:Nonstop_Mutation"
        assert key in prots, "TNK2 p.*1039Sext*15 entry missing (GRCh38)"
        ref = _translate_ref_cds(fa, "NM_005781")
        mut = prots[key]
        assert ref[1038] == "*", "TNK2 ref position 1039 should be stop (*) (GRCh38)"
        assert mut[1038] != "*", "TNK2 mut position 1039 should be a non-stop AA (GRCh38)"
        assert mut[:1038] == ref[:1038], "TNK2 prefix 1-1038 should match reference (GRCh38)"
