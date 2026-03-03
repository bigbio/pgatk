"""ClinVar VCF-to-protein pipeline.

Processes ClinVar VCF records against NCBI RefSeq GTF annotations to produce
a FASTA file of variant protein sequences.  Uses BedTools interval overlap
(not VEP) to find which transcripts are affected by each variant.
"""
from __future__ import annotations

import logging
import sqlite3
import tempfile
from pathlib import Path
from typing import Optional

import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pybedtools import BedTool

from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
from pgatk.config.registry import load_config
from pgatk.toolbox.vcf_utils import (
    check_overlap,
    get_altseq,
    get_orfs_vcf,
    write_output,
)

logger = logging.getLogger(__name__)


class ClinVarService:
    """Main ClinVar VCF-to-protein service.

    Parameters
    ----------
    vcf_file : str
        Path to ClinVar VCF file (numeric chromosome names).
    gtf_file : str
        Path to NCBI RefSeq GTF file (NC_ chromosome names).
    fasta_file : str
        Path to transcript nucleotide FASTA with ``CDS=start-end`` headers.
    assembly_report : str
        Path to NCBI assembly report for chromosome name mapping.
    output_file : str, optional
        Output FASTA file path.  Defaults to config value.
    config_file : str, optional
        Path to custom YAML config.  Defaults to bundled ``clinvar_config.yaml``.
    """

    def __init__(
        self,
        vcf_file: str,
        gtf_file: str,
        fasta_file: str,
        assembly_report: str,
        output_file: Optional[str] = None,
        config_file: Optional[str] = None,
    ) -> None:
        self._vcf_file = vcf_file
        self._gtf_file = gtf_file
        self._fasta_file = fasta_file
        self._assembly_report = assembly_report

        cfg = load_config("clinvar", config_file)
        self._cfg = cfg.get("clinvar_translation", {})

        self._output_file = output_file or self._cfg.get(
            "proteindb_output_file", "clinvar-peptide-database.fa"
        )
        self._translation_table = self._cfg.get("translation_table", 1)
        self._mito_translation_table = self._cfg.get("mito_translation_table", 2)
        self._protein_prefix = self._cfg.get("protein_prefix", "clinvar_")
        self._report_ref_seq = self._cfg.get("report_ref_seq", False)
        self._num_orfs = self._cfg.get("num_orfs", 1)
        self._clnsig_exclude = self._cfg.get(
            "clinical_significance_exclude",
            ["Benign", "Likely_benign", "Benign/Likely_benign"],
        )
        self._include_consequences = self._cfg.get(
            "include_consequences", ["all"]
        )

    # ------------------------------------------------------------------
    # Static helper methods — ClinVar INFO field parsers
    # ------------------------------------------------------------------

    @staticmethod
    def _get_info_field(info_str: str, key: str) -> str:
        """Extract the value of *key* from a VCF INFO string.

        INFO fields are semicolon-separated ``KEY=VALUE`` pairs.
        Returns an empty string when the key is absent.
        """
        for part in info_str.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                if k == key:
                    return v
        return ""

    @staticmethod
    def passes_clnsig_filter(clnsig: str, exclude_list: list[str]) -> bool:
        """Return True when *clnsig* should **not** be filtered out.

        An empty or missing CLNSIG always passes.
        """
        if not clnsig:
            return True
        return clnsig not in exclude_list

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

    @staticmethod
    def parse_mc_consequence(mc_field: str) -> str:
        """Return the first molecular consequence from the ClinVar MC field.

        MC values look like ``SO:0001583|missense_variant``.
        Returns the consequence string after the pipe, or empty string.
        """
        if not mc_field:
            return ""
        first = mc_field.split(",")[0]
        if "|" in first:
            return first.split("|", 1)[1]
        return first

    @staticmethod
    def parse_mc_consequences(mc_field: str) -> list[str]:
        """Return all molecular consequences from the ClinVar MC field."""
        if not mc_field:
            return []
        results = []
        for entry in mc_field.split(","):
            if "|" in entry:
                results.append(entry.split("|", 1)[1])
            else:
                results.append(entry)
        return results

    @staticmethod
    def parse_geneinfo(geneinfo: str) -> tuple[str, str]:
        """Parse the ClinVar GENEINFO field.

        Format is ``SYMBOL:GENEID`` (or ``SYM1:ID1|SYM2:ID2`` for multi-gene).
        Returns ``(gene_symbol, gene_id)`` for the first gene.
        """
        if not geneinfo:
            return ("", "")
        first = geneinfo.split("|")[0]
        if ":" in first:
            parts = first.split(":", 1)
            return (parts[0], parts[1])
        return (first, "")

    # ------------------------------------------------------------------
    # GTF parsing (gffutils, same pattern as EnsemblDataService)
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_gtf(gtf_file: str) -> gffutils.FeatureDB:
        """Parse a GTF file into a gffutils FeatureDB.

        The database is stored alongside the GTF as ``<name>.db``.
        """
        db_file = str(Path(gtf_file).with_suffix(".db"))
        try:
            gffutils.create_db(
                gtf_file,
                db_file,
                merge_strategy="create_unique",
                keep_order=True,
                disable_infer_transcripts=True,
                disable_infer_genes=True,
                verbose=False,
                force=False,
            )
        except (ValueError, sqlite3.OperationalError):
            logger.debug("gffutils DB already exists: %s", db_file)

        return gffutils.FeatureDB(db_file)

    @staticmethod
    def _get_features(
        db: gffutils.FeatureDB,
        feature_id: str,
        feature_types: Optional[list[str]] = None,
    ) -> tuple:
        """Retrieve chromosome, strand and coding feature intervals for a transcript.

        Returns ``(chrom, strand, coding_features)`` where *coding_features*
        is a list of ``[start, end, type]`` triples.  Returns ``(None, None, None)``
        when the feature is not found.
        """
        if feature_types is None:
            feature_types = ["CDS"]
        try:
            feature = db[feature_id]
        except gffutils.exceptions.FeatureNotFoundError:
            try:
                feature = db[feature_id.split(".")[0]]
            except gffutils.exceptions.FeatureNotFoundError:
                logger.warning(
                    "Feature %s not found in GTF database.", feature_id
                )
                return None, None, None

        coding_features = []
        for f in db.children(feature, featuretype=feature_types, order_by="end"):
            coding_features.append([f.start, f.end, f.featuretype])
        return feature.chrom, feature.strand, coding_features

    # ------------------------------------------------------------------
    # VCF reading (pandas, same pattern as EnsemblDataService)
    # ------------------------------------------------------------------

    @staticmethod
    def _read_vcf(vcf_file: str) -> tuple[list, pd.DataFrame]:
        """Read a VCF file and return metadata lines and a DataFrame of records."""
        HEADERS = {
            "CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        }

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

        vcf_df = pd.DataFrame(data, columns=HEADERS)
        return metadata, vcf_df

    # ------------------------------------------------------------------
    # BedTools-based transcript overlap annotation
    # ------------------------------------------------------------------

    @staticmethod
    def _build_overlap_map(
        vcf_df: pd.DataFrame,
        gtf_file: str,
        chrom_mapper: ChromosomeMapper,
    ) -> dict[str, list[str]]:
        """Find transcripts overlapping each VCF variant via BedTools.

        Builds a BED from an already-loaded DataFrame so the VCF file does not
        need to be read a second time.

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

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def run(self) -> str:
        """Execute the ClinVar VCF-to-protein pipeline.

        Returns the path to the output FASTA file.
        """
        logger.info("Starting ClinVar pipeline")

        # 1. Load chromosome mapper
        chrom_mapper = ChromosomeMapper.from_assembly_report(self._assembly_report)

        # 2. Parse GTF
        db = self._parse_gtf(self._gtf_file)

        # 3. Load transcript FASTA
        transcripts_dict = SeqIO.index(
            self._fasta_file,
            "fasta",
            key_function=lambda h: h.split("|")[0].split(" ")[0],
        )
        # Build mapping without version for fallback lookup
        transcript_id_mapping = {
            k.split(".")[0]: k for k in transcripts_dict.keys()
        }

        # 4. Read VCF once into DataFrame
        _metadata, vcf_df = self._read_vcf(self._vcf_file)

        # 5. Find overlapping transcripts via BedTools (from DataFrame)
        overlap_map = self._build_overlap_map(vcf_df, self._gtf_file, chrom_mapper)
        logger.info("Found %d variants with transcript overlaps", len(overlap_map))

        stats = {
            "variants_processed": 0,
            "variants_filtered_clnsig": 0,
            "variants_filtered_mc": 0,
            "variants_no_overlap": 0,
            "variants_translated": 0,
        }

        with open(self._output_file, "w", encoding="utf-8") as prots_fn:
            for _, record in vcf_df.iterrows():
                stats["variants_processed"] += 1

                # --- Validate alleles ---
                ref = str(record.REF)
                if any(c not in "ACGT" for c in ref):
                    continue
                alts = []
                for a in str(record.ALT).split(","):
                    a = a.strip()
                    if a and all(c in "ACGT" for c in a):
                        alts.append(a)
                if not alts:
                    continue

                # --- CLNSIG filter ---
                info = str(record.INFO)
                clnsig = self._get_info_field(info, "CLNSIG")
                if not self.passes_clnsig_filter(clnsig, self._clnsig_exclude):
                    stats["variants_filtered_clnsig"] += 1
                    continue

                # --- MC consequence filter ---
                mc_field = self._get_info_field(info, "MC")
                if not self.passes_mc_filter(mc_field, self._include_consequences):
                    stats["variants_filtered_mc"] += 1
                    continue

                # --- Parse gene symbol and CLNSIG for description ---
                gene_symbol, _ = self.parse_geneinfo(
                    self._get_info_field(info, "GENEINFO")
                )
                desc_str = f"{clnsig}|{gene_symbol}" if gene_symbol else clnsig

                # --- Find overlapping transcripts ---
                chrom = str(record.CHROM)
                pos = int(record.POS)

                # Translation table (mito vs standard)
                trans_table = self._translation_table
                if chrom.upper() in ("M", "MT"):
                    trans_table = self._mito_translation_table

                for alt in alts:
                    variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                    transcript_ids = overlap_map.get(variant_key, [])
                    if not transcript_ids:
                        stats["variants_no_overlap"] += 1
                        continue

                    for transcript_id in transcript_ids:
                        # Resolve transcript in FASTA
                        tid = transcript_id
                        if tid not in transcripts_dict:
                            tid = transcript_id_mapping.get(
                                transcript_id.split(".")[0], transcript_id
                            )
                        try:
                            fasta_record = transcripts_dict[tid]
                        except KeyError:
                            logger.debug(
                                "Transcript %s not found in FASTA", transcript_id
                            )
                            continue

                        ref_seq = fasta_record.seq
                        desc = str(fasta_record.description)

                        # Determine CDS info and feature types
                        cds_info: list[int] = []
                        feature_types = ["exon"]
                        num_orfs = 3
                        if "CDS=" in desc:
                            try:
                                cds_str = [
                                    p
                                    for p in desc.split()
                                    if p.startswith("CDS=")
                                ][0]
                                cds_info = [
                                    int(x)
                                    for x in cds_str.split("=")[1].split("-")
                                ]
                                feature_types = ["CDS", "stop_codon"]
                                num_orfs = self._num_orfs
                            except (ValueError, IndexError):
                                logger.debug(
                                    "Could not extract CDS info from header: %s",
                                    desc,
                                )

                        # Get features from GTF
                        feat_chrom, strand, features_info = self._get_features(
                            db, tid, feature_types
                        )
                        if feat_chrom is None:
                            continue

                        # Check overlap at feature level
                        var_end = pos + len(ref) - 1
                        if not check_overlap(pos, var_end, features_info):
                            continue

                        # Apply variant
                        coding_ref_seq, coding_alt_seq = get_altseq(
                            ref_seq,
                            Seq(ref),
                            Seq(alt),
                            pos,
                            strand,
                            features_info,
                            cds_info,
                        )

                        if coding_alt_seq == "":
                            continue

                        # Translate
                        ref_orfs, alt_orfs = get_orfs_vcf(
                            coding_ref_seq,
                            coding_alt_seq,
                            trans_table,
                            num_orfs,
                        )

                        # Build sequence ID
                        record_id = ""
                        if record.ID and str(record.ID) != ".":
                            record_id = str(record.ID)

                        seq_id = "_".join(
                            [
                                self._protein_prefix + record_id,
                                ".".join([chrom, str(pos), ref, alt]),
                                tid,
                            ]
                        )

                        write_output(
                            seq_id=seq_id,
                            desc=desc_str,
                            seqs=alt_orfs,
                            prots_fn=prots_fn,
                            seqs_filter=ref_orfs,
                        )

                        stats["variants_translated"] += 1

                        if self._report_ref_seq:
                            write_output(
                                seq_id=tid,
                                desc="",
                                seqs=ref_orfs,
                                prots_fn=prots_fn,
                            )

        logger.info("ClinVar pipeline complete. Stats: %s", stats)
        return self._output_file


def _extract_transcript_id(attrs: str) -> str:
    """Extract transcript_id value from a GTF attributes string."""
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith("transcript_id"):
            # transcript_id "NM_000001.1"
            value = part.split(" ", 1)
            if len(value) > 1:
                return value[1].strip().strip('"')
    return ""
