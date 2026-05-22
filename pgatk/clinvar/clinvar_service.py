"""ClinVar VCF-to-protein pipeline.

Processes ClinVar VCF records against NCBI RefSeq GTF annotations to produce
a FASTA file of variant protein sequences.  Uses BedTools interval overlap
(not VEP) to find which transcripts are affected by each variant.
"""
from __future__ import annotations

import logging
import multiprocessing
import os
import re
import shutil
import sqlite3
import time
import tempfile
from pathlib import Path
from typing import Optional

import gffutils
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from pybedtools import BedTool

from pgatk.clinvar.chromosome_mapper import ChromosomeMapper
from pgatk.config.registry import load_config
from pgatk.toolbox.general import open_vcf
from pgatk.toolbox.vcf_utils import (
    check_overlap,
    get_altseq,
    get_orfs_vcf,
    write_output,
)

logger = logging.getLogger(__name__)


class _FeatureCache:
    """Per-run memoization of _get_features() results keyed by (tid, feature_types)."""

    def __init__(self, maxsize: int = 50000) -> None:
        self._cache: dict[tuple, tuple] = {}
        self._maxsize = maxsize

    def get(self, key: tuple) -> Optional[tuple]:
        return self._cache.get(key)

    def put(self, key: tuple, value: tuple) -> None:
        if len(self._cache) >= self._maxsize:
            for k in list(self._cache.keys())[: self._maxsize // 5]:
                del self._cache[k]
        self._cache[key] = value


_CLINVAR_BATCH_SIZE = 50_000

# Per-worker state populated once by _clinvar_worker_init and reused across all
# tasks assigned to that worker process.
_clinvar_worker_state: dict = {}


def _fasta_key_fn(header: str) -> str:
    """Key function for SeqIO FASTA indexing — strips rna- prefix and version."""
    return header.split("|")[0].split(" ")[0].removeprefix("rna-")


def _ensure_fasta_index(fasta_file: str) -> str:
    """Return path to a SQLite SeqIO index for fasta_file, building it if absent or stale."""
    idx_path = fasta_file + ".idx"
    if os.path.exists(idx_path):
        if os.path.getmtime(idx_path) >= os.path.getmtime(fasta_file):
            return idx_path
        try:
            os.remove(idx_path)
        except OSError:
            pass
    SeqIO.index_db(idx_path, [fasta_file], "fasta", key_function=_fasta_key_fn)
    return idx_path


def _split_vcf_into_batches(
    vcf_file: str, output_dir: str, batch_size: int = _CLINVAR_BATCH_SIZE
) -> list:
    """Stream vcf_file once, writing fixed-size variant-count batches into output_dir.

    Returns an ordered list of batch VCF paths.
    """
    header: list = []
    batch_paths: list = []
    handle = None
    count = 0
    try:
        with open_vcf(vcf_file) as f:
            for line in f:
                if line.startswith("#"):
                    header.append(line)
                    continue
                if not line.strip():
                    continue
                if handle is None or count >= batch_size:
                    if handle is not None:
                        handle.close()
                    chunk_path = os.path.join(
                        output_dir, f"batch_{len(batch_paths):04d}.vcf"
                    )
                    handle = open(chunk_path, "w", encoding="utf-8")
                    handle.writelines(header)
                    batch_paths.append(chunk_path)
                    count = 0
                handle.write(line)
                count += 1
    finally:
        if handle is not None:
            handle.close()
    return batch_paths


def _clinvar_worker_init(
    vcf_file: str,
    gff_file: str,
    fasta_file: str,
    assembly_report: str,
    output_file: str,
    config_file: Optional[str],
    overlap_map: dict,
) -> None:
    """Pool initializer: reconstruct service and open DB + FASTA index once per worker.

    The overlap_map (built in the main process) is pickled once per worker at
    pool startup via initargs, not once per task.
    """
    svc = ClinVarService(
        vcf_file=vcf_file,
        gff_file=gff_file,
        fasta_file=fasta_file,
        assembly_report=assembly_report,
        output_file=output_file,
        config_file=config_file,
    )
    db = ClinVarService._parse_gtf(gff_file)
    idx_path = _ensure_fasta_index(fasta_file)
    transcripts_dict = SeqIO.index_db(
        idx_path, [fasta_file], "fasta", key_function=_fasta_key_fn
    )
    transcript_id_mapping = {k.split(".")[0]: k for k in transcripts_dict.keys()}
    _clinvar_worker_state.update(
        {
            "svc": svc,
            "db": db,
            "transcripts_dict": transcripts_dict,
            "transcript_id_mapping": transcript_id_mapping,
            "overlap_map": overlap_map,
        }
    )


def _clinvar_worker(vcf_batch_path: str, output_path: str) -> dict:
    """Process one VCF batch chunk using per-worker pre-initialized state."""
    state = _clinvar_worker_state
    _, vcf_records = ClinVarService._read_vcf(vcf_batch_path)
    return state["svc"]._process_batch(
        vcf_records=vcf_records,
        overlap_map=state["overlap_map"],
        db=state["db"],
        transcripts_dict=state["transcripts_dict"],
        transcript_id_mapping=state["transcript_id_mapping"],
        output_file=output_path,
    )


class ClinVarService:
    """Main ClinVar VCF-to-protein service.

    Parameters
    ----------
    vcf_file : str
        Path to ClinVar VCF file (numeric chromosome names).
    gff_file : str
        Path to NCBI RefSeq GFF3 annotation file (NC_ chromosome names).
        Previously accepted a GTF; GFF3 is now required because the NCBI GTF
        lacks the parent hierarchy that gffread needs to generate CDS= headers.
    fasta_file : str
        Path to transcript nucleotide FASTA with ``CDS=start-end`` headers.
        Generate with ``ncbi-downloader --generate-transcripts``.
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
        gff_file: str,
        fasta_file: str,
        assembly_report: str,
        output_file: Optional[str] = None,
        config_file: Optional[str] = None,
    ) -> None:
        self._vcf_file = vcf_file
        self._gtf_file = gff_file
        self._fasta_file = fasta_file
        self._assembly_report = assembly_report
        self._config_file = config_file

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
        biotypes_raw = self._cfg.get("include_biotypes", "all")
        if isinstance(biotypes_raw, str):
            self._include_biotypes = [b.strip() for b in biotypes_raw.split(",")]
        else:
            self._include_biotypes = list(biotypes_raw)
        try:
            self._workers = max(1, int(self._cfg.get("workers", 1)))
        except (TypeError, ValueError):
            self._workers = 1

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

        An empty or missing CLNSIG always passes.  Compound values like
        ``Benign/Likely_benign`` or ``Benign,_risk_factor`` are split on
        ``/`` and ``,`` so that each component is checked individually.
        """
        if not clnsig:
            return True
        # Exact match first (covers simple values and compound entries
        # that are explicitly listed, e.g. "Benign/Likely_benign").
        if clnsig in exclude_list:
            return False
        # Split compound values and check each component.
        parts = re.split(r"[/,]", clnsig)
        non_empty = [p.strip() for p in parts if p.strip()]
        if not non_empty:
            return True  # delimiter-only or whitespace-only → treat as empty
        return not all(p.replace("_", " ").lower() in
                       [e.replace("_", " ").lower() for e in exclude_list]
                       for p in non_empty)

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
        if not os.path.exists(db_file):
            logger.warning(
                "Building gffutils database from %s — this takes 10–30 minutes "
                "on first run and is cached as %s for subsequent runs.",
                gtf_file,
                db_file,
            )
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
                try:
                    # GFF3: NCBI prefixes transcript IDs with "rna-"
                    feature = db[f"rna-{feature_id}"]
                except gffutils.exceptions.FeatureNotFoundError:
                    logger.warning(
                        "Feature %s not found in annotation database.", feature_id
                    )
                    return None, None, None

        coding_features = []
        for f in db.children(feature, featuretype=feature_types, order_by="end"):
            coding_features.append([f.start, f.end, f.featuretype])
        return feature.chrom, feature.strand, coding_features

    @staticmethod
    def _get_transcript_biotype(db: gffutils.FeatureDB, transcript_id: str) -> str:
        """Extract gene_biotype from a gffutils transcript/mRNA feature.

        Returns empty string if the transcript or attribute is not found.
        In GFF3, gene_biotype lives on the parent gene feature; this method
        traverses up when the transcript feature itself lacks the attribute.
        """
        feature = None
        try:
            feature = db[transcript_id]
        except gffutils.exceptions.FeatureNotFoundError:
            try:
                feature = db[transcript_id.split(".")[0]]
            except gffutils.exceptions.FeatureNotFoundError:
                try:
                    feature = db[f"rna-{transcript_id}"]  # GFF3 NCBI prefix
                except gffutils.exceptions.FeatureNotFoundError:
                    base = transcript_id.split(".")[0]
                    for ftype in ("transcript", "mRNA"):
                        for f in db.all_features(featuretype=ftype):
                            # Strip GFF3 "rna-" prefix before version comparison
                            if f.id.removeprefix("rna-").split(".")[0] == base:
                                feature = f
                                break
                        if feature is not None:
                            break
                    if feature is None:
                        return ""
        try:
            biotype = feature.attributes.get("gene_biotype", [""])[0]
            if not biotype:
                # GFF3: gene_biotype is on the parent gene feature, not the mRNA
                for parent in db.parents(feature, featuretype="gene"):
                    biotype = parent.attributes.get("gene_biotype", [""])[0]
                    if biotype:
                        break
            return biotype
        except (IndexError, AttributeError):
            return ""

    # ------------------------------------------------------------------
    # VCF reading (pandas, same pattern as EnsemblDataService)
    # ------------------------------------------------------------------

    _VCFRecord = collections.namedtuple(
        "VCFRecord", ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    )

    @staticmethod
    def _read_vcf(vcf_file: str) -> tuple[list, list]:
        """Read a VCF file and return metadata lines and a list of VCFRecord namedtuples."""
        metadata: list[str] = []
        records: list = []
        with open_vcf(vcf_file) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    metadata.append(line)
                else:
                    parts = line.split("\t", 7)
                    if len(parts) < 8:
                        continue
                    records.append(ClinVarService._VCFRecord(
                        CHROM=parts[0], POS=int(parts[1]), ID=parts[2],
                        REF=parts[3], ALT=parts[4], QUAL=parts[5],
                        FILTER=parts[6], INFO=parts[7],
                    ))
        return metadata, records

    # ------------------------------------------------------------------
    # BedTools-based transcript overlap annotation
    # ------------------------------------------------------------------

    @staticmethod
    def _build_overlap_map(
        vcf_records: list,
        gtf_file: str,
        chrom_mapper: ChromosomeMapper,
    ) -> dict[str, list[str]]:
        """Find transcripts overlapping each VCF variant via BedTools.

        Returns a dict mapping ``"CHROM:POS:REF:ALT"`` variant keys to lists
        of overlapping transcript IDs.
        """
        bed_lines: list[str] = []
        for row in vcf_records:
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
    # Batch processing (used by both sequential and parallel paths)
    # ------------------------------------------------------------------

    def _process_batch(
        self,
        vcf_records: list,
        overlap_map: dict,
        db: gffutils.FeatureDB,
        transcripts_dict,
        transcript_id_mapping: dict,
        output_file: str,
    ) -> dict:
        """Translate one DataFrame of VCF records, writing protein sequences to output_file.

        Returns a stats dict with per-category counts.
        """
        stats = {
            "variants_processed": 0,
            "variants_filtered_clnsig": 0,
            "variants_filtered_mc": 0,
            "variants_no_overlap": 0,
            "transcripts_filtered_biotype": 0,
            "variants_translated": 0,
        }

        processed_pairs: set = set()
        feature_cache = _FeatureCache()
        biotype_cache: dict = {}
        missing_transcripts: set = set()

        with open(output_file, "w", encoding="utf-8") as prots_fn:
            for record in vcf_records:
                stats["variants_processed"] += 1
                if stats["variants_processed"] % 1000 == 0:
                    logger.info(
                        "Progress: %d variants processed | %d filtered (CLNSIG) | "
                        "%d filtered (MC) | %d translated",
                        stats["variants_processed"],
                        stats["variants_filtered_clnsig"],
                        stats["variants_filtered_mc"],
                        stats["variants_translated"],
                    )

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

                info = str(record.INFO)
                clnsig = self._get_info_field(info, "CLNSIG")
                if not self.passes_clnsig_filter(clnsig, self._clnsig_exclude):
                    stats["variants_filtered_clnsig"] += 1
                    continue

                mc_field = self._get_info_field(info, "MC")
                if not self.passes_mc_filter(mc_field, self._include_consequences):
                    stats["variants_filtered_mc"] += 1
                    continue

                gene_symbol, _ = self.parse_geneinfo(
                    self._get_info_field(info, "GENEINFO")
                )
                sig_label = clnsig if clnsig else "not_provided"
                desc_str = f"{sig_label}|{gene_symbol}" if gene_symbol else sig_label

                chrom = str(record.CHROM)
                pos = int(record.POS)

                trans_table = self._translation_table
                chrom_bare = chrom.lstrip("chr").upper()
                if chrom_bare in ("M", "MT"):
                    trans_table = self._mito_translation_table

                for alt in alts:
                    variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                    transcript_ids = overlap_map.get(variant_key, [])
                    if not transcript_ids:
                        stats["variants_no_overlap"] += 1
                        continue

                    for transcript_id in transcript_ids:
                        pair_key = f"{variant_key}:{transcript_id}"
                        if pair_key in processed_pairs:
                            continue
                        processed_pairs.add(pair_key)

                        tid = transcript_id
                        if tid not in transcripts_dict:
                            tid = transcript_id_mapping.get(
                                transcript_id.split(".")[0], transcript_id
                            )

                        if tid in missing_transcripts:
                            continue

                        try:
                            fasta_record = transcripts_dict[tid]
                        except KeyError:
                            missing_transcripts.add(tid)
                            logger.debug(
                                "Transcript %s not found in FASTA", transcript_id
                            )
                            continue

                        if self._include_biotypes != ["all"]:
                            if tid not in biotype_cache:
                                biotype_cache[tid] = self._get_transcript_biotype(
                                    db, tid
                                )
                            biotype = biotype_cache[tid]
                            if biotype and biotype not in self._include_biotypes:
                                stats["transcripts_filtered_biotype"] += 1
                                continue

                        ref_seq = fasta_record.seq
                        desc = str(fasta_record.description)

                        cds_info: list = []
                        feature_types = ["exon"]
                        num_orfs = 3
                        if "CDS=" in desc:
                            try:
                                cds_str = [
                                    p.strip("[]")
                                    for p in desc.split()
                                    if p.strip("[]").startswith("CDS=")
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

                        cache_key = (tid, tuple(feature_types))
                        cached = feature_cache.get(cache_key)
                        if cached is None:
                            cached = self._get_features(db, tid, feature_types)
                            feature_cache.put(cache_key, cached)
                        feat_chrom, strand, features_info = cached
                        if feat_chrom is None:
                            continue

                        var_end = pos + len(ref) - 1
                        if not check_overlap(pos, var_end, features_info):
                            continue

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

                        ref_orfs, alt_orfs = get_orfs_vcf(
                            coding_ref_seq,
                            coding_alt_seq,
                            trans_table,
                            num_orfs,
                        )

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

        return stats

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def run(self, workers: Optional[int] = None) -> str:
        """Execute the ClinVar VCF-to-protein pipeline.

        Parameters
        ----------
        workers : int, optional
            Number of parallel worker processes.  Defaults to the ``workers``
            key in the config file (default 1).  Pass ``1`` to force
            sequential execution regardless of config.

        Returns the path to the output FASTA file.
        """
        if workers is None:
            workers = self._workers

        logger.info("Starting ClinVar pipeline (workers=%d)", workers)

        chrom_mapper = ChromosomeMapper.from_assembly_report(self._assembly_report)

        logger.info("Loading GFF annotation database from %s", self._gtf_file)
        t0 = time.perf_counter()
        db = self._parse_gtf(self._gtf_file)
        logger.info("GFF annotation database ready (%.1f s)", time.perf_counter() - t0)

        logger.info("Indexing transcript FASTA (%s)", self._fasta_file)
        t0 = time.perf_counter()
        if workers > 1:
            # Build a persistent SQLite index so all worker processes can open
            # the same file simultaneously without racing on an in-memory dict.
            idx_path = _ensure_fasta_index(self._fasta_file)
            transcripts_dict = SeqIO.index_db(
                idx_path, [self._fasta_file], "fasta", key_function=_fasta_key_fn
            )
        else:
            transcripts_dict = SeqIO.index(
                self._fasta_file, "fasta", key_function=_fasta_key_fn
            )
        transcript_id_mapping = {k.split(".")[0]: k for k in transcripts_dict.keys()}
        logger.info(
            "FASTA indexed: %d transcripts (%.1f s)",
            len(transcripts_dict),
            time.perf_counter() - t0,
        )

        _metadata, vcf_records = self._read_vcf(self._vcf_file)

        logger.info("Building BedTools transcript overlap map ...")
        t0 = time.perf_counter()
        overlap_map = self._build_overlap_map(vcf_records, self._gtf_file, chrom_mapper)
        logger.info(
            "Found %d variants with transcript overlaps (%.1f s)",
            len(overlap_map),
            time.perf_counter() - t0,
        )

        # ---- Sequential path ----
        if workers <= 1:
            stats = self._process_batch(
                vcf_records, overlap_map, db, transcripts_dict,
                transcript_id_mapping, self._output_file,
            )
            logger.info("ClinVar pipeline complete. Stats: %s", stats)
            return self._output_file

        # ---- Parallel path ----
        with tempfile.TemporaryDirectory(prefix="pgatk_clinvar_") as tmpdir:
            logger.info(
                "Splitting VCF into batches (batch_size=%d) ...", _CLINVAR_BATCH_SIZE
            )
            batch_paths = _split_vcf_into_batches(self._vcf_file, tmpdir)
            logger.info("Created %d batch(es)", len(batch_paths))

            if len(batch_paths) <= 1:
                # Entire VCF fits in one batch — skip pool overhead.
                stats = self._process_batch(
                    vcf_records, overlap_map, db, transcripts_dict,
                    transcript_id_mapping, self._output_file,
                )
                logger.info("ClinVar pipeline complete. Stats: %s", stats)
                return self._output_file

            n_workers = min(workers, len(batch_paths))
            tasks = [
                (bp, os.path.join(tmpdir, f"out_{i:04d}.fa"))
                for i, bp in enumerate(batch_paths)
            ]

            logger.info(
                "Dispatching %d batch(es) across %d worker(s)",
                len(batch_paths),
                n_workers,
            )

            # overlap_map is passed as an initarg — pickled once per worker
            # at pool startup, not once per task.
            with multiprocessing.get_context("spawn").Pool(
                n_workers,
                initializer=_clinvar_worker_init,
                initargs=(
                    self._vcf_file,
                    self._gtf_file,
                    self._fasta_file,
                    self._assembly_report,
                    self._output_file,
                    self._config_file,
                    overlap_map,
                ),
            ) as pool:
                all_stats = pool.starmap(_clinvar_worker, tasks)

            # Concatenate per-batch output FASTAs into the final file.
            with open(self._output_file, "wb") as out:
                for _, batch_out in tasks:
                    if os.path.exists(batch_out):
                        with open(batch_out, "rb") as f:
                            shutil.copyfileobj(f, out)

        # Aggregate per-batch stats and emit a single summary.
        combined: dict = {}
        for s in all_stats:
            for k, v in s.items():
                combined[k] = combined.get(k, 0) + v
        logger.info(
            "ClinVar pipeline complete (%d batches). Stats: %s",
            len(all_stats),
            combined,
        )
        return self._output_file


def _extract_transcript_id(attrs: str) -> str:
    """Extract transcript_id from a GTF attribute string or GFF3 Parent field.

    Handles three formats:
    - GTF:   ``transcript_id "NM_000001.1"``   (space-quoted)
    - GFF3:  ``transcript_id=NM_000001.1``      (key=value)
    - GFF3 fallback: ``Parent=rna-NM_000001.1`` (strips ``rna-`` prefix)
    """
    parent_val = ""
    for part in attrs.split(";"):
        part = part.strip()
        if not part:
            continue
        # GTF: transcript_id "NM_000001.1"
        if part.startswith("transcript_id "):
            return part.split(" ", 1)[1].strip().strip('"')
        # GFF3: transcript_id=NM_000001.1
        if part.startswith("transcript_id="):
            return part[len("transcript_id="):].strip()
        # Save GFF3 Parent for fallback
        if part.startswith("Parent=") and not parent_val:
            parent_val = part[len("Parent="):]
    # GFF3 CDS features reference their transcript via Parent=rna-NM_000001.1
    if parent_val:
        return parent_val.split(",")[0].removeprefix("rna-").strip()
    return ""
