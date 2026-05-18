from __future__ import annotations

import gzip
import logging
import os
import sqlite3
from pathlib import Path
from typing import Any, Optional

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from pybedtools import BedTool
import pandas as pd
from pgatk.toolbox.general import ParameterConfiguration
from pgatk.gnomad.data_downloader import check_gencode_version_compatibility
from pgatk.toolbox.vcf_utils import (
    check_overlap as _check_overlap,
    get_altseq as _get_altseq,
    get_orfs_vcf as _get_orfs_vcf,
    write_output as _write_output,
)


class _FeatureCache:

    """Per-run memoization of get_features() results, keyed by transcript id.

    Bounded by a soft maxsize to avoid unbounded growth on whole-genome runs.
    Defaults to 50000, which comfortably covers human ENSEMBL (~250k transcripts
    across all chromosomes; a single chr is ~10-20k).
    """

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


_MISSING = object()

# Number of VCF data lines per variant-batch chunk.  Tune up for very large
# VCFs (whole-genome) or down when memory per worker is constrained.
_VARIANT_BATCH_SIZE = 50_000

# Per-worker state populated once by _worker_init and reused across all tasks
# assigned to that worker process.
_worker_state: dict = {}


def _worker_init(default_params: dict, pipeline_args: dict,
                 input_fasta: str, gtf_path: str,
                 annotation_indices: Optional[tuple] = None) -> None:
    """Pool initializer: open the GTF database and FASTA index once per worker
    process rather than once per task, amortising startup cost across all
    variant-batch chunks handled by that worker.
    annotation_indices is pre-computed by the main process to avoid re-parsing
    the VCF ##INFO header for every batch.
    """
    svc = EnsemblDataService(default_params, pipeline_args)
    db_path = str(Path(gtf_path).with_suffix('.db'))
    db = svc.parse_gtf(gtf_path, db_path)
    idx_path = _ensure_fasta_index(input_fasta)
    transcripts_dict = SeqIO.index_db(
        idx_path, [input_fasta], "fasta", key_function=EnsemblDataService.get_key
    )
    _worker_state.update({'svc': svc, 'db': db, 'transcripts_dict': transcripts_dict,
                          'annotation_indices': annotation_indices})


def _ensure_fasta_index(input_fasta: str) -> str:
    """Return the path to a `.idx` SQLite index for input_fasta, building it
    if absent or stale. Stale means: the .idx exists but the FASTA's mtime
    is newer than the .idx's. The index is keyed by `EnsemblDataService.get_key`.
    """
    idx_path = input_fasta + ".idx"
    if os.path.exists(idx_path):
        if os.path.getmtime(idx_path) >= os.path.getmtime(input_fasta):
            return idx_path
        try:
            os.remove(idx_path)
        except OSError:
            pass
    # Build the index. SeqIO.index_db materialises the SQLite file on disk
    # as a side effect; we don't need to keep the returned handle here.
    SeqIO.index_db(idx_path, [input_fasta], "fasta", key_function=EnsemblDataService.get_key)
    return idx_path


def _split_vcf_into_batches(vcf_file: str, output_dir: str,
                             batch_size: int = _VARIANT_BATCH_SIZE) -> list[str]:
    """Stream `vcf_file` once, writing fixed-size variant-count batches into
    `output_dir`.

    Each output file is `<output_dir>/batch_<NNNN>.vcf` and contains the full
    VCF header followed by at most `batch_size` data lines.  Batches may span
    chromosome boundaries — the per-variant chromosome check inside
    `_vcf_to_proteindb_chunk` handles cross-chromosome records correctly.

    Returns an ordered list of chunk paths.  Constant memory — holds exactly
    one open file handle at a time.
    """
    header: list[str] = []
    batch_paths: list[str] = []
    handle = None
    count = 0
    try:
        with open(vcf_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    header.append(line)
                    continue
                if not line.strip():
                    continue
                if handle is None or count >= batch_size:
                    if handle is not None:
                        handle.close()
                    chunk_path = os.path.join(output_dir, f"batch_{len(batch_paths):04d}.vcf")
                    handle = open(chunk_path, 'w', encoding='utf-8')
                    handle.writelines(header)
                    batch_paths.append(chunk_path)
                    count = 0
                handle.write(line)
                count += 1
    finally:
        if handle is not None:
            handle.close()
    return batch_paths


def _vcf_to_proteindb_worker(vcf_file: str, output_path: str) -> None:
    """Module-level worker for multiprocessing.Pool.starmap.

    Uses the GTF DB, FASTA index, and EnsemblDataService instance opened once
    per worker process by _worker_init rather than re-constructing them per task.
    """
    svc = _worker_state['svc']
    _, stats = svc._vcf_to_proteindb_chunk(
        vcf_file, None, None, output_path,
        db=_worker_state['db'],
        transcripts_dict=_worker_state['transcripts_dict'],
        log_summary=False,
        annotation_indices=_worker_state.get('annotation_indices'),
    )
    return stats


class EnsemblDataService(ParameterConfiguration):
    CONFIG_KEY_VCF = "ensembl_translation"
    INPUT_FASTA = "input_fasta"
    TRANSLATION_TABLE = "translation_table"
    MITO_TRANSLATION_TABLE = "mito_translation_table"
    HEADER_VAR_PREFIX = "protein_prefix"
    REPORT_REFERENCE_SEQ = "report_ref_seq"
    PROTEIN_DB_OUTPUT = "proteindb_output_file"
    ANNOTATION_FIELD_NAME = "annotation_field_name"
    AF_FIELD = "af_field"
    AF_THRESHOLD = "af_threshold"
    TRANSCRIPT_STR = 'transcript_str'
    CONSEQUENCE_STR = "consequence_str"
    EXCLUDE_BIOTYPES = "exclude_biotypes"
    EXCLUDE_CONSEQUENCES = "exclude_consequences"
    SKIP_INCLUDING_ALL_CDS = "skip_including_all_cds"
    INCLUDE_BIOTYPES = "include_biotypes"
    INCLUDE_CONSEQUENCES = "include_consequences"
    BIOTYPE_STR = "biotype_str"
    TRANSCRIPT_DESCRIPTION_SEP = "transcript_description_sep"
    SKIP_INCLUDING_ALL_CDSS = "skip_including_all_CDSs"
    CONFIG_KEY_DATA = "ensembl_translation"
    NUM_ORFS = "num_orfs"
    NUM_ORFS_COMPLEMENT = "num_orfs_complement"
    EXPRESSION_STR = "expression_str"
    EXPRESSION_THRESH = "expression_thresh"
    IGNORE_FILTERS = "ignore_filters"
    ACCEPTED_FILTERS = "accepted_filters"
    WORKERS = "workers"

    def __init__(self, config_file: dict, pipeline_arguments: dict) -> None:
        """
    Init the class with the specific parameters.
    :param config_file configuration file
    :param pipeline_arguments pipelines arguments
    """
        super(EnsemblDataService, self).__init__(self.CONFIG_KEY_DATA, config_file,
                                                 pipeline_arguments)

        self._proteindb_output = self.get_config_value(self.PROTEIN_DB_OUTPUT, 'peptide-database.fa')
        self._translation_table = self.get_config_value(self.TRANSLATION_TABLE, 1)

        self._mito_translation_table = self.get_translation_properties(variable=self.MITO_TRANSLATION_TABLE,
                                                                       default_value=2)
        self._header_protein_prefix = self.get_translation_properties(variable=self.HEADER_VAR_PREFIX, default_value="var_")
        self._report_reference_seq = self.get_translation_properties(variable=self.REPORT_REFERENCE_SEQ,
                                                                     default_value=False)
        self._annotation_field_name = self.get_translation_properties(variable=self.ANNOTATION_FIELD_NAME,
                                                                      default_value='CSQ')
        self._transcript_str = self.get_translation_properties(variable=self.TRANSCRIPT_STR, default_value='FEATURE')
        self._consequence_str = self.get_translation_properties(variable=self.CONSEQUENCE_STR,
                                                                default_value='CONSEQUENCE')
        self._af_field = self.get_translation_properties(variable=self.AF_FIELD, default_value='')

        self._af_threshold = self.get_translation_properties(variable=self.AF_THRESHOLD, default_value=0.01)
        self._af_threshold = float(self._af_threshold)

        self._exclude_biotypes = self.get_multiple_options(
            self.get_translation_properties(variable=self.EXCLUDE_BIOTYPES, default_value=''))
        self._exclude_consequences = self.get_multiple_options(
            self.get_translation_properties(variable=self.EXCLUDE_CONSEQUENCES,
                                            default_value='downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant, regulatory_region_variant'))

        self._skip_including_all_cds = self.get_translation_properties(variable=self.SKIP_INCLUDING_ALL_CDS,
                                                                       default_value=False)
        self._include_biotypes = self.get_multiple_options(
            self.get_translation_properties(variable=self.INCLUDE_BIOTYPES,
                                            default_value='protein_coding,protein_coding_CDS_not_defined,protein_coding_LoF,nonsense_mediated_decay,non_stop_decay,translated_processed_pseudogene,IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,TEC'))

        self._include_consequences = self.get_multiple_options(
            self.get_translation_properties(variable=self.INCLUDE_CONSEQUENCES, default_value='all'))
        self._biotype_str = self.get_translation_properties(variable=self.BIOTYPE_STR, default_value='transcript_biotype')

        self._transcript_description_sep = self.get_translation_properties(variable=self.TRANSCRIPT_DESCRIPTION_SEP,
                                                                           default_value=';')
        self._num_orfs = self.get_translation_properties(variable=self.NUM_ORFS, default_value=3)
        self._num_orfs_complement = self.get_translation_properties(variable=self.NUM_ORFS_COMPLEMENT, default_value=0)
        self._expression_str = self.get_translation_properties(variable=self.EXPRESSION_STR, default_value="")
        self._expression_thresh = self.get_translation_properties(variable=self.EXPRESSION_THRESH, default_value=5.0)

        self._ignore_filters = self.get_translation_properties(variable=self.IGNORE_FILTERS, default_value=False)

        self._accepted_filters = self.get_multiple_options(
            self.get_translation_properties(variable=self.ACCEPTED_FILTERS, default_value='PASS'))

        raw_workers = self.get_translation_properties(variable=self.WORKERS, default_value=1)
        try:
            self._workers = max(1, int(raw_workers))
        except (TypeError, ValueError):
            self._workers = 1

    def get_translation_properties(self, variable: str, default_value: Any) -> Any:
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_DATA in self.get_default_parameters() and \
                self.CONFIG_KEY_VCF in self.get_default_parameters()[self.CONFIG_KEY_DATA] and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][variable]
        return value_return

    def three_frame_translation(self, input_fasta: str) -> None:
        """
    This function translate a transcriptome into a 3'frame translation protein sequence database
    :param input_fasta: fasta input file
    :return:
    """

        with open(input_fasta, 'r', encoding='utf-8') as input_handle:
            with open(self._proteindb_output, 'w', encoding='utf-8') as output_handle:
                for record in SeqIO.parse(input_handle, 'fasta'):
                    seq = record.seq
                    RF1 = seq.translate(table=self._translation_table)
                    RF2 = seq[1::].translate(table=self._translation_table)
                    RF3 = seq[2::].translate(table=self._translation_table)

                    if record.id == "":
                        self.get_logger().warning("skip entries without id %s", record.description)
                        continue
                    output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF1', RF1))
                    output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF2', RF2))
                    output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF3', RF3))

    @staticmethod
    def get_multiple_options(options_str: str) -> list[str]:
        """
        This method takes an String like option1, option2, ... and produce and array [option1, option2,... ]
        :param options_str:
        :return: Array
        """
        return list(map(lambda x: x.strip(), options_str.split(",")))

    @staticmethod
    def check_overlap(var_start: int, var_end: int, features_info: Optional[list] = None) -> bool:
        """Check if a variant overlaps any of the given features.  Delegates to :func:`pgatk.toolbox.vcf_utils.check_overlap`."""
        return _check_overlap(var_start, var_end, features_info)

    @staticmethod
    def get_altseq(ref_seq: str, ref_allele: str, var_allele: str, var_pos: int, strand: str, features_info: list, cds_info: Optional[list] = None) -> tuple:
        """Modify a reference sequence based on a variant allele.  Delegates to :func:`pgatk.toolbox.vcf_utils.get_altseq`."""
        return _get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info)

    @staticmethod
    def parse_gtf(gene_annotations_gtf: str, gtf_db_file: str) -> gffutils.FeatureDB:
        """
    Convert GTF file into a FeatureDB
    :param gene_annotations_gtf:
    :param gtf_db_file:
    :return:
    """
        try:
            gffutils.create_db(gene_annotations_gtf, gtf_db_file, merge_strategy="create_unique",
                               keep_order=True, disable_infer_transcripts=True, disable_infer_genes=True,
                               verbose=True,
                               force=False)
        except (ValueError, sqlite3.OperationalError) as e:  # already exists
            logging.getLogger(__name__).warning("Database already exists: %s %s", e, gtf_db_file)

        db = gffutils.FeatureDB(gtf_db_file)
        return db

    @staticmethod
    def get_features(db: gffutils.FeatureDB, feature_id: str, feature_types: Optional[list] = None) -> tuple:
        """
    Get chr, genomic positions, strand and biotype for feature_id
    also genomic positions for all its elements (exons/cds&start_codon)
    :param db:
    :param feature_id:
    :param feature_types:
    :return:
    """

        if feature_types is None:
            feature_types = ['exon']
        try:
            feature = db[feature_id]
        except gffutils.exceptions.FeatureNotFoundError:  # remove version number from the ID
            try:
                feature = db[feature_id.split('.')[0]]
            except gffutils.exceptions.FeatureNotFoundError:
                logging.getLogger(__name__).warning(
                    "Feature %s found in fasta file but not in gtf file. Check that the fasta file and the "
                    "gtf files match. A common issue is when the fasta file have chromosome patches but not the gtf",
                    feature_id)
                return None, None, None
        coding_features = []
        features = db.children(feature, featuretype=feature_types, order_by='end')
        for f in features:
            f_type = f.featuretype
            coding_features.append([f.start, f.end, f_type])
        return feature.chrom, feature.strand, coding_features

    @staticmethod
    def get_orfs_vcf(ref_seq: str, alt_seq: str, translation_table: int, num_orfs: int = 1) -> tuple[list, list]:
        """Translate coding ref/alt sequences into ORFs.  Delegates to :func:`pgatk.toolbox.vcf_utils.get_orfs_vcf`."""
        return _get_orfs_vcf(ref_seq, alt_seq, translation_table, num_orfs)

    @staticmethod
    def get_orfs_dna(ref_seq: str, translation_table: int, num_orfs: int, num_orfs_complement: int, to_stop: bool) -> list:
        """
    translate the coding_ref into ORFs
    """

        ref_orfs = []
        for n in range(0, num_orfs):
            ref_orfs.append(ref_seq[n::].translate(translation_table, to_stop=to_stop))

        rev_ref_seq = ref_seq.reverse_complement()
        for n in range(0, num_orfs_complement):
            ref_orfs.append(rev_ref_seq[n::].translate(translation_table, to_stop=to_stop))

        return ref_orfs

    def dnaseq_to_proteindb(self, input_fasta: str) -> str:
        """
    translates DNA sequences to protein sequences
    :param input_fasta: input fasta file
    :return:
    """

        seq_dict = SeqIO.index(input_fasta, "fasta")

        with open(self._proteindb_output, 'w', encoding='utf-8') as prots_fn:
            for record_id in seq_dict.keys():

                ref_seq = seq_dict[record_id].seq
                desc = str(seq_dict[record_id].description)

                key_values = {}
                sep = self._transcript_description_sep
                desc = desc.replace(' ', sep)
                for value in desc.split(sep):
                    if value.split('=')[0] == 'cds' or value.split(':')[0] == 'cds':
                        value = value.replace('cds', 'CDS')
                    if '=' in value:
                        key_values[value.split('=')[0]] = value.split('=')[1]
                    elif ':' in value:
                        key_values[value.split(':')[0]] = value.split(':')[1]
                    elif value.split('=')[
                        0] == 'CDS':  # when only it is specified to be a CDS, it means the whole sequence to be used
                        key_values[value.split('=')[0]] = '{}-{}'.format(1, len(ref_seq))

                feature_biotype = ""
                if self._biotype_str:
                    try:
                        feature_biotype = key_values[self._biotype_str]
                    except KeyError:
                        msg = "Biotype info was not found in the header using {} for record {} {}".format(
                            self._biotype_str,
                            record_id, desc)
                        self.get_logger().debug(msg)

                if 'CDS' in key_values.keys() and (
                        not self._skip_including_all_cds or 'altORFs' in self._include_biotypes):
                    pass
                elif self._biotype_str and (feature_biotype == "" or (feature_biotype in self._exclude_biotypes or
                                                                      (feature_biotype not in self._include_biotypes and self._include_biotypes != [
                                                                          'all']))):
                    continue

                if self._expression_str:
                    try:
                        if float(key_values[self._expression_str]) < self._expression_thresh:
                            continue
                    except KeyError:
                        msg = "Expression information not found in the fasta header with expression_str: {} for record {} {}".format(
                            self._expression_str, record_id, desc)
                        self.get_logger().debug(msg)
                        continue
                    except TypeError:
                        msg = "Expression value is not of valid type (float) at record: {} {}".format(record_id,
                                                                                                      key_values[
                                                                                                          self._expression_str])
                        self.get_logger().debug(msg)
                        continue

                # translate the whole sequences (3 ORFs) for non CDS sequences and not take alt_ORFs for CDSs
                if 'CDS' not in key_values.keys() or ('CDS' in key_values.keys() and
                                                      ('altORFs' in self._include_biotypes or
                                                       self._include_biotypes == ['all'])):
                    ref_orfs = self.get_orfs_dna(ref_seq, self._translation_table, self._num_orfs,
                                                 self._num_orfs_complement, to_stop=False)
                    self.write_output(seq_id=self._header_protein_prefix + record_id, desc=desc, seqs=ref_orfs,
                                      prots_fn=prots_fn)

                # also allow for direct translation of the CDS, when the cds info exists in the fasta header skip_including_all_cds is false
                if 'CDS' in key_values.keys() and not self._skip_including_all_cds:
                    try:
                        cds_info = [int(x) for x in key_values['CDS'].split('-')]
                        ref_seq = ref_seq[cds_info[0] - 1:cds_info[1]]
                        ref_orfs = self.get_orfs_dna(ref_seq, self._translation_table, 1, 0, to_stop=True)
                        self.write_output(seq_id=record_id, desc=desc, seqs=ref_orfs, prots_fn=prots_fn)
                    except (ValueError, IndexError, KeyError):
                        self.get_logger().warning("Could not extract cds position from fasta header for: %s %s", record_id, desc)

        return self._proteindb_output

    @staticmethod
    def get_key(fasta_header: str) -> str:
        return fasta_header.split('|')[0].split(' ')[0]

    @staticmethod
    def annoate_vcf(vcf_file: str, gtf_file: str,
                    vcf_info_field_index: int = 7, record_type_index: int = 2,
                    gene_info_index: int = 8, gene_info_sep: str = ';',
                    transcript_str: str = 'transcript_id', transcript_info_sep: str = ' ',
                    annotation_str: str = 'transcriptOverlaps') -> str:
        """
    intersect vcf and a gtf, add ID of the overlapping transcript to the vcf INFO field
    """
        vcf_stem = str(Path.cwd() / Path(vcf_file).stem)

        BedTool(gtf_file).intersect(BedTool(vcf_file), wo=True).saveas(f"{vcf_stem}_all.bed")

        muts_dict: dict[str, list[str]] = {}
        # VCF columns start after the GTF columns (gene_info_index + 1)
        vcf_start_col = gene_info_index + 1
        with open(f"{vcf_stem}_all.bed", 'r', encoding='utf-8') as an:
            for line in an.readlines():
                sl = line.strip().split('\t')
                if sl[record_type_index].strip() != 'CDS':
                    continue
                transcript_id = 'NO_OVERLAP'

                gene_info = sl[gene_info_index].split(gene_info_sep)
                for info in gene_info:  # extract transcript id from the gtf info column
                    if info.strip().startswith(transcript_str):
                        transcript_id = info.strip().split(transcript_info_sep)[1].strip('"')
                        continue

                if transcript_id == 'NO_OVERLAP':
                    continue

                # Use canonical positional key (CHROM:POS:REF:ALT) for reliable
                # matching between BED output and raw VCF lines.
                try:
                    vcf_chrom = sl[vcf_start_col]
                    vcf_pos = sl[vcf_start_col + 1]
                    vcf_ref = sl[vcf_start_col + 3]
                    vcf_alt = sl[vcf_start_col + 4]
                    key = f"{vcf_chrom}\t{vcf_pos}\t{vcf_ref}\t{vcf_alt}"
                except IndexError:
                    continue
                muts_dict.setdefault(key, []).append(transcript_id)

        with open(f"{vcf_stem}_annotated.vcf", 'w', encoding='utf-8') as ann, open(vcf_file, 'r', encoding='utf-8') as v:
            for line in v.readlines():
                if line.startswith('#'):
                    ann.write(line)
                else:
                    sl = line.strip().split('\t')
                    if len(sl) < 8:
                        ann.write(line)
                        continue
                    key = f"{sl[0]}\t{sl[1]}\t{sl[3]}\t{sl[4]}"
                    if key in muts_dict:
                        sl[vcf_info_field_index] = '{};{}={}'.format(
                            sl[vcf_info_field_index].strip(), annotation_str,
                            ','.join(set(muts_dict[key])))
                        ann.write('\t'.join(sl) + '\n')
                    else:
                        ann.write(line)

        return f"{vcf_stem}_annotated.vcf"

    @staticmethod
    def vcf_from_file(vcf_file: str) -> tuple[list, pd.DataFrame]:
        """
    Read a VCF file and return a dataframe for the records
    as well as a list for the metadata
    """

        HEADERS = {
            'CHROM': str,
            'POS': int,
            'ID': str,
            'REF': str,
            'ALT': str,
            'QUAL': str,
            'FILTER': str,
            'INFO': str,
        }

        metadata = []
        data = []
        _compressed = vcf_file.endswith('.gz') or vcf_file.endswith('.bgz')
        _opener = gzip.open(vcf_file, 'rt', encoding='utf-8') if _compressed else open(vcf_file, 'r', encoding='utf-8')
        with _opener as vcf:
            line = vcf.readline().strip()
            while line:
                if line.startswith('#'):
                    metadata.append(line)
                else:
                    data.append(line.split('\t')[0:8])
                line = vcf.readline().strip()

        vcf_df = pd.DataFrame(data, columns=HEADERS)

        return metadata, vcf_df

    def _parse_annotation_indices(self, metadata: list, vcf_file: str = "") -> tuple:
        """Resolve transcript/consequence/biotype column positions from VCF ##INFO metadata.

        Takes already-read metadata lines (from vcf_from_file) so the file is not
        re-opened.  Returns (transcript_index, consequence_index, biotype_index); any
        value may be None if the corresponding field is absent from the CSQ FORMAT.
        transcript_index falls back to 0 when no structured FORMAT declaration exists.
        """
        annotation_cols: list[str] = []
        try:
            raw_cols = \
                [x for x in metadata if x.startswith('##INFO=<ID={}'.format(self._annotation_field_name))][
                    0].upper().split('FORMAT')[1].strip(' ').split(':')[-1].split('=')[-1].strip(' ').split('|')
            annotation_cols = [c.strip(' ">\n\r') for c in raw_cols]
        except IndexError:
            pass

        ti = ci = bi = None
        try:
            ti = annotation_cols.index(self._transcript_str.upper())
        except ValueError:
            self.get_logger().debug("Unable to find %s in metadata header %s of VCF file: %s",
                                    self._transcript_str, annotation_cols, vcf_file)
        try:
            ci = annotation_cols.index(self._consequence_str.upper())
        except ValueError:
            self.get_logger().debug("Unable to find %s in metadata header %s of VCF file: %s",
                                    self._consequence_str, annotation_cols, vcf_file)
        try:
            bi = annotation_cols.index(self._biotype_str.upper())
        except ValueError:
            self.get_logger().debug("Unable to find %s in metadata header %s of VCF file: %s",
                                    self._biotype_str, annotation_cols, vcf_file)
        if ti is None:
            ti = 0
        return ti, ci, bi

    def _vcf_to_proteindb_chunk(self, vcf_file: str, input_fasta: Optional[str],
                                gene_annotations_gtf: Optional[str], output_path: str,
                                *, db=None, transcripts_dict=None,
                                log_summary: bool = True,
                                annotation_indices: Optional[tuple] = None) -> tuple[str, dict]:
        """
    Generate proteins for variants by modifying sequences of affected transcripts.
    In case of already annotated variants it only considers variants within
    potential coding regions of the transcript (CDSs & stop codons for protein-coding genes,
    exons for non-protein coding genes)
    In case of not annotated variants, it considers all variants overlapping
    transcripts from the selected biotypes.
    :param vcf_file:
    :param input_fasta: path to the FASTA file (required when db/transcripts_dict are None)
    :param gene_annotations_gtf: path to the GTF file (required when db is None)
    :param output_path: path for writing the output FASTA
    :param db: pre-opened gffutils.FeatureDB (supplied by pool initializer in parallel runs)
    :param transcripts_dict: pre-opened SeqIO index (supplied by pool initializer in parallel runs)
    :return:
    """
        if db is None:
            db = self.parse_gtf(gene_annotations_gtf, str(Path(gene_annotations_gtf).with_suffix('.db')))
        if transcripts_dict is None:
            idx_path = _ensure_fasta_index(input_fasta)
            transcripts_dict = SeqIO.index_db(idx_path, [input_fasta], "fasta",
                                              key_function=self.get_key)
        # handle cases where the transcript has version in the GTF but not in the VCF
        # Built lazily on the first KeyError to avoid iterating 207k keys up-front.
        transcript_id_mapping: Optional[dict[str, str]] = None
        feature_cache = _FeatureCache()
        # Value is (ref_seq, desc) for a known transcript, or None for a transcript
        # we've already looked up and confirmed isn't in the FASTA index (avoids re-trying
        # the disk seek). _MISSING sentinel below distinguishes "not yet looked up".
        seq_cache: dict[str, Optional[tuple]] = {}

        transcript_index, consequence_index, biotype_index = None, None, None
        if self._annotation_field_name:
            if annotation_indices is not None:
                # Pre-computed by the parallel orchestrator — skip per-batch header re-parsing.
                transcript_index, consequence_index, biotype_index = annotation_indices
                _, vcf_reader = self.vcf_from_file(vcf_file)
            else:
                metadata, vcf_reader = self.vcf_from_file(vcf_file)
                transcript_index, consequence_index, biotype_index = \
                    self._parse_annotation_indices(metadata, vcf_file)

        else:
            # in case the given VCF is not annotated, annotate it by identifying the overlapping transcripts
            vcf_file = self.annoate_vcf(vcf_file, gene_annotations_gtf)
            metadata, vcf_reader = self.vcf_from_file(vcf_file)
            self._annotation_field_name = 'transcriptOverlaps'
            transcript_index = 0

        invalid_records = {'# variants with invalid record': 0,
                           '# variants not passing Filter': 0,
                           '# variants not passing AF threshold': 0,
                           '# feature IDs from VCF that are not found in the given FASTA file': 0,
                           '# variants successfully translated': 0}

        self._accepted_filters = [x.upper() for x in self._accepted_filters]

        with open(output_path, 'w', buffering=1 << 20, encoding='utf-8') as prots_fn:
            for record in vcf_reader.itertuples(index=False, name='VCFRecord'):
                trans = False
                if [x for x in str(record.REF) if x not in 'ACGT']:
                    invalid_records['# variants with invalid record'] += 1
                    if self.get_logger().isEnabledFor(logging.DEBUG):
                        self.get_logger().debug("Invalid VCF record, skipping: %s", record)
                    continue

                alts = []
                for alt in record.ALT.split(','):
                    if alt is None:
                        continue
                    elif [x for x in str(alt) if x not in 'ACGT']:
                        continue
                    alts.append(alt)
                if not alts:
                    invalid_records['# variants with invalid record'] += 1
                    if self.get_logger().isEnabledFor(logging.DEBUG):
                        self.get_logger().debug("Invalid VCF record, skipping: %s", record)
                    continue

                # Parse INFO once; avoids repeated split-and-search list-comprehensions per variant.
                info_kv: dict[str, str] = {}
                for entry in record.INFO.split(';'):
                    k, _, v = entry.partition('=')
                    info_kv[k] = v

                if not self._ignore_filters and self._accepted_filters != ['ALL']:
                    if record.FILTER and record.FILTER != '.' and record.FILTER != 'NA' and record.FILTER != '':  # None and empty means PASS
                        filters = set(record.FILTER.upper().split(','))
                        if ';' in record.FILTER and len(filters) <= 1:
                            filters = set(record.FILTER.upper().split(';'))
                        if not filters <= set(self._accepted_filters):
                            invalid_records['# variants not passing Filter'] += 1
                            continue
                if self._af_field:
                    try:
                        af = float(info_kv[self._af_field])
                    except (ValueError, KeyError):
                        invalid_records['# variants with invalid record'] += 1
                        continue

                    if af < self._af_threshold:
                        invalid_records['# variants not passing AF threshold'] += 1
                        continue

                trans_table = self._translation_table
                if str(record.CHROM).lstrip('chr').upper() in ['M', 'MT']:
                    trans_table = self._mito_translation_table

                processed_transcript_allele = set()
                try:
                    transcript_records = info_kv[self._annotation_field_name]
                except KeyError:
                    invalid_records['# variants with invalid record'] += 1
                    if self.get_logger().isEnabledFor(logging.DEBUG):
                        self.get_logger().debug(
                            "skipped record %s, no annotation feature was found", record)
                    continue

                for transcript_record in transcript_records.split(','):
                    transcript_info = transcript_record.split('|')
                    if consequence_index is not None:
                        try:
                            consequence = transcript_info[consequence_index]
                        except IndexError:
                            invalid_records['# variants with invalid record'] += 1
                            if self.get_logger().isEnabledFor(logging.DEBUG):
                                self.get_logger().debug(
                                    "Give a valid index for the consequence in the INFO field for: %s",
                                    transcript_record)
                            continue
                        except TypeError:
                            pass
                    if biotype_index is not None:
                        try:
                            biotype = transcript_info[biotype_index]
                        except IndexError:
                            invalid_records['# variants with invalid record'] += 1
                            if self.get_logger().isEnabledFor(logging.DEBUG):
                                self.get_logger().debug(
                                    "Give a valid index for the biotype in the INFO field for: %s",
                                    transcript_record)
                            continue
                        except TypeError:
                            pass

                    try:
                        transcript_id = transcript_info[transcript_index]
                    except IndexError:
                        invalid_records['# variants with invalid record'] += 1
                        if self.get_logger().isEnabledFor(logging.DEBUG):
                            self.get_logger().debug(
                                "Give a valid index for the Transcript IDs in the INFO field for: %s",
                                transcript_record)
                        continue
                    if transcript_id == "":
                        continue

                    try:
                        transcript_id_v = transcript_id_mapping[transcript_id]  # type: ignore[index]
                    except (KeyError, TypeError):
                        if transcript_id_mapping is None:
                            transcript_id_mapping = {k.split('.')[0]: k for k in transcripts_dict.keys()}
                            try:
                                transcript_id_v = transcript_id_mapping[transcript_id]
                            except KeyError:
                                transcript_id_v = transcript_id
                        else:
                            transcript_id_v = transcript_id

                    # Apply consequence/biotype filters before any I/O (FASTA lookup,
                    # get_features DB query).  Both values come from transcript_info
                    # which is already parsed above.
                    if consequence_index is not None:
                        if (consequence in self._exclude_consequences or
                                (consequence not in self._include_consequences and
                                 self._include_consequences != ['all'])):
                            continue

                    if biotype_index is not None:
                        if (biotype in self._exclude_biotypes or
                                (biotype not in self._include_biotypes and
                                 self._include_biotypes != ['all'])):
                            continue

                    cached_row = seq_cache.get(transcript_id_v, _MISSING)
                    if cached_row is _MISSING:
                        try:
                            row = transcripts_dict[transcript_id_v]
                            ref_seq = row.seq
                            desc = str(row.description)
                            seq_cache[transcript_id_v] = (ref_seq, desc)
                        except KeyError:
                            invalid_records['# feature IDs from VCF that are not found in the given FASTA file'] += 1
                            if self.get_logger().isEnabledFor(logging.DEBUG):
                                self.get_logger().debug(
                                    "Feature %s not found in fasta of the GTF file %s",
                                    transcript_id_v, record)
                            seq_cache[transcript_id_v] = None
                            continue
                    elif cached_row is None:
                        # already-known-missing; skip without re-counting in the stats
                        continue
                    else:
                        ref_seq, desc = cached_row

                    feature_types = ['exon']
                    cds_info = []
                    num_orfs = 3
                    if 'CDS=' in desc:
                        try:
                            cds_token = next(
                                (t.strip('[]') for t in desc.split(' ')
                                 if t.strip('[]').startswith('CDS=')),
                                None,
                            )
                            if cds_token:
                                cds_info = [int(x) for x in cds_token.split('=')[1].split('-')]
                                feature_types = ['CDS', 'stop_codon']
                                num_orfs = 1
                        except (ValueError, IndexError):
                            if self.get_logger().isEnabledFor(logging.DEBUG):
                                self.get_logger().debug(
                                    "Could not extract cds position from fasta header for: %s", desc)

                    cache_key = (transcript_id_v, tuple(feature_types))
                    cached = feature_cache.get(cache_key)
                    if cached is None:
                        chrom, strand, features_info = self.get_features(
                            db, transcript_id_v, feature_types)
                        feature_cache.put(cache_key, (chrom, strand, features_info))
                    else:
                        chrom, strand, features_info = cached
                    if chrom is None:
                        continue

                    for alt in alts:
                        dedup_key = transcript_id + str(record.REF) + str(alt)
                        if dedup_key in processed_transcript_allele:
                            continue
                        processed_transcript_allele.add(dedup_key)

                        # for non-CDSs, only consider the exon that actually overlaps the variant
                        try:
                            overlap_flag = self.check_overlap(int(record.POS),
                                                              int(record.POS) + len(str(record.REF)) - 1,
                                                              features_info)
                        except TypeError:
                            invalid_records['# variants with invalid record'] += 1
                            if self.get_logger().isEnabledFor(logging.DEBUG):
                                self.get_logger().debug("Wrong VCF record in %s", record)
                            continue

                        if (chrom.lstrip("chr") == str(record.CHROM).lstrip("chr") and
                                overlap_flag):
                            coding_ref_seq, coding_alt_seq = self.get_altseq(ref_seq, Seq(str(record.REF)),
                                                                             Seq(str(alt)), int(record.POS), strand,
                                                                             features_info, cds_info)
                            if coding_alt_seq != "":
                                ref_orfs, alt_orfs = self.get_orfs_vcf(coding_ref_seq, coding_alt_seq, trans_table,
                                                                       num_orfs)
                                record_id = ""
                                if record.ID:
                                    record_id = '_' + str(record.ID)
                                self.write_output(seq_id='_'.join([self._header_protein_prefix + str(record_id),
                                                                   '.'.join([str(record.CHROM), str(record.POS),
                                                                             str(record.REF), str(alt)]),
                                                                   transcript_id_v]),
                                                  desc='',
                                                  seqs=alt_orfs,
                                                  prots_fn=prots_fn,
                                                  seqs_filter=ref_orfs)
                                trans = True

                                if self._report_reference_seq:
                                    self.write_output(seq_id=transcript_id_v,
                                                      desc='',
                                                      seqs=ref_orfs,
                                                      prots_fn=prots_fn)
                if trans:
                    invalid_records['# variants successfully translated'] += 1

        if log_summary:
            msg = "Translation summary:\n {}".format(
                '\n'.join([x + ":" + str(invalid_records[x]) for x in invalid_records.keys()]))
            self.get_logger().info(msg)

        return output_path, invalid_records

    def vcf_to_proteindb(self, vcf_file: str, input_fasta: str, gene_annotations_gtf: str, workers=None) -> str:
        """Generate proteins for variants by modifying sequences of affected transcripts.

        If workers is None, falls back to self._workers (config) which defaults
        to 1 (sequential, backward-compatible). Pass workers > 1 to split the VCF
        into fixed-size variant batches and fan out via multiprocessing.Pool.
        :param vcf_file:
        :param input_fasta:
        :param gene_annotations_gtf:
        :param workers: number of parallel worker processes (None => use config default)
        :return: path to the output proteindb FASTA
        """
        if workers is None:
            workers = self._workers if self._workers else 1

        check_gencode_version_compatibility(vcf_file, gene_annotations_gtf)

        # Fast path: sequential — single call, identical behaviour to the original implementation.
        # For sequential runs we do NOT pre-annotate here; _vcf_to_proteindb_chunk handles that
        # in its else-branch so that transcript_index=0 is set correctly for unannotated VCFs.
        if workers <= 1:
            output_path, _ = self._vcf_to_proteindb_chunk(vcf_file, input_fasta, gene_annotations_gtf, self._proteindb_output)
            return output_path

        # Parallel path: pre-annotate unannotated VCFs in the main process. This avoids each
        # worker racing on the same bedtools-output bed file (which annoate_vcf writes to cwd)
        # and amortises the bedtools intersect across workers.
        if not self._annotation_field_name:
            vcf_file = self.annoate_vcf(vcf_file, gene_annotations_gtf)
            self._annotation_field_name = 'transcriptOverlaps'

        # Build the FASTA index and GTF database ONCE in the main process so all
        # workers find them ready.  Without this, N workers spawned simultaneously
        # would race to create the same .db file via gffutils.create_db, with losers
        # catching a sqlite3.OperationalError on a partially-written file.
        _ensure_fasta_index(input_fasta)
        self.parse_gtf(gene_annotations_gtf, str(Path(gene_annotations_gtf).with_suffix('.db')))

        # Parallel — split into fixed-size variant-batch VCFs, fan out to a Pool.
        import multiprocessing as mp
        import shutil
        import tempfile

        # Pipeline args forwarded to the pool initializer (annotation field already set above).
        pa = dict(self.get_pipeline_parameters())
        pa[EnsemblDataService.ANNOTATION_FIELD_NAME] = self._annotation_field_name

        # Parse annotation column indices ONCE from the (now-annotated) VCF so workers
        # skip re-parsing the identical header for every batch.
        _ann_meta, _ = self.vcf_from_file(vcf_file)
        _annotation_indices = self._parse_annotation_indices(_ann_meta, vcf_file)

        with tempfile.TemporaryDirectory(prefix='pgatk_v2p_') as tmpdir:
            # Stream-split VCF into fixed-size variant-count batches (constant memory).
            # Batches may span chromosome boundaries; _vcf_to_proteindb_chunk handles this.
            batch_paths = _split_vcf_into_batches(vcf_file, tmpdir)

            if len(batch_paths) <= 1:
                # Entire VCF fits in one batch — run sequentially without pool overhead.
                output_path, _ = self._vcf_to_proteindb_chunk(vcf_file, input_fasta, gene_annotations_gtf,
                                                               self._proteindb_output)
                return output_path

            tasks = [(bp, os.path.join(tmpdir, f"out_{i:04d}.fa"))
                     for i, bp in enumerate(batch_paths)]

            self.get_logger().info(
                "vcf-to-proteindb: dispatching %d variant-batch chunk(s) across %d worker(s)",
                len(tasks), min(workers, len(tasks)))

            # Pool initializer opens the GTF DB and FASTA index once per worker process,
            # not once per task — all batches handled by the same worker reuse them.
            with mp.get_context('spawn').Pool(
                min(workers, len(tasks)),
                initializer=_worker_init,
                initargs=(self.get_default_parameters(), pa, input_fasta, gene_annotations_gtf,
                          _annotation_indices),
            ) as pool:
                all_stats = pool.starmap(_vcf_to_proteindb_worker, tasks)

            # Concatenate the per-batch FASTAs into the final output.
            with open(self._proteindb_output, 'wb') as out:
                for _, batch_out in tasks:
                    if os.path.exists(batch_out):
                        with open(batch_out, 'rb') as f:
                            shutil.copyfileobj(f, out)

            # Aggregate per-batch counters and emit a single combined summary.
            combined: dict[str, int] = {}
            for stats in all_stats:
                for k, v in stats.items():
                    combined[k] = combined.get(k, 0) + v
            msg = "Translation summary ({} batch(es)):\n {}".format(
                len(all_stats),
                '\n'.join([x + ":" + str(combined[x]) for x in combined]))
            self.get_logger().info(msg)

        return self._proteindb_output

    @staticmethod
    def add_protein_to_map(seq: str, new_desc_string: str, protein_id: str, proteins: list, output_handle: Any) -> list:
        protein = {'description': new_desc_string, 'sequence': seq, 'accession': protein_id}
        proteins.append(protein)
        output_handle.write(">{}\t{}\n{}\n".format(protein_id, new_desc_string, seq))
        return proteins

    def check_proteindb(self, input_fasta: Optional[str] = None, add_stop_codon: bool = False, num_aa: int = 6) -> None:

        with open(input_fasta, 'r', encoding='utf-8') as input_handle:
            with open(self._proteindb_output, 'w', encoding='utf-8') as output_handle:
                proteins = []
                pcount = 0
                stop_count = 0
                gap_count = 0
                no_met = 0
                less = 0

                for record in SeqIO.parse(input_handle, 'fasta'):

                    seq = str(record.seq)
                    pcount += 1

                    # parse the description string into a dictionary
                    new_desc_string = record.description
                    new_desc_string = new_desc_string[new_desc_string.find(' ') + 1:]
                    # test for odd amino acids, stop codons, gaps
                    if not seq.startswith('M'):
                        no_met += 1
                    if seq.endswith('*'):
                        seq = seq[:-1]
                    if '-' in seq:
                        gap_count += 1
                        new_desc_string = new_desc_string + ' (Contains gaps)'
                    if '*' in seq:
                        stop_count += 1
                        if add_stop_codon:
                            seq_list = seq.split("*")
                            codon_index = 1
                            for codon in seq_list:
                                codon_description = new_desc_string + ' codon ' + str(codon_index)
                                protein_id = record.id + '_codon_' + str(codon_index)
                                seq = codon
                                if len(seq) > num_aa:
                                    proteins = self.add_protein_to_map(seq, codon_description, protein_id, proteins,
                                                                       output_handle)
                                codon_index = codon_index + 1
                        else:
                            cut = seq.index('*')
                            string = ' (Premature stop %s/%s)' % (cut, len(seq))
                            new_desc_string = new_desc_string + string
                            seq = seq[:cut]
                            # save the protein in list
                            if len(seq) > num_aa:
                                protein_id = record.id
                                proteins = self.add_protein_to_map(seq, new_desc_string, protein_id, proteins, output_handle)
                            else:
                                less += 1
                    else:
                        if len(seq) > num_aa:
                            protein_id = record.id
                            proteins = self.add_protein_to_map(seq, new_desc_string, protein_id, proteins, output_handle)
                        else:
                            less += 1

                self.get_logger().info("   translations that do not start with Met: %s", no_met)
                self.get_logger().info("   translations that have premature stop codons: %s", stop_count)
                self.get_logger().info("   translations that contain gaps: %s", gap_count)
                self.get_logger().info("   total number of input sequences was: %s", pcount)
                self.get_logger().info("   total number of sequences written was: %s", len(proteins))
        self.get_logger().info("   total number of proteins less than %s aminoacids: %s", num_aa, less)

    @staticmethod
    def write_output(seq_id: str, desc: str, seqs: list, prots_fn: Any, seqs_filter: Optional[list] = None) -> None:
        """Write ORFs to a FASTA output file handle.  Delegates to :func:`pgatk.toolbox.vcf_utils.write_output`."""
        _write_output(seq_id, desc, seqs, prots_fn, seqs_filter)


if __name__ == '__main__':
    logging.getLogger(__name__).error("This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
