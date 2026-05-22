from __future__ import annotations

import gzip
import logging
import re
import sqlite3
from pathlib import Path
from typing import Any, Optional

import gffutils
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.Seq import Seq

from pgatk.cgenomes.models import SNP
from pgatk.toolbox.general import ParameterConfiguration
from pgatk.toolbox.vcf_utils import get_altseq, get_orfs_vcf


# NCBI RefSeq accession → human chromosome name (GRCh37 and GRCh38).
# NC_000001–NC_000022 are autosomes 1–22; NC_000023=X, NC_000024=Y, NC_012920=MT.
_NCBI_CHR_MAP: dict[str, str] = {
    f"NC_{i:06d}": str(j) for i, j in zip(range(1, 23), range(1, 23))
}
_NCBI_CHR_MAP.update({"NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT"})

# Some MAF files (e.g. TCGA GDC) use numeric aliases for sex chromosomes.
_MAF_CHR_ALIASES: dict[str, str] = {"23": "X", "24": "Y", "25": "MT"}

_CBIO_BATCH_SIZE = 2000

# Per-worker state populated once by _cbio_worker_init and reused across all
# batches handled by the same worker process.
_cbio_worker_state: dict = {}


def _cbio_worker_init(
    seq_dic: dict,
    gene_tids: dict,
    gff_db_path: Optional[str],
    config: dict,
) -> None:
    """Pool initializer: cache FASTA data and open the GFF DB once per worker.

    Called once per spawned process.  Populates _cbio_worker_state so that
    every batch handled by that process reuses the same in-memory structures
    rather than re-opening or re-pickling them for each task.
    gff_db_path is the pre-built .db file path (built in the main process
    before the pool starts to avoid N workers racing on the same file).
    """
    _cbio_worker_state['seq_dic'] = seq_dic
    _cbio_worker_state['gene_tids'] = gene_tids
    _cbio_worker_state['gff_cache'] = {}
    _cbio_worker_state['gff_db'] = (
        gffutils.FeatureDB(gff_db_path) if gff_db_path else None
    )
    _cbio_worker_state.update(config)


def _cbio_translate_batch(rows: list, batch_idx: int) -> list:
    """Translate a batch of pre-parsed mutation rows.

    Each element of `rows` is a dict produced by
    CancerGenomesService._parse_maf_rows with keys:
        idx, gene, enst, pos, aa_mut, vartype, varclass, group,
        var_chrom, var_start, ref_allele, alt_allele

    Returns [(header, protein_str, group), ...] — one tuple per translated
    reading frame.  The caller writes these to the output FASTA.
    """
    seq_dic: dict = _cbio_worker_state['seq_dic']
    gene_tids: dict = _cbio_worker_state['gene_tids']
    gff_db = _cbio_worker_state.get('gff_db')
    gff_cache: dict = _cbio_worker_state['gff_cache']
    translation_table: int = _cbio_worker_state['translation_table']
    include_biotypes: list = _cbio_worker_state['include_biotypes']
    exclude_biotypes: list = _cbio_worker_state['exclude_biotypes']
    skip_including_all_cds: bool = _cbio_worker_state['skip_including_all_cds']

    nucleotide = ["A", "T", "C", "G"]
    log = logging.getLogger(__name__)
    results: list = []

    for row in rows:
        idx = row['idx']
        gene = row['gene']
        enst = row['enst']
        pos = row['pos']
        aa_mut = row['aa_mut']
        vartype = row['vartype']
        varclass = row['varclass']
        group = row.get('group')

        entry = seq_dic.get(enst)
        _gene_name_candidates: Optional[list] = None
        seq = None
        cds_info: list = []

        if entry is None:
            if not pos and gff_db is not None:
                _gene_name_candidates = gene_tids.get(gene, [])
                if not _gene_name_candidates:
                    log.debug(
                        "No FASTA record or GFF candidates for gene '%s' "
                        "(Transcript_ID=%s); skipping line %d", gene, enst, idx)
                    continue
            else:
                log.warning(
                    "No matching record for gene (%s) from row %s in FASTA file",
                    enst, idx)
                continue
        else:
            cds_info = entry['cds']
            biotype = entry['biotype']

            if not (bool(cds_info) and not skip_including_all_cds):
                if biotype:
                    if (biotype in exclude_biotypes or
                            (biotype not in include_biotypes and
                             include_biotypes != ['all'])):
                        continue

            seq = entry['seq'][cds_info[0] - 1: cds_info[1]] if cds_info else entry['seq']

        seq_mut = ""

        if pos:
            if ":" in pos:
                cdna_pos = pos.split(":")[1]
            else:
                cdna_pos = pos

            if vartype == "SNP":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos)[0])
                except IndexError:
                    log.warning("Incorrect SNP format or record %s %s", idx, pos)
                    continue
                if ">" not in pos:
                    log.warning("SNP position string missing '>' (line %s): %s", idx, pos)
                    continue
                arrow = pos.index(">")
                ref_dna = pos[arrow - 1]
                mut_dna = pos[arrow + 1]
                if mut_dna not in nucleotide:
                    log.warning("%s is not a nucleotide base %s", mut_dna, pos)
                    continue
                try:
                    if ref_dna == seq[enst_pos - 1]:
                        seq_mut = seq[:enst_pos - 1] + mut_dna + seq[enst_pos:]
                    else:
                        log.warning(
                            "incorrect substitution, unmatched nucleotide %s %s",
                            pos, enst)
                except IndexError:
                    log.warning("incorrect substitution, out of index %s", pos)
            elif vartype == "DEL":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                except IndexError:
                    log.warning("incorrect del format or record %s %s", idx, pos)
                    continue
                del_dna = pos.split("del")[1]
                if del_dna == seq[enst_pos - 1:enst_pos - 1 + len(del_dna)]:
                    seq_mut = seq[:enst_pos - 1] + seq[enst_pos - 1 + len(del_dna):]
                else:
                    log.warning("incorrect deletion, unmatched nucleotide %s", pos)
            elif vartype == "INS":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                except IndexError:
                    log.warning("incorrect ins/dup format or record %s %s", idx, pos)
                    continue
                if "ins" in pos:
                    ins_dna = pos.split("ins")[1]
                elif "dup" in pos:
                    ins_dna = pos.split("dup")[1]
                    if len(ins_dna) > 1:
                        enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[1])[0])
                else:
                    log.warning("unexpected insertion format at line %s", idx)
                    continue
                seq_mut = seq[:enst_pos] + ins_dna + seq[enst_pos:]

        elif gff_db is not None:
            var_chrom = row.get('var_chrom')
            var_start = row.get('var_start')
            raw_ref = row.get('ref_allele')
            raw_alt = row.get('alt_allele')

            if any(v is None for v in (var_chrom, var_start, raw_ref, raw_alt)):
                log.warning(
                    "Coordinate columns absent — cannot apply coordinate fallback "
                    "at line %d", idx)
                continue

            ref_allele = Seq("" if raw_ref == "-" else raw_ref)
            alt_allele = Seq("" if raw_alt == "-" else raw_alt)
            maf_chrom_raw = str(var_chrom).lstrip("chr")
            maf_chrom_norm = _MAF_CHR_ALIASES.get(maf_chrom_raw, maf_chrom_raw)

            candidates = _gene_name_candidates if _gene_name_candidates is not None else [enst]
            resolved = None
            for cand in candidates:
                cand_entry = seq_dic.get(cand) if _gene_name_candidates is not None else entry
                if cand_entry is None:
                    continue
                if cand in gff_cache:
                    feat_chrom, feat_strand, features_info = gff_cache[cand]
                else:
                    feat_chrom, feat_strand, features_info = \
                        CancerGenomesService._get_gff_features(gff_db, cand)
                    gff_cache[cand] = (feat_chrom, feat_strand, features_info)
                if feat_chrom is None:
                    continue
                gff_chrom_norm = _ncbi_accession_to_chrom(feat_chrom).lstrip("chr")
                if gff_chrom_norm != maf_chrom_norm:
                    if _gene_name_candidates is None:
                        log.warning(
                            "Chromosome mismatch for %s: GFF=%s MAF=%s",
                            cand, feat_chrom, var_chrom)
                    continue
                cand_cds = cand_entry['cds']
                coding_ref, coding_alt = get_altseq(
                    cand_entry['seq'], ref_allele, alt_allele,
                    var_start, feat_strand, features_info, cand_cds,
                )
                if coding_alt == "":
                    continue
                resolved = (cand, cand_entry, cand_cds, coding_ref, coding_alt)
                break

            if resolved is None:
                continue
            enst, entry, cds_info, seq, seq_mut = resolved

        else:
            log.debug(
                "Skipping %s %s (line %d): no HGVSc and no GFF annotation file",
                enst, aa_mut, idx)
            continue

        if seq_mut == "":
            continue

        num_orfs = 1 if cds_info else 3
        _, alt_orfs = get_orfs_vcf(
            Seq(str(seq)), Seq(str(seq_mut)), translation_table, num_orfs)

        for frame_idx, prot in enumerate(alt_orfs):
            if len(prot) <= 6:
                continue
            suffix = f'_RF{frame_idx + 1}' if num_orfs > 1 else ''
            header = f'cbiomut:{enst}:{gene}:{aa_mut}:{varclass}{suffix}'
            results.append((header, str(prot), group))

    return results


def _ncbi_accession_to_chrom(acc: str) -> str:
    """Convert 'NC_000001.10' → '1', 'NC_000023.11' → 'X', etc.
    Returns the original string when no mapping is found.
    """
    return _NCBI_CHR_MAP.get(acc.split(".")[0], acc)


def _open_text(path: str, mode: str = 'r', encoding: str = 'utf-8', **kwargs):
    """Open a file for text I/O, decompressing automatically if path ends with .gz."""
    enc = encoding or "utf-8"
    if path.endswith('.gz'):
        text_mode = mode if 't' in mode else mode + 't'
        return gzip.open(path, text_mode, encoding=enc, **kwargs)
    if enc == "utf-8":
        return open(path, mode, encoding="utf-8", **kwargs)
    if enc == "latin-1":
        return open(path, mode, encoding="latin-1", **kwargs)
    # Only utf-8 and latin-1 are used; fallback to utf-8 for any other value
    return open(path, mode, encoding="utf-8", **kwargs)


def _three_to_one(aa_str: str) -> str:
    """Convert a string of 3-letter amino acid codes to 1-letter codes.

    If the string is already in 1-letter codes (all single uppercase letters),
    it is returned unchanged.  For 3-letter codes (e.g. ``AlaGlyVal``), each
    triplet is converted using Biopython's mapping.
    """
    # Quick check: if all characters are uppercase single letters and the
    # string length is not a multiple of 3, or it contains lowercase, try
    # 3-letter conversion.
    if len(aa_str) >= 3 and any(c.islower() for c in aa_str):
        result = []
        for i in range(0, len(aa_str), 3):
            triplet = aa_str[i:i + 3]
            # Title-case the triplet for lookup (e.g. "ala" -> "Ala")
            triplet_title = triplet.capitalize()
            if triplet_title in protein_letters_3to1:
                result.append(protein_letters_3to1[triplet_title])
            else:
                return aa_str  # not a valid 3-letter sequence, return as-is
        return "".join(result)
    return aa_str


def _substitution_mismatch_detail(snp, seqs) -> str:
    """Return a human-readable string explaining a REF mismatch for substitution variants."""
    if ">" not in (snp.dna_mut or ""):
        return ""
    try:
        positions = re.findall(r'\d+', snp.dna_mut)
        if not positions:
            return ""
        tmplist = snp.dna_mut.split(">")
        expected_ref = re.sub("[^A-Z]+", "", tmplist[0])
        index = int(positions[0]) - 1
        found_bases = set()
        for seq in seqs:
            if index < len(seq):
                found_bases.add(str(seq[index]).upper())
        if found_bases:
            return f" [expected REF={expected_ref}, found={'/'.join(sorted(found_bases))} at pos {index + 1}]"
    except Exception:  # noqa: BLE001
        pass
    return ""


class CancerGenomesService(ParameterConfiguration):
    CONFIG_CANCER_GENOMES_MUTATION_FILE = 'mutation_file'
    CONFIG_COMPLETE_GENES_FILE = "all_cds_genes_file"
    CONFIG_OUTPUT_FILE = "output_file"
    CONFIG_COSMIC_DATA = "cosmic_data"
    CONFIG_KEY_DATA = 'proteindb'
    CONFIG_FILTER_INFO = 'filter_info'
    FILTER_COLUMN = "filter_column"
    ACCEPTED_VALUES = "accepted_values"
    SPLIT_BY_FILTER_COLUMN = "split_by_filter_column"
    WORKERS = 'workers'
    CBIO_BATCH_SIZE = 'cbio_batch_size'
    CLINICAL_SAMPLE_FILE = 'clinical_sample_file'
    CONFIG_COSMIC_SERVER = 'cosmic_server'
    CONFIG_GFF_FILE = 'gff_file'
    INCLUDE_BIOTYPES = 'include_biotypes'
    EXCLUDE_BIOTYPES = 'exclude_biotypes'
    SKIP_INCLUDING_ALL_CDS = 'skip_including_all_cds'
    BIOTYPE_STR = 'biotype_str'
    INCLUDE_VARIANT_CLASSIFICATIONS = 'include_variant_classifications'
    EXCLUDE_VARIANT_CLASSIFICATIONS = 'exclude_variant_classifications'
    TRANSLATION_TABLE = 'translation_table'

    def __init__(self, config_file: dict, pipeline_arguments: dict) -> None:
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CancerGenomesService, self).__init__(self.CONFIG_COSMIC_DATA, config_file, pipeline_arguments)

        self._filter_column = self.get_mutations_default_options(variable=self.FILTER_COLUMN,
                                                                 default_value='CANCER_TYPE')
        self._accepted_values = self.get_multiple_options(self.get_mutations_default_options(variable=self.ACCEPTED_VALUES, default_value="all"))

        self._split_by_filter_column = self.get_mutations_default_options(variable=self.SPLIT_BY_FILTER_COLUMN,
                                                                          default_value=False)
        self._local_clinical_sample_file = self.get_mutations_default_options(variable=self.CLINICAL_SAMPLE_FILE,
                                                                              default_value='')

        self._local_mutation_file = 'CosmicMutantExport.tsv.gz'
        if self.CONFIG_CANCER_GENOMES_MUTATION_FILE in self.get_pipeline_parameters():
            self._local_mutation_file = self.get_pipeline_parameters()[self.CONFIG_CANCER_GENOMES_MUTATION_FILE]
        elif self.CONFIG_COSMIC_DATA in self.get_default_parameters() and \
                self.CONFIG_COSMIC_SERVER in self.get_default_parameters()[self.CONFIG_COSMIC_DATA] and \
                self.CONFIG_CANCER_GENOMES_MUTATION_FILE in self.get_default_parameters()[self.CONFIG_COSMIC_DATA][
            self.CONFIG_COSMIC_SERVER]:
            self._local_mutation_file = \
                self.get_default_parameters()[self.CONFIG_COSMIC_DATA][self.CONFIG_COSMIC_SERVER][
                    self.CONFIG_CANCER_GENOMES_MUTATION_FILE]

        self._local_complete_genes = 'All_COSMIC_Genes.fasta.gz'
        if self.CONFIG_COMPLETE_GENES_FILE in self.get_pipeline_parameters():
            self._local_complete_genes = self.get_pipeline_parameters()[self.CONFIG_COMPLETE_GENES_FILE]
        elif self.CONFIG_COSMIC_DATA in self.get_default_parameters() and \
                self.CONFIG_COSMIC_SERVER in self.get_default_parameters()[self.CONFIG_COSMIC_DATA] and \
                self.CONFIG_COMPLETE_GENES_FILE in self.get_default_parameters()[self.CONFIG_COSMIC_DATA][
            self.CONFIG_COSMIC_SERVER]:
            self._local_complete_genes = \
                self.get_default_parameters()[self.CONFIG_COSMIC_DATA][self.CONFIG_COSMIC_SERVER][
                    self.CONFIG_COMPLETE_GENES_FILE]

        self._local_output_file = 'output_database.fasta'
        if self.CONFIG_OUTPUT_FILE in self.get_pipeline_parameters():
            self._local_output_file = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_FILE]

        self._local_gff_file = ''
        if self.CONFIG_GFF_FILE in self.get_pipeline_parameters():
            self._local_gff_file = self.get_pipeline_parameters()[self.CONFIG_GFF_FILE]

        self._biotype_str = self.get_mutations_default_options(self.BIOTYPE_STR, 'gene_biotype')
        self._include_biotypes = self.get_multiple_options(
            self.get_mutations_default_options(self.INCLUDE_BIOTYPES, 'protein_coding'))
        self._exclude_biotypes = self.get_multiple_options(
            self.get_mutations_default_options(self.EXCLUDE_BIOTYPES, ''))
        self._skip_including_all_cds = self.get_mutations_default_options(self.SKIP_INCLUDING_ALL_CDS, False)
        self._include_variant_classifications = self.get_multiple_options(
            self.get_mutations_default_options(self.INCLUDE_VARIANT_CLASSIFICATIONS, 'all'))
        self._exclude_variant_classifications = self.get_multiple_options(
            self.get_mutations_default_options(self.EXCLUDE_VARIANT_CLASSIFICATIONS, 'Nonsense_Mutation'))
        self._translation_table = int(self.get_mutations_default_options(self.TRANSLATION_TABLE, 1))
        try:
            self._workers = max(1, int(
                self.get_mutations_default_options(self.WORKERS, 1) or 1))
        except (TypeError, ValueError):
            self._workers = 1
        try:
            self._cbio_batch_size = max(1, int(
                self.get_mutations_default_options(self.CBIO_BATCH_SIZE, _CBIO_BATCH_SIZE) or _CBIO_BATCH_SIZE))
        except (TypeError, ValueError):
            self._cbio_batch_size = _CBIO_BATCH_SIZE

    def get_mutations_default_options(self, variable: str, default_value: Any) -> Any:
        return_value = default_value
        if variable in self.get_pipeline_parameters():
            return_value = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_DATA in self.get_default_parameters() \
                and self.CONFIG_FILTER_INFO in self.get_default_parameters()[self.CONFIG_KEY_DATA] \
                and variable in self.get_default_parameters()[
            self.CONFIG_KEY_DATA][self.CONFIG_FILTER_INFO]:
            return_value = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_FILTER_INFO][variable]
        return return_value

    @staticmethod
    def get_multiple_options(options_str: str) -> list[str]:
        """
        This method takes an String like option1, option2, ... and produce and array [option1, option2,... ]
        :param options_str:
        :return: Array
        """
        return list(map(lambda x: x.strip(), options_str.split(",")))

    @staticmethod
    def get_mut_pro_seq(snp: SNP, seq: Seq) -> Optional[str]:
        nucleotide = ["A", "T", "C", "G"]
        mut_pro_seq = ""
        # A small number of COSMIC gene FASTA entries have a leading N (masked
        # nucleotide) before the real CDS start.  HGVS c. positions count from
        # the first coding base, so strip any leading non-ACGT characters to
        # keep position lookups consistent with the annotation.
        leading = len(seq) - len(str(seq).lstrip('Nn'))
        if leading:
            seq = seq[leading:]
        if (snp.dna_mut is not None and "?" not in snp.dna_mut and
                snp.aa_mut is not None and snp.aa_mut != 'p.?'):  # unambiguous DNA change known in CDS sequence
            positions = re.findall(r'\d+', snp.dna_mut)
            if ">" in snp.dna_mut and len(positions) == 1:  # Substitution
                tmplist = snp.dna_mut.split(">")
                ref_dna = re.sub("[^A-Z]+", "", tmplist[0])
                mut_dna = re.sub("[^A-Z]+", "", tmplist[1])
                index = int(positions[0]) - 1
                if ref_dna == str(seq[index]).upper() and mut_dna in nucleotide:
                    seq_mut = seq[:index] + mut_dna + seq[index + 1:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))
                else:
                    return None  # REF base doesn't match FASTA — transcript version mismatch
            elif "delins" in snp.dna_mut:
                # HGVS delins: one or more nucleotides replaced by one or more other nucleotides.
                coord_part, insert_raw = snp.dna_mut.split("delins", 1)
                insert_dna = insert_raw.upper()
                if re.search(r'\d[+-]|\*|-\d', coord_part):
                    # Intronic (c.N+X, c.N-X), 5'UTR (c.-N) or 3'UTR (c.*N) offsets
                    # cannot be mapped onto a CDS-only FASTA sequence; skip.
                    pass
                elif not insert_dna.isalpha():
                    # Unknown ('?'), N[n], or conversion-notation insertions cannot be resolved.
                    pass
                elif len(positions) >= 2:
                    # Range delins: delete positions[0]..positions[1] (1-indexed, inclusive).
                    del_index1 = int(positions[0]) - 1  # 0-indexed start
                    del_index2 = int(positions[1])      # 0-indexed exclusive end
                    seq_mut = seq[:del_index1] + insert_dna + seq[del_index2:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))
                elif len(positions) == 1:
                    # Single-position delins: replace the one nucleotide at positions[0].
                    del_index1 = int(positions[0]) - 1
                    seq_mut = seq[:del_index1] + insert_dna + seq[del_index1 + 1:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))

            elif "ins" in snp.dna_mut:
                index = snp.dna_mut.index("ins")
                insert_dna = snp.dna_mut[index + 3:]
                if insert_dna.isalpha():
                    ins_index1 = int(positions[0])
                    seq_mut = seq[:ins_index1] + insert_dna + seq[ins_index1:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))

            elif "dup" in snp.dna_mut:
                # Tandem duplication: re-insert the duplicated range after its end position.
                if len(positions) == 2:
                    dup_start = int(positions[0]) - 1
                    dup_end = int(positions[1])
                    dup_seq = str(seq[dup_start:dup_end])
                    seq_mut = seq[:dup_end] + dup_seq + seq[dup_end:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))
                elif len(positions) == 1:
                    dup_pos = int(positions[0])
                    dup_seq = str(seq[dup_pos - 1:dup_pos])
                    seq_mut = seq[:dup_pos] + dup_seq + seq[dup_pos:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))

            elif "del" in snp.dna_mut:
                if len(positions) == 2:
                    del_index1 = int(positions[0]) - 1
                    del_index2 = int(positions[1])
                    seq_mut = seq[:del_index1] + seq[del_index2:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))
                elif len(positions) == 1:
                    del_index1 = int(positions[0]) - 1
                    seq_mut = seq[:del_index1] + seq[del_index1 + 1:]
                    mut_pro_seq = str(seq_mut.translate(to_stop=False))
        else:
            if snp.aa_mut is not None and "?" not in snp.aa_mut:  # unambiguous aa change known in protein sequence
                positions = re.findall(r'\d+', snp.aa_mut)
                protein_seq = str(seq.translate(to_stop=False))

                if "Missense" in snp.mutation_type or "missense_variant" in snp.mutation_type:
                    # Extract the mutant residue from HGVS like p.V600E (1-letter)
                    # or p.Val600Glu (3-letter).  For 1-letter the last char is the
                    # AA; for 3-letter we need to convert the last triplet.
                    aa_suffix = re.sub(r'^p\..*\d+', '', snp.aa_mut)
                    mut_aa = _three_to_one(aa_suffix) if aa_suffix else ''
                    if not mut_aa or len(mut_aa) != 1 or not mut_aa.isalpha():
                        return ''
                    index = int(positions[0]) - 1
                    mut_pro_seq = protein_seq[:index] + mut_aa + protein_seq[index + 1:]
                elif "Nonsense" in snp.mutation_type or "stop_gained" in snp.mutation_type:
                    index = int(positions[0]) - 1
                    mut_pro_seq = protein_seq[:index]
                elif "Insertion - In frame" in snp.mutation_type or "inframe_insertion" in snp.mutation_type:
                    if "dup" in snp.aa_mut:
                        # Protein-level tandem dup: re-insert the duplicated residues after the range.
                        dup_start = int(positions[0]) - 1
                        dup_end = int(positions[-1])
                        dup_aa = protein_seq[dup_start:dup_end]
                        mut_pro_seq = protein_seq[:dup_end] + dup_aa + protein_seq[dup_end:]
                    else:
                        try:
                            index = snp.aa_mut.index("ins")
                        except ValueError:
                            return ''
                        insert_aa_raw = snp.aa_mut[index + 3:]
                        insert_aa = _three_to_one(insert_aa_raw)
                        if insert_aa.isalpha():
                            ins_index1 = int(positions[0])
                            mut_pro_seq = protein_seq[:ins_index1] + insert_aa + protein_seq[ins_index1:]
                elif "Deletion - In frame" in snp.mutation_type or "inframe_deletion" in snp.mutation_type:
                    if len(positions) == 2:
                        del_index1 = int(positions[0]) - 1
                        del_index2 = int(positions[1])
                        mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index2:]
                    elif len(positions) == 1:
                        del_index1 = int(positions[0]) - 1
                        mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index1 + 1:]
                elif ("Complex" in snp.mutation_type or "protein_altering_variant" in snp.mutation_type) \
                        and "frameshift" not in snp.mutation_type:
                    try:
                        index = snp.aa_mut.index(">")
                    except ValueError:
                        return ''
                    mut_aa_raw = snp.aa_mut[index + 1:]
                    mut_aa = _three_to_one(mut_aa_raw.replace('*', ''))
                    if '*' in mut_aa_raw:
                        mut_aa += '*'
                    if not mut_aa.replace('*','').isalpha():
                        return ''
                    if "deletion" in snp.mutation_type:
                        del_index1 = int(positions[0]) - 1
                        del_index2 = int(positions[1])
                        mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]

                    elif "insertion" in snp.mutation_type:
                        ins_index1 = int(positions[0])
                        mut_pro_seq = protein_seq[:ins_index1] + mut_aa + protein_seq[ins_index1:]
                    elif "compound substitution" in snp.mutation_type:
                        if "*" not in mut_aa:
                            del_index1 = int(positions[0]) - 1
                            del_index2 = int(positions[1])
                            mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                        else:
                            del_index1 = int(positions[0]) - 1
                            mut_pro_seq = protein_seq[:del_index1] + mut_aa.replace("*", "")

        return mut_pro_seq

    def cosmic_to_proteindb(self) -> None:
        """
        This function translates the mutation file + COSMIC genes into a protein Fasta database. The
        method writes into the file system the output Fasta.
        :return:
        """
        self.get_logger().debug("Starting reading All cosmic genes")
        COSMIC_CDS_DB = {}
        with _open_text(self._local_complete_genes, encoding='utf-8') as genes_handle:
            for record in SeqIO.parse(genes_handle, 'fasta'):
                try:
                    COSMIC_CDS_DB[record.id].append(record)
                except KeyError:
                    COSMIC_CDS_DB[record.id] = [record]

        # Build phenotype-to-group mapping from the classification file if provided.
        # COSMIC_PHENOTYPE_ID is present directly in the mutation file (column 6) and
        # is the join key to the Cosmic_Classification file where PRIMARY_SITE lives.
        sample_groups_dict = {}
        if self._local_clinical_sample_file and self._filter_column:
            sample_groups_dict = self.get_value_per_sample(
                self._local_clinical_sample_file,
                self._filter_column,
                sample_id_column='COSMIC_PHENOTYPE_ID',
            )
            if not sample_groups_dict:
                self.get_logger().warning(
                    "clinical_sample_file '%s' produced no phenotype mappings; "
                    "falling back to direct filter-column lookup in the mutation file.",
                    self._local_clinical_sample_file,
                )

        regex = re.compile('[^a-zA-Z]')
        mutation_dic = {}
        groups_mutations_dict = {}
        ref_mismatch_count = 0
        unsupported_count = 0
        self.get_logger().debug("Reading input CosmicMutantExport.tsv ...")
        line_counter = 1

        required_columns = ["GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "MUTATION_CDS", "MUTATION_AA", "MUTATION_DESCRIPTION"]
        with _open_text(self._local_mutation_file, encoding="latin-1") as cosmic_input, \
             open(self._local_output_file, 'w', encoding='utf-8') as output:
            col_header = cosmic_input.readline().strip().split("\t")
            try:
                gene_col = col_header.index("GENE_SYMBOL")
                enst_col = col_header.index("TRANSCRIPT_ACCESSION")
                cds_col = col_header.index("MUTATION_CDS")
                aa_col = col_header.index("MUTATION_AA")
                muttype_col = col_header.index("MUTATION_DESCRIPTION")
            except ValueError as e:
                self.get_logger().error(
                    "COSMIC file missing required columns. Expected %s, got: %s",
                    required_columns, col_header
                )
                raise ValueError(
                    f"COSMIC mutation file missing required column: {e}. "
                    f"Expected columns include: {required_columns}"
                ) from e

            # COSMIC_PHENOTYPE_ID is the direct join key to the classification file.
            try:
                phenotype_id_col = col_header.index("COSMIC_PHENOTYPE_ID")
            except ValueError:
                phenotype_id_col = None
                if sample_groups_dict:
                    self.get_logger().warning(
                        "COSMIC_PHENOTYPE_ID not found in mutation file header; classification file join will be skipped."
                    )

            # Fall back to a direct filter column in the mutation file when no classification file is provided.
            filter_col = None
            if not sample_groups_dict and self._filter_column:
                try:
                    filter_col = col_header.index(self._filter_column)
                except ValueError:
                    self.get_logger().warning(
                        "Filter column '%s' not found in COSMIC header. Filtering disabled.",
                        self._filter_column,
                    )

            max_col = max(gene_col, enst_col, cds_col, aa_col, muttype_col)
            if filter_col is not None:
                max_col = max(max_col, filter_col)
            if phenotype_id_col is not None and sample_groups_dict:
                max_col = max(max_col, phenotype_id_col)

            for line in cosmic_input:
                if line_counter % 10000 == 0:
                    self.get_logger().debug("Number of lines finished -- '%s'", line_counter)
                line_counter += 1
                row = line.strip().split("\t")
                if len(row) <= max_col:
                    self.get_logger().debug("Skipping malformed row (insufficient columns) at line %s: %s", line_counter, row[:5])
                    continue

                # Determine the group (e.g. primary site) for this mutation.
                group = None
                if sample_groups_dict and phenotype_id_col is not None:
                    group = sample_groups_dict.get(row[phenotype_id_col])
                    if group is None and (self._accepted_values != ['all'] or self._split_by_filter_column):
                        self.get_logger().warning(
                            "No classification found for COSMIC_PHENOTYPE_ID '%s'; skipping row.",
                            row[phenotype_id_col],
                        )
                        continue
                elif filter_col is not None:
                    group = row[filter_col]

                if group is not None and group not in self._accepted_values and self._accepted_values != ['all']:
                    continue

                if "coding silent" in row[muttype_col] or "synonymous_variant" in row[muttype_col]:
                    continue

                # Skip non-coding mutations with unknown protein consequence (UTR,
                # intronic, splice-region variants etc.). These can never produce
                # a mutant protein sequence so there is nothing useful to warn about.
                if row[aa_col] == 'p.?':
                    continue

                snp = SNP(gene=row[gene_col], mrna=row[enst_col], dna_mut=row[cds_col], aa_mut=row[aa_col],
                          mutation_type=row[muttype_col])
                fasta_header = "COSMIC:%s:%s:%s" % (snp.gene, snp.aa_mut, snp.mutation_type.replace(" ", ""))
                try:
                    this_gene_records = COSMIC_CDS_DB[snp.gene]
                    seqs = []
                    for record in this_gene_records:
                        seqs.append(record.seq)

                except KeyError:  # geneID is not in All_COSMIC_Genes.fasta
                    continue

                failure_reason = "unsupported"
                mut_pro_seq = ""
                for seq in seqs:
                    try:
                        result = self.get_mut_pro_seq(snp, seq)
                    except (IndexError, ValueError):
                        failure_reason = "ref_mismatch"
                        continue
                    if result is None:
                        failure_reason = "ref_mismatch"
                        continue
                    if result:
                        mut_pro_seq = result
                        break

                if mut_pro_seq:
                    entry = ">%s\n%s\n" % (fasta_header, mut_pro_seq)
                    if fasta_header not in mutation_dic:
                        output.write(entry)
                        mutation_dic[fasta_header] = 1

                    if self._split_by_filter_column and group is not None:
                        try:
                            groups_mutations_dict[group][fasta_header] = entry
                        except KeyError:
                            groups_mutations_dict[group] = {fasta_header: entry}
                else:
                    if failure_reason == "ref_mismatch":
                        ref_mismatch_count += 1
                        if self.get_logger().isEnabledFor(logging.DEBUG):
                            detail = _substitution_mismatch_detail(snp, seqs)
                            self.get_logger().debug(
                                "Skipped (REF mismatch): gene=%s, dna_mut=%s, "
                                "aa_mut=%s, mutation_type=%s%s",
                                snp.gene, snp.dna_mut, snp.aa_mut, snp.mutation_type,
                                detail,
                            )
                    else:
                        unsupported_count += 1
                        self.get_logger().debug(
                            "Skipped (unsupported HGVS): gene=%s, dna_mut=%s, "
                            "aa_mut=%s, mutation_type=%s",
                            snp.gene, snp.dna_mut, snp.aa_mut, snp.mutation_type,
                        )

        if ref_mismatch_count:
            self.get_logger().warning(
                "%d mutation records skipped: REF base in COSMIC does not match the "
                "gene FASTA (transcript version mismatch); run with DEBUG logging to see details.",
                ref_mismatch_count,
            )
        if unsupported_count:
            self.get_logger().warning(
                "%d mutation records skipped: HGVS notation not supported by the parser; "
                "run with DEBUG logging to see details.",
                unsupported_count,
            )

        for group_name in groups_mutations_dict.keys():
            output_base = str(Path(self._local_output_file).with_suffix(''))
            with open(f"{output_base}_{regex.sub('', group_name)}.fa", 'w', encoding='utf-8') as fn:
                for fasta_hdr in groups_mutations_dict[group_name].keys():
                    fn.write(groups_mutations_dict[group_name][fasta_hdr])

        self.get_logger().debug("COSMIC contains in total {} non redundant mutations".format(len(mutation_dic)))

    @staticmethod
    def get_sample_headers(header_line: list, filter_column: str,
                           sample_id_column: str = 'SAMPLE_ID') -> tuple[Optional[int], Optional[int]]:
        _logger = logging.getLogger(__name__)
        try:
            filter_col = header_line.index(filter_column)
        except ValueError:
            _logger.warning('%s was not found in the header row: %s', filter_column, header_line)
            return None, None
        try:
            sample_id_col = header_line.index(sample_id_column)
        except ValueError:
            _logger.warning('%s was not found in the header row: %s', sample_id_column, header_line)
            return None, None
        return filter_col, sample_id_col

    def get_value_per_sample(self, local_clinical_sample_file: str, filter_column: str,
                              sample_id_column: str = 'SAMPLE_ID') -> dict:
        sample_value = {}
        if local_clinical_sample_file:
            with _open_text(local_clinical_sample_file, encoding='utf-8') as clin_fn:
                filter_column_col, sample_id_col = None, None
                for line in clin_fn.readlines():
                    if line.startswith('#'):
                        continue
                    sl = line.strip().split('\t')
                    if sample_id_column in sl and filter_column in sl:
                        filter_column_col, sample_id_col = self.get_sample_headers(
                            sl, filter_column, sample_id_column
                        )
                        continue
                    if filter_column_col is not None and sample_id_col is not None:
                        if (sample_id_col < len(sl) and filter_column_col < len(sl)):
                            sample_value[sl[sample_id_col]] = sl[filter_column_col].strip().replace(' ', '_')
                    else:
                        self.get_logger().warning("No column was found for %s, %s in %s", filter_column,
                                                  sample_id_column, local_clinical_sample_file)
        return sample_value

    @staticmethod
    def get_mut_header_cols(header_cols: dict, row: list) -> dict:
        for col in header_cols.keys():
            header_cols[col] = row.index(col)

        return header_cols

    @staticmethod
    def _load_gff_db(gff_file: str) -> gffutils.FeatureDB:
        """Build or load a gffutils SQLite database from a GFF3/GTF file.

        The DB is stored as ``<gff_file>.db`` alongside the source file.
        Building is slow on first run (NCBI GFF3 ~10–30 min); subsequent
        runs load the cached SQLite file in seconds.
        """
        log = logging.getLogger(__name__)
        db_file = str(Path(gff_file).with_suffix(".db"))

        def _build(force: bool) -> None:
            gffutils.create_db(
                gff_file, db_file,
                merge_strategy="create_unique",
                keep_order=True,
                disable_infer_transcripts=True,
                disable_infer_genes=True,
                verbose=False,
                force=force,
            )

        if not Path(db_file).exists():
            log.warning(
                "Building gffutils database from %s — this takes 10–30 minutes on "
                "first run and is cached as %s for subsequent runs.", gff_file, db_file
            )
            _build(force=False)
        else:
            try:
                return gffutils.FeatureDB(db_file)
            except Exception:
                log.warning(
                    "gffutils DB at %s appears corrupt — rebuilding from %s "
                    "(this may take 10–30 minutes).", db_file, gff_file
                )
                _build(force=True)

        return gffutils.FeatureDB(db_file)

    @staticmethod
    def _get_gff_features(
        db: gffutils.FeatureDB, transcript_id: str, feature_types: Optional[list] = None
    ) -> tuple:
        """Return (chrom, strand, features_info) for a transcript.

        Tries the ID directly, then without version, then with the GFF3
        NCBI ``rna-`` prefix.  Returns (None, None, None) when not found.
        """
        if feature_types is None:
            feature_types = ["CDS"]
        feature = None
        bare = transcript_id.split(".")[0]
        for candidate in (
            transcript_id,
            bare,
            f"rna-{transcript_id}",
            f"rna-{bare}",
        ):
            try:
                feature = db[candidate]
                break
            except gffutils.exceptions.FeatureNotFoundError:
                continue
        if feature is None:
            # Last resort: NCBI GFF3 uses versioned IDs like 'rna-NM_001255.3'.
            # Query the SQLite backing store with a LIKE pattern on the ID column.
            try:
                rows = db.conn.execute(
                    "SELECT id FROM features "
                    "WHERE featuretype IN ('mRNA','transcript','primary_transcript') "
                    "AND id LIKE ? LIMIT 1",
                    [f"rna-{bare}.%"],
                ).fetchall()
                if rows:
                    feature = db[rows[0][0]]
            except Exception:
                pass
        if feature is None:
            logging.getLogger(__name__).warning(
                "Transcript %s not found in GFF annotation database.", transcript_id
            )
            return None, None, None
        coding_features = [
            [f.start, f.end, f.featuretype]
            for f in db.children(feature, featuretype=feature_types, order_by="end")
        ]
        return feature.chrom, feature.strand, coding_features

    @staticmethod
    def _parse_cbio_fasta_header(description: str, biotype_str: str) -> tuple[list, str, str]:
        """Parse CDS=[start]-[end], biotype, and gene name from a gffread -F -w FASTA description.

        Returns (cds_info, biotype, gene_name). cds_info is [] when no CDS= token is present.
        Mirrors the CDS= parsing in ensembl.py:802-820.
        """
        cds_info = []
        biotype = ''
        gene_name = ''
        for token in description.split(' '):
            stripped = token.strip('[]')
            if stripped.startswith('CDS='):
                try:
                    cds_info = [int(x) for x in stripped.split('=')[1].split('-')]
                except (ValueError, IndexError):
                    pass
                break
        for attr in description.split(';'):
            attr = attr.strip()
            if attr.startswith(f'{biotype_str}='):
                biotype = attr.split('=', 1)[1]
            elif attr.startswith('gene='):
                gene_name = attr.split('=', 1)[1]
        return cds_info, biotype, gene_name

    def cbioportal_to_proteindb(self) -> None:
        """Translate cBioportal mutation MAF into a protein FASTA database.

        When workers > 1 the translation phase is parallelised: the main
        process parses and pre-filters all rows (Phase 1), fans out fixed-size
        batches to a process pool (Phase 2), then writes the merged results and
        per-group split files (Phase 3).  workers=1 (default) runs Phase 2 in
        the main process with zero pool overhead, preserving the original
        single-process behaviour exactly.
        """
        regex = re.compile('[^a-zA-Z]')

        # ── FASTA loading ─────────────────────────────────────────────────────
        seq_dic: dict = {}
        _MAX_GENE_CANDIDATES = 10
        gene_tids_raw: dict[str, list] = {}
        with _open_text(self._local_complete_genes, encoding='utf-8') as genes_handle:
            for record in SeqIO.parse(genes_handle, "fasta"):
                tid = record.id.removeprefix('rna-').split('.')[0]
                if tid not in seq_dic:
                    cds_info, biotype, gene_name = self._parse_cbio_fasta_header(
                        record.description, self._biotype_str)
                    seq_dic[tid] = {'seq': record.seq, 'cds': cds_info, 'biotype': biotype}
                    if gene_name and cds_info:
                        cds_len = cds_info[1] - cds_info[0]
                        gene_tids_raw.setdefault(gene_name, []).append((cds_len, tid))
        gene_tids: dict[str, list] = {
            g: [tid for _, tid in sorted(tids, reverse=True)[:_MAX_GENE_CANDIDATES]]
            for g, tids in gene_tids_raw.items()
        }

        # ── Clinical sample / group filter ────────────────────────────────────
        sample_groups_dict: dict = {}
        if self._accepted_values != ['all'] or self._split_by_filter_column:
            if self._local_clinical_sample_file:
                sample_groups_dict = self.get_value_per_sample(
                    self._local_clinical_sample_file, self._filter_column)
                self.get_logger().info('sample_groups_dict %s %s',
                                       self._local_clinical_sample_file, self._filter_column)
                if sample_groups_dict == {}:
                    return
            else:
                self.get_logger().warning(
                    'No clinical sample file is given therefore no filter could be applied.')
                return

        # ── GFF DB: build once before workers start ───────────────────────────
        gff_db_path: Optional[str] = None
        if self._local_gff_file:
            self._load_gff_db(self._local_gff_file)
            gff_db_path = str(Path(self._local_gff_file).with_suffix('.db'))

        # ── Phase 1: parse + pre-filter all rows (serial, I/O-bound) ─────────
        parsed_rows = self._parse_maf_rows(sample_groups_dict)
        if not parsed_rows:
            open(self._local_output_file, 'w', encoding='utf-8').close()
            return

        worker_config = {
            'translation_table': self._translation_table,
            'include_biotypes': self._include_biotypes,
            'exclude_biotypes': self._exclude_biotypes,
            'skip_including_all_cds': self._skip_including_all_cds,
            'split_by_filter_column': self._split_by_filter_column,
        }

        # ── Phase 2: translate (sequential or parallel) ───────────────────────
        if self._workers <= 1:
            _cbio_worker_init(seq_dic, gene_tids, gff_db_path, worker_config)
            all_results = _cbio_translate_batch(parsed_rows, 0)
        else:
            all_results = self._run_cbio_parallel(
                parsed_rows, seq_dic, gene_tids, gff_db_path, worker_config)

        # ── Phase 3: write main output + per-group split files ────────────────
        # Collect into dicts first so duplicate headers (same mutation in multiple
        # patients) are collapsed to one entry before writing.
        main_mutations: dict = {}
        group_mutations_dict: dict = {}
        for header, prot, group in all_results:
            main_mutations[header] = prot
            if self._split_by_filter_column and group is not None:
                group_mutations_dict.setdefault(group, {})[header] = prot

        with open(self._local_output_file, 'w', encoding='utf-8') as output:
            for hdr, seq in main_mutations.items():
                output.write(f'>{hdr}\n{seq}\n')

        for group, mutations in group_mutations_dict.items():
            output_base = str(Path(self._local_output_file).with_suffix(''))
            with open(f"{output_base}_{regex.sub('', group)}.fa", 'w', encoding='utf-8') as fn:
                for hdr, seq in mutations.items():
                    fn.write(f'>{hdr}\n{seq}\n')

    def _parse_maf_rows(self, sample_groups_dict: dict) -> list:
        """Phase 1: read the mutation MAF, apply pre-translation filters, return row dicts.

        Applies: header detection, sample-group lookup, variant-classification
        filter.  Biotype filtering is deferred to the worker because it needs
        seq_dic which is not loaded here.

        Each returned dict has keys:
            idx, gene, enst, pos, aa_mut, vartype, varclass, group,
            var_chrom, var_start, ref_allele, alt_allele
        """
        required_cols = {
            "HGVSc": None, "Transcript_ID": None, "Variant_Classification": None,
            "Variant_Type": None, "HGVSp_Short": None, "Tumor_Sample_Barcode": None,
        }
        coord_cols = {
            "Chromosome": None, "Start_Position": None,
            "Reference_Allele": None, "Tumor_Seq_Allele2": None, "Strand": None,
        }
        header_cols = {**required_cols, **coord_cols}
        header_detected = False
        rows: list = []

        with _open_text(self._local_mutation_file, encoding='utf-8') as mutfile:
            for i, line in enumerate(mutfile):
                row = line.strip().split("\t")
                if row[0] == '#':
                    self.get_logger().info("skipping line (%s): %s", i, row)
                    continue
                if not header_detected and set(required_cols.keys()).issubset(set(row)):
                    for col in required_cols:
                        header_cols[col] = row.index(col)
                    for col in coord_cols:
                        header_cols[col] = row.index(col) if col in row else None
                    header_detected = True
                    continue
                if not header_detected:
                    self.get_logger().error("Incorrect header column is given")
                    continue

                group = None
                if self._accepted_values != ['all'] or self._split_by_filter_column:
                    try:
                        group = sample_groups_dict[row[header_cols['Tumor_Sample_Barcode']]]
                    except KeyError:
                        self.get_logger().warning(
                            "No clinical info was found for sample %s. Skipping (line %s).",
                            row[header_cols['Tumor_Sample_Barcode']], i)
                        continue
                    except IndexError:
                        self.get_logger().warning(
                            "No sampleID was found in (line %s): %s", i, row)
                if group not in self._accepted_values and self._accepted_values != ['all']:
                    continue

                try:
                    gene = row[0]
                    pos = row[header_cols["HGVSc"]]
                    enst = row[header_cols["Transcript_ID"]]
                    aa_mut = row[header_cols["HGVSp_Short"]]
                    vartype = row[header_cols["Variant_Type"]]
                    varclass = row[header_cols["Variant_Classification"]]
                except IndexError:
                    self.get_logger().warning("Incorrect line (%s): %s", i, row)
                    continue

                if self._include_variant_classifications != ['all']:
                    if varclass not in self._include_variant_classifications:
                        continue
                if varclass in self._exclude_variant_classifications:
                    continue

                # Extract coordinate columns (may be None if column absent from header)
                def _get_col(col_name):
                    ci = header_cols.get(col_name)
                    return row[ci] if ci is not None and ci < len(row) else None

                var_start = None
                raw_start = _get_col("Start_Position")
                if raw_start is not None:
                    try:
                        var_start = int(raw_start)
                    except (ValueError, TypeError):
                        pass

                rows.append({
                    'idx': i,
                    'gene': gene,
                    'enst': enst,
                    'pos': pos,
                    'aa_mut': aa_mut,
                    'vartype': vartype,
                    'varclass': varclass,
                    'group': group,
                    'var_chrom': _get_col("Chromosome"),
                    'var_start': var_start,
                    'ref_allele': _get_col("Reference_Allele"),
                    'alt_allele': _get_col("Tumor_Seq_Allele2"),
                })

        return rows

    def _run_cbio_parallel(
        self,
        rows: list,
        seq_dic: dict,
        gene_tids: dict,
        gff_db_path: Optional[str],
        config: dict,
    ) -> list:
        """Phase 2 (parallel): fan rows out to a process pool in fixed-size batches.

        Splits `rows` into batches of self._cbio_batch_size and dispatches them
        via multiprocessing.Pool.starmap.  Follows the same spawn-context +
        pool-initializer pattern as the VCF parallel path in ensembl.py so that
        seq_dic and the GFF DB are loaded once per worker, not once per batch.

        Falls back to a single in-process call when rows fit in one batch or
        worker count collapses to 1.
        """
        import multiprocessing as mp

        batch_size = self._cbio_batch_size
        batches = [rows[i:i + batch_size] for i in range(0, len(rows), batch_size)]
        n_workers = min(self._workers, len(batches))

        if n_workers <= 1 or len(batches) <= 1:
            _cbio_worker_init(seq_dic, gene_tids, gff_db_path, config)
            return _cbio_translate_batch(rows, 0)

        self.get_logger().info(
            "cbioportal-to-proteindb: dispatching %d batch(es) of ≤%d rows "
            "across %d worker(s)",
            len(batches), batch_size, n_workers)

        with mp.get_context('spawn').Pool(
            n_workers,
            initializer=_cbio_worker_init,
            initargs=(seq_dic, gene_tids, gff_db_path, config),
        ) as pool:
            batch_results = pool.starmap(
                _cbio_translate_batch,
                [(batch, idx) for idx, batch in enumerate(batches)],
            )

        return [item for batch in batch_results for item in batch]
