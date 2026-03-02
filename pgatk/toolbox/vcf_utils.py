"""
Shared VCF-to-protein utility functions.

These functions were extracted from ``pgatk.ensembl.ensembl.EnsemblDataService``
so that they can be reused by other pipelines (e.g. ClinVar/NCBI RefSeq)
without pulling in Ensembl-specific dependencies.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

from Bio.Seq import Seq

logger = logging.getLogger(__name__)


def check_overlap(var_start: int, var_end: int, features_info: Optional[list] = None) -> bool:
    """Return *True* when the variant overlaps any of the features.

    Parameters
    ----------
    var_start : int
        Genomic start position of the variant.
    var_end : int
        Genomic end position of the variant.
    features_info : list, optional
        List of ``[start, end, type]`` triples describing feature regions.
        Defaults to ``[[0, 1, 'type']]``.
    """
    if features_info is None:
        features_info = [[0, 1, 'type']]
    if var_start == -1:
        return True
    # check if the var overlaps any of the features
    for feature_pos in features_info:
        pep_start = feature_pos[0]
        pep_end = feature_pos[1]
        if var_start <= pep_start <= var_end:  # fully contained or partial overlap from the end
            return True
        elif var_start <= pep_end <= var_end:  # partial overlap in the begining
            return True
        elif pep_start <= var_start and pep_end >= var_end:  # fully covered
            return True
    return False


def get_altseq(
    ref_seq: str,
    ref_allele: str,
    var_allele: str,
    var_pos: int,
    strand: str,
    features_info: list,
    cds_info: Optional[list] = None,
) -> tuple:
    """Modify a reference sequence based on a variant allele.

    The given sequence in the FASTA file represents all exons of the transcript
    combined.  For protein-coding genes the CDS is specified, so the sequence
    position has to be calculated based on CDS positions.  For non-protein-coding
    genes the whole sequence is used.

    Parameters
    ----------
    ref_seq : str or Bio.Seq.Seq
        Reference transcript sequence.
    ref_allele : str or Bio.Seq.Seq
        Reference allele.
    var_allele : str or Bio.Seq.Seq
        Variant (alternate) allele.
    var_pos : int
        Genomic position of the variant.
    strand : str
        ``'+'`` or ``'-'``.
    features_info : list
        List of ``[start, end, type]`` triples.
    cds_info : list, optional
        Two-element list ``[cds_start, cds_end]`` (1-based).

    Returns
    -------
    tuple
        ``(coding_ref_seq, coding_alt_seq)``
    """
    if cds_info is None:
        cds_info = []
    alt_seq = ""
    if len(cds_info) == 2:
        start_coding_index = cds_info[0] - 1  # it should be index not pos
        stop_coding_index = cds_info[1]  # get end position of the last cds
    else:
        start_coding_index = 0
        total_len = 0
        for x in features_info:
            total_len += x[1] - x[0] + 1
        stop_coding_index = total_len

    if strand == '-':
        ref_seq = ref_seq[::-1]
        ref_allele = ref_allele.complement()
        var_allele = var_allele.complement()

    if strand == '-' and len(cds_info) == 2:
        n = len(ref_seq)
        ref_seq = ref_seq[n - stop_coding_index:n - start_coding_index]
    else:
        ref_seq = ref_seq[start_coding_index:stop_coding_index]

    nc_index = 0
    if len(ref_allele) == len(var_allele) or ref_allele[0] == var_allele[0]:
        for feature in features_info:
            if var_pos in range(feature[0], feature[1] + 1):
                var_index_in_cds = nc_index + (var_pos - feature[0])
                c = len(ref_allele)
                alt_seq = ref_seq[0:var_index_in_cds] + var_allele + ref_seq[var_index_in_cds + c::]
                if strand == '-':
                    return ref_seq[::-1], alt_seq[::-1]
                else:
                    return ref_seq, alt_seq

            nc_index += (feature[1] - feature[0] + 1)

    return ref_seq, alt_seq


def get_orfs_vcf(
    ref_seq: str,
    alt_seq: str,
    translation_table: int,
    num_orfs: int = 1,
) -> tuple[list, list]:
    """Translate the coding ref and alt sequences into ORFs.

    Parameters
    ----------
    ref_seq : str or Bio.Seq.Seq
        Coding reference sequence.
    alt_seq : str or Bio.Seq.Seq
        Coding alternate sequence.
    translation_table : int
        NCBI translation table number.
    num_orfs : int
        Number of reading frames to translate (1 or 3).

    Returns
    -------
    tuple[list, list]
        ``(ref_orfs, alt_orfs)``
    """
    ref_orfs = []
    alt_orfs = []
    for n in range(0, num_orfs):
        ref_orfs.append(ref_seq[n::].translate(translation_table))
        alt_orfs.append(alt_seq[n::].translate(translation_table))

    return ref_orfs, alt_orfs


def write_output(
    seq_id: str,
    desc: str,
    seqs: list,
    prots_fn: Any,
    seqs_filter: Optional[list] = None,
) -> None:
    """Write ORFs to a FASTA output file handle.

    Parameters
    ----------
    seq_id : str
        Sequence accession / identifier.
    desc : str
        Sequence description.
    seqs : list
        List of translated ORF sequences.
    prots_fn : file-like
        Writable file handle.
    seqs_filter : list, optional
        ORFs present in this list will be skipped (used to exclude
        alt ORFs identical to reference).
    """
    if seqs_filter is None:
        seqs_filter = []
    write_i = False
    if len(seqs) > 1:  # only add _num when multiple ORFs are generated (e.g in 3 ORF)
        write_i = True

    formatted_desc = " " + desc if desc else ""
    for i, orf in enumerate(seqs):
        if orf in seqs_filter:
            continue
        if write_i:
            prots_fn.write('>{}{}\n{}\n'.format(seq_id + "_" + str(i + 1), formatted_desc, orf))
        else:
            prots_fn.write('>{}{}\n{}\n'.format(seq_id, formatted_desc, orf))
