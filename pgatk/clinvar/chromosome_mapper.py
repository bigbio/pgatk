"""Bidirectional chromosome name mapping: RefSeq (NC_) <-> numeric <-> UCSC (chr).

Parses an NCBI assembly report to build the mapping tables.
"""
from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_VALID_CONVENTIONS = {"numeric", "refseq", "ucsc"}


class ChromosomeMapper:
    """Maps chromosome names between naming conventions.

    Supports three conventions:
    - 'numeric': e.g. '1', '2', 'X', 'Y', 'MT'
    - 'refseq': e.g. 'NC_000001.11'
    - 'ucsc': e.g. 'chr1', 'chrX', 'chrM'
    """

    def __init__(
        self,
        numeric_to_refseq: dict[str, str],
        numeric_to_ucsc: dict[str, str],
    ) -> None:
        self._num_to_ref = numeric_to_refseq
        self._num_to_ucsc = numeric_to_ucsc
        self._ref_to_num = {v: k for k, v in numeric_to_refseq.items()}
        self._ucsc_to_num = {v: k for k, v in numeric_to_ucsc.items()}

    @classmethod
    def from_assembly_report(cls, report_path: str) -> ChromosomeMapper:
        """Parse an NCBI assembly report file and build chromosome maps.

        The assembly report is a TSV with comment lines starting with '#'.
        Relevant columns (0-indexed):
        - 0: Sequence-Name (e.g. '1', 'X', 'MT')
        - 6: RefSeq-Accn (e.g. 'NC_000001.11')
        - 9: UCSC-style-name (e.g. 'chr1')
        """
        numeric_to_refseq: dict[str, str] = {}
        numeric_to_ucsc: dict[str, str] = {}

        with open(report_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                cols = line.split("\t")
                if len(cols) < 10:
                    continue
                seq_name = cols[0]
                refseq_accn = cols[6]
                ucsc_name = cols[9]
                if refseq_accn and refseq_accn != "na":
                    numeric_to_refseq[seq_name] = refseq_accn
                if ucsc_name and ucsc_name != "na":
                    numeric_to_ucsc[seq_name] = ucsc_name

        logger.info("Loaded chromosome mapping: %d entries from %s",
                     len(numeric_to_refseq), Path(report_path).name)
        return cls(numeric_to_refseq, numeric_to_ucsc)

    def _to_numeric(self, name: str) -> str | None:
        """Convert any convention to numeric. Returns None if unknown."""
        if name in self._num_to_ref or name in self._num_to_ucsc:
            return name
        if name in self._ref_to_num:
            return self._ref_to_num[name]
        if name in self._ucsc_to_num:
            return self._ucsc_to_num[name]
        return None

    def map_chrom(self, name: str, target_convention: str) -> str:
        """Map a chromosome name to the target convention.

        :param name: Input chromosome name in any convention.
        :param target_convention: One of 'numeric', 'refseq', 'ucsc'.
        :return: Mapped name, or the original name if mapping is unknown.
        :raises ValueError: If target_convention is not valid.
        """
        if target_convention not in _VALID_CONVENTIONS:
            raise ValueError(
                f"Invalid convention '{target_convention}'. "
                f"Must be one of: {', '.join(sorted(_VALID_CONVENTIONS))}"
            )
        numeric = self._to_numeric(name)
        if numeric is None:
            logger.warning("Unknown chromosome name: %s — returning as-is", name)
            return name
        if target_convention == "numeric":
            return numeric
        elif target_convention == "refseq":
            return self._num_to_ref.get(numeric, name)
        elif target_convention == "ucsc":
            return self._num_to_ucsc.get(numeric, name)
        return name
