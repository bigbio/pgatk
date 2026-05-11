#!/usr/bin/env python3
"""One-shot benchmark for ``pgatk vcf-to-proteindb``.

Times one full run end-to-end. Compare two runs of this script (one on a
pre-fix build, one on a post-fix build) to measure the speedup achieved
by issue #99. Not invoked by CI.

Usage:
    python scripts/benchmark_vcf_to_proteindb.py \
        --vcf    /path/to/sample.vcf \
        --fasta  /path/to/transcripts.fa \
        --gtf    /path/to/annotation.gtf \
        --output /tmp/benchmark_proteindb.fa
"""
import argparse
import time
from pathlib import Path

from pgatk.ensembl.ensembl import EnsemblDataService
from pgatk.config.registry import load_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="Input VCF (annotated or raw)")
    parser.add_argument("--fasta", required=True, help="Transcript FASTA matching the GTF")
    parser.add_argument("--gtf", required=True, help="Gene annotations GTF")
    parser.add_argument("--output", required=True, help="Output protein DB FASTA")
    parser.add_argument(
        "--annotation-field-name",
        default="CSQ",
        help="VCF INFO field that carries per-transcript annotation (default: CSQ). "
        "Use empty string to force re-annotation via bedtools intersect.",
    )
    args = parser.parse_args()

    config_data = load_config("ensembl_config", None)
    pipeline_arguments = {
        EnsemblDataService.PROTEIN_DB_OUTPUT: args.output,
        EnsemblDataService.ANNOTATION_FIELD_NAME: args.annotation_field_name,
    }
    svc = EnsemblDataService(config_data, pipeline_arguments)

    start = time.perf_counter()
    svc.vcf_to_proteindb(args.vcf, args.fasta, args.gtf)
    elapsed = time.perf_counter() - start

    output_path = Path(args.output)
    output_size = output_path.stat().st_size if output_path.exists() else 0
    seq_count = 0
    if output_path.exists():
        with open(output_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    seq_count += 1

    print()
    print("=== BENCHMARK RESULT ===")
    print(f"  VCF:        {args.vcf}")
    print(f"  Elapsed:    {elapsed:.2f} s  ({elapsed / 60:.2f} min)")
    print(f"  Output:     {args.output}  ({output_size} bytes, {seq_count} sequences)")


if __name__ == "__main__":
    main()
