#!/usr/bin/env python3
"""One-shot benchmark / profiler for ``pgatk vcf-to-proteindb``.

Times one full run end-to-end. With ``--profile-out PATH``, also wraps the
run in cProfile and writes a ``.prof`` file; print the top-N hotspots with
``--print-top``. Not invoked by CI.

Basic usage:
    python scripts/benchmark_vcf_to_proteindb.py \
        --vcf    /path/to/sample.vcf \
        --fasta  /path/to/transcripts.fa \
        --gtf    /path/to/annotation.gtf \
        --output /tmp/benchmark_proteindb.fa

Profile a chr22 run and print the 30 hottest functions:
    python scripts/benchmark_vcf_to_proteindb.py \
        --vcf chr22.vcf --fasta tx.fa --gtf chr22.gtf \
        --output /tmp/chr22.fa \
        --profile-out /tmp/chr22.prof --print-top 30
"""
import argparse
import cProfile
import io
import pstats
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
    parser.add_argument(
        "--profile-out",
        help="If set, wrap the run in cProfile and write the .prof file here.",
    )
    parser.add_argument(
        "--print-top",
        type=int,
        default=0,
        help="With --profile-out, also print the top-N functions by cumulative time "
             "(default: 0, no printing). Recommended N=30 for a first look.",
    )
    args = parser.parse_args()

    config_data = load_config("ensembl_config", None)
    pipeline_arguments = {
        EnsemblDataService.PROTEIN_DB_OUTPUT: args.output,
        EnsemblDataService.ANNOTATION_FIELD_NAME: args.annotation_field_name,
    }
    svc = EnsemblDataService(config_data, pipeline_arguments)

    profiler = cProfile.Profile() if args.profile_out else None
    start = time.perf_counter()
    if profiler is not None:
        profiler.enable()
    try:
        svc.vcf_to_proteindb(args.vcf, args.fasta, args.gtf)
    finally:
        if profiler is not None:
            profiler.disable()
    elapsed = time.perf_counter() - start

    if profiler is not None:
        profiler.dump_stats(args.profile_out)

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
    if args.profile_out:
        print(f"  Profile:    {args.profile_out}")

    if args.profile_out and args.print_top > 0:
        print()
        print(f"=== TOP {args.print_top} BY CUMULATIVE TIME ===")
        buf = io.StringIO()
        stats = pstats.Stats(args.profile_out, stream=buf)
        stats.strip_dirs().sort_stats('cumulative').print_stats(args.print_top)
        print(buf.getvalue())
        print(f"=== TOP {args.print_top} BY OWN (tottime) TIME ===")
        buf = io.StringIO()
        stats = pstats.Stats(args.profile_out, stream=buf)
        stats.strip_dirs().sort_stats('tottime').print_stats(args.print_top)
        print(buf.getvalue())


if __name__ == "__main__":
    main()
