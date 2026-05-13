"""Equivalence test: vcf_to_proteindb with workers=2 produces the same sequence
content as workers=1 on a small fixture."""
import os
from pathlib import Path

from click.testing import CliRunner

from pgatk.cli import cli


def _sequence_set(fasta_path) -> set:
    seqs = []
    current = []
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current))
                    current = []
            else:
                current.append(line.strip())
        if current:
            seqs.append(''.join(current))
    return set(seqs)


def test_parallel_matches_sequential(tmp_path):
    # pgatk package root: .../pgatk/pgatk/
    pkg_root = Path(__file__).resolve().parents[1]
    testdata = pkg_root / 'testdata'
    config = pkg_root / 'config' / 'ensembl_config.yaml'
    out_seq = tmp_path / 'seq.fa'
    out_par = tmp_path / 'par.fa'

    common_args = [
        'vcf-to-proteindb',
        '--config_file', str(config),
        '--vcf', str(testdata / 'test.vcf'),
        '--input_fasta', str(testdata / 'test.fa'),
        '--gene_annotations_gtf', str(testdata / 'test.gtf'),
        '--protein_prefix', 'ensvar',
        '--af_field', 'MAF',
        '--annotation_field_name', 'CSQ',
        '--biotype_str', 'feature_type',
        '--include_biotypes', 'mRNA,ncRNA',
    ]

    runner = CliRunner()
    orig_dir = os.getcwd()
    try:
        # Run from the package root so the .db file is created next to the .gtf.
        os.chdir(str(pkg_root))

        r1 = runner.invoke(cli, common_args + ['--output_proteindb', str(out_seq)])
        assert r1.exit_code == 0, r1.output + (str(r1.exception) if r1.exception else '')

        r2 = runner.invoke(cli, common_args + ['--output_proteindb', str(out_par), '--workers', '2'])
        assert r2.exit_code == 0, r2.output + (str(r2.exception) if r2.exception else '')
    finally:
        os.chdir(orig_dir)

    # FASTA header order may differ across runs (per-worker ordering); compare
    # sequence content as a set.
    assert _sequence_set(out_seq) == _sequence_set(out_par)
