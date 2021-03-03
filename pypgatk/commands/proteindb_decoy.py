import logging
import os

import click

from pypgatk.proteomics.db.protein_database_decoy import ProteinDBDecoyService
from pypgatk.proteomics.models import PYGPATK_ENZYMES
from pypgatk.toolbox.exceptions import AppConfigException

this_dir, this_filename = os.path.split(__file__)


@click.command('generate-decoy', short_help='Create decoy protein sequences using multiple methods DecoyPYrat, Reverse/Shuffled Proteins.')
@click.option('-c', '--config_file', help='Configuration file for the protein database decoy generation',
              default=this_dir + '/../config/protein_decoy.yaml')
@click.option('-out', '--output_database', help='Output file for decoy database', default="protein-decoy.fa")
@click.option('-in', '--input_database', help='FASTA file of target proteins sequences for which to create decoys (*.fasta|*.fa)')
@click.option('-m', '--method', type=click.Choice(['protein-reverse', 'protein-shuffle', 'decoypyrat', 'pgdbdeep'], case_sensitive=False),
              help='The method that would be used to generate the decoys:\n'
                   'protein-reverse: reverse protein sequences\n'
                   'protein-shuffle: shuffle protein sequences\n'
                   'decoypyrat: method developed for proteogenomics that shuffle redundant peptides in decoy databases\n'
                   'pgdbdeep: method developed for proteogenomics developed by pypgatk',
              default = 'protein-reverse')
@click.option('-d', '--decoy_prefix', help='Set accession prefix for decoy proteins in output. Default=DECOY_', default='DECOY_')
@click.option('-e', '--enzyme', default='Trypsin', help='Enzyme used for clevage the protein sequence (Default: Trypsin)', type = click.Choice(PYGPATK_ENZYMES.enzymes.keys()))
@click.option('--cleavage_position', default='c', type=click.Choice(['c','n']),help='Set cleavage to be c or n terminal of specified cleavage sites. Options [c, n], Default = c')
@click.option('-s', '--max_missed_cleavages', default = 0, type = int, help = 'Number of allowed missed cleavages in the protein sequence when digestion is performed')
@click.option('--min_peptide_length', type=int, help='Set minimum length of peptides (Default = 5)', default = 6)
@click.option('--max_peptide_length', type=int, help='Set maximum length of peptides (Default = 100)', default = 1000)
@click.option('--max_iterations', type=int, help='Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100', default = 100)
@click.option('--do_not_shuffle', help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false', is_flag = True, default = False)
@click.option('--do_not_switch', help='Turn OFF switching of cleavage site with preceding amino acid. Default=false', is_flag = True, default = False)
@click.option('--temp_file', help='Set temporary file to write decoys prior to shuffling. Default=tmp.fa', default = 'tmp.fa')
@click.option('--no_isobaric', help='Do not make decoy peptides isobaric. Default=false', is_flag = True, default = False)
@click.option('--keep_target_hits', help='Keep peptides duplicate in target and decoy databases', is_flag = True, default = False)
@click.option('--memory_save', help='Slower but uses less memory (does not store decoy peptide list). Default=false', is_flag = True, default = False)
@click.pass_context
def generate_database(ctx, config_file: str, output_database: str, input_database: str, method: str,
                      decoy_prefix: str , enzyme: str, cleavage_position: str,
                      max_missed_cleavages: int, min_peptide_length: int, max_peptide_length : int,
                      max_iterations: int, do_not_shuffle: bool , do_not_switch: bool, temp_file: str,
                      no_isobaric: bool, keep_target_hits: bool, memory_save: bool):
  if config_file is None:
    msg = "The config file for the pipeline is missing, please provide one "
    logging.error(msg)
    raise AppConfigException(msg)

  pipeline_arguments = {}

  if output_database is None or input_database is None or method is None:
    help(ctx)

  if output_database is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_PROTEINDB_OUTPUT] = output_database

  if input_database is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_INPUT_FILE] = input_database

  if method is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_DECOY_METHOD] = method

  if enzyme is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_ENZYME] = enzyme

  if cleavage_position is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_CLEAVAGE_POSITION] = cleavage_position

  if max_missed_cleavages is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_MAX_MISSED_CLEAVAGES] = max_missed_cleavages

  if min_peptide_length is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_MIN_PEPTIDE_LENGTH] = min_peptide_length

  if max_peptide_length is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_PEPTIDE_LENGTH_MAX] = max_peptide_length

  if max_iterations is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_MAX_ITERATIONS] = max_iterations

  if do_not_shuffle is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_DO_NOT_SUFFLE] = do_not_shuffle

  if do_not_switch is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_DO_NOT_SWITCH] = do_not_switch

  if decoy_prefix is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_DECOY_PREFIX] = decoy_prefix

  if temp_file is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_TEMP_FILE] = temp_file

  if no_isobaric is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_NO_ISOBARIC] = no_isobaric

  if keep_target_hits is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_KEEP_TARGET_HITS] = keep_target_hits

  if memory_save is not None:
    pipeline_arguments[ProteinDBDecoyService.CONFIG_MEMORY_SAVE] = memory_save

  proteindb_decoy = ProteinDBDecoyService(config_file, pipeline_arguments)
  proteindb_decoy.decoy_database()
  proteindb_decoy.print_target_decoy_composition()
