import logging
import os

import click

from pypgatk.proteomics.db.decoy_pyrat import ProteinDBService
from pypgatk.toolbox.exceptions import AppConfigException

this_dir, this_filename = os.path.split(__file__)


@click.command('generate-decoy', short_help='Create decoy protein sequences. Each protein '
                                             'is reversed and the cleavage sites switched with preceding amino acid. Peptides are checked for existence in target sequences if found'
                                             'the tool will attempt to shuffle them. James.Wright@sanger.ac.uk 2015')
@click.option('-c', '--config_file', help='Configuration file for the protein database decoy generation', default= this_dir + '/../config/protein_decoy.yaml')
@click.option('-o', '--output',      help='Output file for decoy database', default="protein-decoy.fa")
@click.option( '-i', '--input',               help='FASTA file of target proteins sequences for which to create decoys (*.fasta|*.fa)')
@click.option('-s', '--cleavage_sites',      default='KR', help='A list of amino acids at which to cleave during digestion. Default = KR')
@click.option('-a', '--anti_cleavage_sites',  help='A list of amino acids at which not to cleave if following cleavage site ie. Proline. Default = none')
@click.option('-p', '--cleavage_position',    default='c', help='Set cleavage to be c or n terminal of specified cleavage sites. Options [c, n], Default = c')
@click.option('-l', '--min_peptide_length',  type=int, help='Set minimum length of peptides to compare between target and decoy. Default = 5')
@click.option('-n', '--max_iterations',      type=int, help='Set maximum number of times to shuffle a peptide to make it non-target before failing. Default=100')
@click.option('-x', '--do_not_shuffle',       help='Turn OFF shuffling of decoy peptides that are in the target database. Default=false')
@click.option('-w', '--do_not_switch',      help='Turn OFF switching of cleavage site with preceding amino acid. Default=false')
@click.option('-d', '--decoy_prefix',       help='Set accession prefix for decoy proteins in output. Default=DECOY_')
@click.option('-t', '--temp_file',          help='Set temporary file to write decoys prior to shuffling. Default=protein-decoy.fa')
@click.option('-b', '--no_isobaric',         help='Do not make decoy peptides isobaric. Default=false')
@click.option( '-m', '--memory_save',        help='Slower but uses less memory (does not store decoy peptide list). Default=false')
@click.pass_context
def generate_database(ctx, config_file, output, input, cleavage_sites, anti_cleavage_sites, cleavage_position, min_peptide_length,
                   max_iterations, do_not_shuffle, do_not_switch, decoy_prefix, temp_file, no_isobaric, memory_save):

    if config_file is None:
        msg = "The config file for the pipeline is missing, please provide one "
        logging.error(msg)
        raise AppConfigException(msg)

    pipeline_arguments = {}

    if output is not None:
        pipeline_arguments[ProteinDBService.CONFIG_PROTEINDB_OUTPUT] = output

    if input is not None:
        pipeline_arguments[ProteinDBService.CONFIG_INPUT_FILE] = input

    if cleavage_position is not None:
        pipeline_arguments[ProteinDBService.CONFIG_CLEAVAGE_SITES] = cleavage_sites

    if cleavage_position is not None:
        pipeline_arguments[ProteinDBService.CONFIG_CLEAVAGE_POSITION] = cleavage_position

    if anti_cleavage_sites is not None:
        pipeline_arguments[ProteinDBService.CONFIG_ANTI_CLEAVAGE_SITES] = anti_cleavage_sites

    if min_peptide_length is not None:
        pipeline_arguments[ProteinDBService.CONFIG_PEPTIDE_LENGTH] = min_peptide_length

    if max_iterations is not None:
        pipeline_arguments[ProteinDBService.CONFIG_MAX_ITERATIONS] = max_iterations

    if do_not_shuffle is not None:
        pipeline_arguments[ProteinDBService.CONFIG_DO_NOT_SUFFLE] = do_not_shuffle

    if do_not_switch is not None:
        pipeline_arguments[ProteinDBService.CONFIG_DO_NOT_SWITCH] = do_not_switch

    if decoy_prefix is not None:
        pipeline_arguments[ProteinDBService.CONFIG_DECOY_PREFIX] = decoy_prefix

    if temp_file is not None:
        pipeline_arguments[ProteinDBService.CONFIG_TEMP_FILE] = temp_file

    if no_isobaric is not None:
        pipeline_arguments[ProteinDBService.CONFIG_NO_ISOBARIC] = no_isobaric

    if memory_save is not None:
        pipeline_arguments[ProteinDBService.CONFIG_MEMORY_SAVE] = memory_save

    proteindb_decoy = ProteinDBService(config_file, pipeline_arguments)
    proteindb_decoy.decoy_database()
