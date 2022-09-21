import logging
from turtle import position

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.spectrumAI.validate_peptides import ValidatePeptidesService
from pypgatk.commands.utils import print_help


log = logging.getLogger(__name__)

@click.command('validate_peptides', short_help='Command to inspect MS2 spectra of single-subsititution peptide identifications')
@click.option('-c', '--config_file', help='Configuration file for the validate peptides pipeline')
@click.option('-p', '--mzml_path', help='The mzml file path.You only need to use either mzml_path or mzml_files')
@click.option('-f', '--mzml_files', help='The mzml files.Different files are separated by ",".You only need to use either mzml_path or mzml_files')
@click.option('-i', '--infile_name', help='Variant peptide PSMs table')
@click.option('-o', '--outfile_name', help='Output file for the results')
@click.option('-ion', '--ions_tolerance', help='MS2 fragment ions mass accuracy')
@click.option('-r', '--relative', help='relative', is_flag=True)
@click.option('-msgf', '--msgf', help='If it is the standard format of MSGF output, please turn on this switch, otherwise it defaults to mzTab format', is_flag=True)
@click.option('-in_psms', '--input_psm_table', help='Input variant peptide PSMs table')
@click.option('-fa', '--input_fasta', help='Protein sequence used')
@click.option('-out_psms', '--output_psm_table', help='Output variant peptide PSMs table')
@click.pass_context
def validate_peptides(ctx, config_file, mzml_path, mzml_files, infile_name, outfile_name, ions_tolerance, relative, input_psm_table, input_fasta, output_psm_table ,msgf):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    validate_flag = bool(infile_name and (mzml_path or mzml_files) and outfile_name)
    position_flag = bool(input_psm_table and input_fasta and output_psm_table)
    if not validate_flag and not position_flag:
        print_help()

    pipeline_arguments = {}  

    if ions_tolerance is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_IONS_TOLERANCE] = ions_tolerance  
    if relative is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_RELATIVE] = relative
    if msgf is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_MSGF] = msgf
    
    validate_peptides_service = ValidatePeptidesService(config_data, pipeline_arguments)
    if validate_flag:
        validate_peptides_service.validate(infile_name, outfile_name, mzml_path , mzml_files)
    elif position_flag:
        validate_peptides_service.get_position(input_psm_table, input_fasta, output_psm_table)


