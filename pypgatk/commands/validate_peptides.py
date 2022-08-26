import logging
from turtle import position

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.spectrumAI.validate_peptides import ValidatePeptidesService
from pypgatk.spectrumAI.get_position import GetPosition
from pypgatk.commands.utils import print_help


log = logging.getLogger(__name__)

@click.command('validate_peptides', short_help='Command to inspect MS2 spectra of single-subsititution peptide identifications')
@click.option('-c', '--config_file', help='Configuration file for the validate peptides pipeline')
@click.option('-p', '--mzml_path', help='The mzml file path.You only need to use either mzml_path or mzml_files')
@click.option('-f', '--mzml_files', help='The mzml files.Different files are separated by ",".You only need to use either mzml_path or mzml_files')
@click.option('-i', '--infile_name', help='Variant peptide PSMs table')
@click.option('-o', '--outfile_name', help='Output file for the results')
@click.option('-ion', '--ions_tolerance', help='MS2 fragment ions mass accuracy', default=0.02)
@click.option('-r', '--relative', help='relative', is_flag=True)

@click.option('-in_psms', '--input_psm_table', help='Input variant peptide PSMs table')
@click.option('-fa', '--input_fasta', help='Protein sequence used')
@click.option('-out_psms', '--output_psm_table', help='Output variant peptide PSMs table')
@click.pass_context
def validate_peptides(ctx, config_file, mzml_path, mzml_files, infile_name, outfile_name, ions_tolerance, relative, input_psm_table, input_fasta, output_psm_table):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    validate_flag = bool(infile_name and (mzml_path or mzml_files) and outfile_name)
    position_flag = bool(input_psm_table and input_fasta and output_psm_table)
    if not validate_flag and not position_flag:
        print_help()

    pipeline_arguments = {}  

    if mzml_path is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_MZML_PATH] = mzml_path
    if mzml_files is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_MZML_FILES] = mzml_files
    if infile_name is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_INFILE_NAME] = infile_name
    if outfile_name is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_OUTFILE_NAME] = outfile_name
    if ions_tolerance is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_IONS_TOLERANCE] = ions_tolerance

    if input_psm_table is not None:
        pipeline_arguments[GetPosition.CONFIG_INPUT_PSM_TABLE] = input_psm_table
    if input_fasta is not None:
        pipeline_arguments[GetPosition.CONFIG_INPUT_FASTA] = input_fasta
    if output_psm_table is not None:
        pipeline_arguments[GetPosition.CONFIG_OUTPUT_PSM_TABLE] = output_psm_table

    if validate_flag:
        validate_peptides_service = ValidatePeptidesService(config_data, pipeline_arguments)
        validate_peptides_service.validate(infile_name, outfile_name, mzml_path , mzml_files, ions_tolerance, relative)
    elif position_flag:
        get_subpos_service = GetPosition(config_data, pipeline_arguments)
        get_subpos_service.get_position(input_psm_table, input_fasta, output_psm_table)


