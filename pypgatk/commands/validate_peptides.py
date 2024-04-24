import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.proteogenomics.spectrumai import SpectrumAIService
from pypgatk.commands.utils import print_help

log = logging.getLogger(__name__)


@click.command('spectrumai',
               short_help='Command to inspect MS2 spectra of single-subsititution peptide identifications')
@click.option('-c', '--config_file', help='Configuration file for the validate peptides pipeline')
@click.option('-p', '--mzml_path', help='The mzml file path.You only need to use either mzml_path or mzml_files')
@click.option('-f', '--mzml_files', help='The mzml files. Different files are separated by ","')
@click.option('-i', '--infile_name', help='Variant peptide PSMs table')
@click.option('-o', '--outfile_name', help='Output file for the results')
@click.option('-ion', '--ions_tolerance', help='MS2 fragment ions mass accuracy')
@click.option('-n', '--number_of_processes', help='Used to specify the number of processes. Default is 40.')
@click.option('-r', '--relative', help='When using ppm as ions_tolerance (not Da), it needs to be turned on',
              is_flag=True)
@click.option('-m', '--mztab',
              help='If the tsv was obtained from mzTab, please enable this option. Default to tsv obtained from parquet',
              is_flag=True)
@click.pass_context
def spectrumai(ctx, config_file, mzml_path, mzml_files, infile_name, outfile_name, ions_tolerance,
               number_of_processes, relative, mztab):
    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    validate_flag = bool(infile_name and (mzml_path or mzml_files) and outfile_name)
    if not validate_flag:
        print_help()

    pipeline_arguments = {}

    if mzml_path is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_MZML_PATH] = mzml_path
    if mzml_files is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_MZML_FILES] = mzml_files
    if ions_tolerance is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_IONS_TOLERANCE] = ions_tolerance
    if number_of_processes is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_NUMBER_OF_PROCESSES] = number_of_processes
    if relative is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_RELATIVE] = relative
    if mztab is not None:
        pipeline_arguments[SpectrumAIService.CONFIG_MZTAB] = mztab

    validate_peptides_service = SpectrumAIService(config_data, pipeline_arguments)
    if validate_flag:
        validate_peptides_service.validate(infile_name, outfile_name)
