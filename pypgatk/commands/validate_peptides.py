import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.spectrumAI.validate_peptides import ValidatePeptidesService
from pypgatk.commands.utils import print_help


log = logging.getLogger(__name__)


@click.command('validate_peptides', short_help='Command to inspect MS2 spectra of single-subsititution peptide identifications')
@click.option('-c', '--config_file', help='Configuration file for the validate peptides pipeline')
@click.option('-m', '--mzml_path', help='The mzml file path')
@click.option('-in', '--infile_name', help='Variant peptide PSMs table')
@click.option('-o', '--outfile_name', help='Output file for the results')
@click.option('-i', '--ions_tolerance', help='MS2 fragment ions mass accuracy', default=0.02)
@click.option('-r', '--relative', help='relative', is_flag=True)
@click.pass_context
def validate_peptides(ctx, config_file, mzml_path, infile_name, outfile_name, ions_tolerance, relative):

    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if infile_name is None or mzml_path is None or outfile_name is None:
        print_help()

    pipeline_arguments = {}  

    if mzml_path is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_MZML_PATH] = mzml_path
    if infile_name is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_INFILE_NAME] = infile_name
    if outfile_name is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_OUTFILE_NAME] = outfile_name
    if ions_tolerance is not None:
        pipeline_arguments[ValidatePeptidesService.CONFIG_IONS_TOLERANCE] = ions_tolerance

    validate_peptides_service = ValidatePeptidesService(config_data, pipeline_arguments)
    validate_peptides_service.validate(infile_name,outfile_name,mzml_path ,ions_tolerance, relative)


