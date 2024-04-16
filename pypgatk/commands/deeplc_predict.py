import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.commands.utils import print_help
from pypgatk.proteogenomics.deeplc_predict import DeepLCPredictService

log = logging.getLogger(__name__)

@click.command('deeplc_predict', short_help='Use deeplc prediction to filter.')
@click.option('-c', '--config_file', help='Configuration file for the fdr peptides pipeline.')
@click.option('-i', '--input_to_deeplc', help='The file name of the input PSM table to deeplc.')
@click.option('-o', '--output', help='DeepLC filtered results file.')
@click.option('-tl', '--transfer_learning', help='Use transfer learning as calibration method. 1 indicates enabled, 0 indicates disabled. Default is 0.')
@click.option('-r', '--filtration_ratio', help='Default is 0.01, filtering out 1% of the data.')
@click.option('-d', '--redundancy_removal_strategy', help='The optional parameters include the "median", "mean", "bestRT". Default is median.')
@click.option('-n', '--number_of_processes', help='Used to specify the number of processes. Default is 40.')

@click.pass_context
def deeplc_predict(ctx, config_file, input_to_deeplc, output, transfer_learning,filtration_ratio ,redundancy_removal_strategy, number_of_processes):  
    config_data = None
    if config_file is not None:
        config_data = read_yaml_from_file(config_file)

    if input_to_deeplc is None or output is None:
        print_help()
    
    pipeline_arguments = {}

    if transfer_learning is not None:
        pipeline_arguments[DeepLCPredictService.CONFIG_TRANSFER_LEARNING] = transfer_learning
    if redundancy_removal_strategy is not None:
        pipeline_arguments[DeepLCPredictService.CONFIG_REDUNDANCY_REMOVAL_STRATEGY] = redundancy_removal_strategy
    if number_of_processes is not None:
        pipeline_arguments[DeepLCPredictService.CONFIG_NUMBER_OF_PROCESSES] = number_of_processes
    if filtration_ratio is not None:
        pipeline_arguments[DeepLCPredictService.CONFIG_FILTRATION_RATIO] = filtration_ratio

    deeplc_predict_service = DeepLCPredictService(config_data, pipeline_arguments)
    deeplc_predict_service.predict(input_to_deeplc, output)