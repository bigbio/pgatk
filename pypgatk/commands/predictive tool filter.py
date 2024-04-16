import logging

import click

from pypgatk.toolbox.general import read_yaml_from_file
from pypgatk.commands.utils import print_help
from pypgatk.proteogenomics.blast_get_position import BlastGetPositionService