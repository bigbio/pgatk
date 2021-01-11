#!/usr/bin/env python3

"""
This is the main tool that give access to all commands and options provided by the py-pgatk

@author ypriverol

"""
import click

from pypgatk.commands import ensembl_downloader as ensembl_downloader_cmd
from pypgatk.commands import ensembl_database as ensembl_database_cmd
from pypgatk.commands import cosmic_downloader as cosmic_downloader_cmd
from pypgatk.commands import cbioportal_downloader as cbioportal_downloader_cmd
from pypgatk.commands import cosmic_to_proteindb as cosmic_to_proteindb_cmd
from pypgatk.commands import cbioportal_to_proteindb as cbioportal_to_proteindb_cmd
from pypgatk.commands import threeframe_translation as threeframe_translation_cmd
from pypgatk.commands import vcf_to_proteindb as vcf_to_proteindb_cmd
from pypgatk.commands import dnaseq_to_proteindb as dnase_to_proteindb_cmd
from pypgatk.commands import proteindb_decoy as proteindb_decoy_cmd

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# Cli returns command line requests
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """This is the main tool that give access to all commands and options provided by the pypgatk"""


cli.add_command(ensembl_downloader_cmd.ensembl_downloader)
cli.add_command(ensembl_database_cmd.ensembl_check)
cli.add_command(cbioportal_downloader_cmd.cbioportal_downloader)
cli.add_command(cosmic_downloader_cmd.cosmic_downloader)
cli.add_command(cosmic_to_proteindb_cmd.cosmic_to_proteindb)
cli.add_command(cbioportal_to_proteindb_cmd.cbioportal_to_proteindb)
cli.add_command(threeframe_translation_cmd.threeframe_translation)
cli.add_command(vcf_to_proteindb_cmd.vcf_to_proteindb)
cli.add_command(dnase_to_proteindb_cmd.dnaseq_to_proteindb)
cli.add_command(proteindb_decoy_cmd.generate_database)


def main():
    cli()


if __name__ == "__main__":
    main()
