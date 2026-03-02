"""
Central config registry: maps config keys to bundled default YAML files.
Enables auto-loading of defaults when user doesn't pass -c.
"""
import os

from pgatk.toolbox.general import read_yaml_from_file

CONFIG_DIR = os.path.dirname(__file__)

COMMAND_CONFIGS = {
    "ensembl_downloader": "ensembl_downloader_config.yaml",
    "ensembl_config": "ensembl_config.yaml",
    "cosmic": "cosmic_config.yaml",
    "cbioportal": "cbioportal_config.yaml",
    "protein_decoy": "protein_decoy.yaml",
    "clinvar": "clinvar_config.yaml",
}


def load_config(config_key, user_config_file=None):
    """
    Load configuration for a command.

    If user_config_file is provided, loads that file.
    Otherwise, loads the bundled default config for the given key.

    :param config_key: key in COMMAND_CONFIGS
    :param user_config_file: optional path to user-provided config file
    :return: parsed YAML config dict, or empty dict if no config found
    """
    if user_config_file is not None:
        return read_yaml_from_file(user_config_file)

    default_file = COMMAND_CONFIGS.get(config_key)
    if default_file is None:
        return {}
    return read_yaml_from_file(os.path.join(CONFIG_DIR, default_file))
