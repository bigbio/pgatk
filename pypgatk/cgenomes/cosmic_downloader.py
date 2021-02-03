import base64
import json

import requests

from pypgatk.toolbox.exceptions import AppConfigException
from pypgatk.toolbox.general import ParameterConfiguration, check_create_folders


class CosmicDownloadService(ParameterConfiguration):
  CONFIG_KEY_DATA_DOWNLOADER = 'cosmic_data'
  CONFIG_OUTPUT_DIRECTORY = 'output_directory'
  CONFIG_COSMIC_SERVER = 'cosmic_server'
  CONFIG_COSMIC_FTP_URL = 'cosmic_ftp'
  CONFIG_COSMIC_FTP_USER = "cosmic_user"
  CONFIG_COSMIC_FTP_PASSWORD = "cosmic_password"
  CONFIG_COSMIC_MUTATIONS_URL = "mutations_url"
  CONFIG_COSMIC_MUTATIONS_FILE = "mutations_file"
  CONFIG_COSMIC_CELLLINE_MUTATIONS_URL = "mutations_cellline_url"
  CONFIG_COSMIC_CELLLINE_MUTATIONS_FILE = "mutations_cellline_file"
  CONFIG_COSMIC_CDS_GENES_FILE = "all_cds_genes_file"
  CONFIG_COSMIC_CELLLINES_GENES_FILE = "all_celllines_genes_file"

  def __init__(self, config_file, pipeline_arguments):
    """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
    super(CosmicDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_file, pipeline_arguments)

    if self.CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
      self._local_path_cosmic = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_DIRECTORY]
    else:
      self._local_path_cosmic = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
        self.CONFIG_OUTPUT_DIRECTORY]

    self._cosmic_token = base64.b64encode("{}:{}".format(self.get_pipeline_parameters()[self.CONFIG_COSMIC_FTP_USER],
                                                         self.get_pipeline_parameters()[
                                                           self.CONFIG_COSMIC_FTP_PASSWORD])
                                          .encode()).decode('utf-8')

    self.prepare_local_cosmic_repository()

  def prepare_local_cosmic_repository(self):
    self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
      self.get_local_path_root_cosmic_repo()))
    check_create_folders([self.get_local_path_root_cosmic_repo()])
    self.get_logger().debug(
      "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cosmic_repo()))

  def get_local_path_root_cosmic_repo(self):
    return self._local_path_cosmic

  def download_mutation_file(self):
    """
        This function will download the mutations file from Cosmic Database.
        :return: None
        """

    mutation_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                          self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                                            self.CONFIG_COSMIC_SERVER][self.CONFIG_COSMIC_MUTATIONS_FILE])
    cds_genes_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                           self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                                             self.CONFIG_COSMIC_SERVER][self.CONFIG_COSMIC_CDS_GENES_FILE])

    mutation_celline_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                                  self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                                                    self.CONFIG_COSMIC_SERVER][
                                                    self.CONFIG_COSMIC_CELLLINE_MUTATIONS_FILE])
    cellines_genes_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                                self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                                                  self.CONFIG_COSMIC_SERVER][self.CONFIG_COSMIC_CELLLINES_GENES_FILE])

    server = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_FTP_URL]

    cosmic_version = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_MUTATIONS_URL]
    mutation_file = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_MUTATIONS_FILE]

    cosmic_cellline_version = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_CELLLINE_MUTATIONS_URL]
    mutation_celline_file = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_CELLLINE_MUTATIONS_FILE]

    all_cds_gene_file = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_CDS_GENES_FILE]

    all_celllines_gene_file = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
      self.CONFIG_COSMIC_CELLLINES_GENES_FILE]

    mutation_url = "{}/{}/{}".format(server, cosmic_version, mutation_file)
    cds_gene_url = "{}/{}/{}".format(server, cosmic_version, all_cds_gene_file)

    celllines_gene_url = "{}/{}/{}".format(server, cosmic_cellline_version, all_celllines_gene_file)
    mutation_cellline_url = "{}/{}/{}".format(server, cosmic_cellline_version, mutation_celline_file)

    token = "Basic {}".format(self._cosmic_token)
    self.download_file_cosmic(mutation_url, mutation_output_file, token)
    self.download_file_cosmic(cds_gene_url, cds_genes_output_file, token)

    self.download_file_cosmic(celllines_gene_url, cellines_genes_output_file, token)
    self.download_file_cosmic(mutation_cellline_url, mutation_celline_output_file, token)

  def download_file_cosmic(self, url, local_file, token):
    """
        Download file from cosmic repository using requests
        :param url: url of the file to be download
        :param local_file: local file
        :param token: token to be used
        :return:
        """

    response = requests.get(url, stream=True, headers={'Authorization': token})
    if response.status_code == 200:
      url = json.loads(response.text)['url']
      msg = "Downloading file from url '{}'".format(url)
      self.get_logger().debug(msg)

      response = requests.get(url, stream=True)
      if response.status_code == 200:
        with open(local_file, 'wb') as f:
          f.write(response.content)
          msg = "Download Finish for file '{}'".format(local_file)
          self.get_logger().debug(msg)
    else:
      msg = "Error downloading the COSMIC data, error code {} , error message '{}'".format(response.status_code,
                                                                                           local_file)
      self.get_logger().debug(msg)
      raise AppConfigException(msg)
