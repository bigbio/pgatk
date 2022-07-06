import base64
import json
import gzip
import requests
import os

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

        self._local_path_cosmic = self.get_configuration_default_params(variable=self.CONFIG_OUTPUT_DIRECTORY,
                                                                        default_value='./database_cosmic/')
        self._cosmic_ftp_url = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_FTP_URL,
                                                                     default_value='https://cancer.sanger.ac.uk')
        self._cosmic_user = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_FTP_USER,
                                                                  default_value='')
        self._cosmic_password = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_FTP_PASSWORD,
                                                                      default_value='')
        self._cosmic_mutation_url = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_MUTATIONS_URL,
                                                                          default_value='cosmic/file_download/GRCh38/cosmic/v94')
        self._cosmic_mutations_file = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_MUTATIONS_FILE,
                                                                            default_value='CosmicMutantExport.tsv.gz')
        self._cosmic_cellline_mutation_url = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_CELLLINE_MUTATIONS_URL,
            default_value='cosmic/file_download/GRCh38/cell_lines/v94')
        self._cosmic_cellline_mutation_file = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_CELLLINE_MUTATIONS_FILE, default_value='CosmicCLP_MutantExport.tsv.gz')
        self._cosmic_cellline_mutation_gene = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_CELLLINES_GENES_FILE, default_value='All_CellLines_Genes.fasta.gz')
        self._cosmic_cdns_file = self.get_configuration_default_params(variable=self.CONFIG_COSMIC_CDS_GENES_FILE,
                                                                       default_value='All_COSMIC_Genes.fasta.gz')

        self._cosmic_token = base64.b64encode("{}:{}".format(self._cosmic_user, self._cosmic_password)
                                              .encode()).decode('utf-8')

        self.prepare_local_cosmic_repository()

    def get_configuration_default_params(self, variable: str, default_value):
        return_value = default_value
        if variable in self.get_pipeline_parameters():
            return_value = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() \
                and self.CONFIG_COSMIC_SERVER in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER] \
                and variable in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
            self.CONFIG_COSMIC_SERVER]:
            return_value = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_COSMIC_SERVER][
                variable]
        return return_value

    def prepare_local_cosmic_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cosmic_repo()))
        check_create_folders([self.get_local_path_root_cosmic_repo()])
        self.get_logger().debug(
            "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cosmic_repo()))

    def get_local_path_root_cosmic_repo(self):
        return self._local_path_cosmic

    def download_mutation_file(self, url_file_name=None):
        """
        This function will download the mutations file from Cosmic Database.
        :return: None
        """

        mutation_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(), self._cosmic_mutations_file)
        cds_genes_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(), self._cosmic_cdns_file)

        mutation_celline_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                                      self._cosmic_cellline_mutation_file)
        cellines_genes_output_file = "{}/{}".format(self.get_local_path_root_cosmic_repo(),
                                                    self._cosmic_cellline_mutation_gene)

        server = self._cosmic_ftp_url

        cosmic_version = self._cosmic_mutation_url
        mutation_file = self._cosmic_mutations_file

        cosmic_cellline_version = self._cosmic_cellline_mutation_url
        mutation_celline_file = self._cosmic_cellline_mutation_file

        all_cds_gene_file = self._cosmic_cdns_file

        all_celllines_gene_file = self._cosmic_cellline_mutation_gene

        mutation_url = "{}/{}/{}".format(server, cosmic_version, mutation_file)
        cds_gene_url = "{}/{}/{}".format(server, cosmic_version, all_cds_gene_file)

        celllines_gene_url = "{}/{}/{}".format(server, cosmic_cellline_version, all_celllines_gene_file)
        mutation_cellline_url = "{}/{}/{}".format(server, cosmic_cellline_version, mutation_celline_file)

        if url_file_name is None:
            token = "Basic {}".format(self._cosmic_token)
            self.download_file_cosmic(mutation_url, mutation_output_file, token)
            self.download_file_cosmic(cds_gene_url, cds_genes_output_file, token)

            self.download_file_cosmic(celllines_gene_url, cellines_genes_output_file, token)
            self.download_file_cosmic(mutation_cellline_url, mutation_celline_output_file, token)

        else:
            if url_file_name is not None:
                with open(url_file_name, 'w') as url_file:
                    url_file.write("{}\t{}\n".format(mutation_url, mutation_output_file))
                    url_file.write("{}\t{}\n".format(cds_gene_url, cds_genes_output_file))
                    url_file.write("{}\t{}\n".format(celllines_gene_url, cellines_genes_output_file))
                    url_file.write("{}\t{}\n".format(mutation_cellline_url, mutation_celline_output_file))

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
                if local_file.endswith('.gz'):
                    extracted_file = local_file.replace('.gz', '')
                    with open(extracted_file, 'w') as outfile:
                        try:
                            outfile.write(gzip.decompress(open(local_file, 'rb').read()).decode('utf-8'))
                        except UnicodeDecodeError:
                            outfile.write(gzip.decompress(open(local_file, 'rb').read()).decode('ISO-8859–1'))
                        os.remove(local_file)
                        local_file = extracted_file
                        msg = "Extracted file '{}'".format(local_file)
                        self.get_logger().debug(msg)
        else:
            msg = "Error downloading the COSMIC data, error code {} , error message '{}'".format(response.status_code,
                                                                                                 local_file)
            self.get_logger().debug(msg)
            raise AppConfigException(msg)
