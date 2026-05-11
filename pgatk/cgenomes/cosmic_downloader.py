import base64
import os
import tarfile

import requests
from tqdm import tqdm

from pgatk.toolbox.exceptions import AppConfigException
from pgatk.toolbox.general import ParameterConfiguration, check_create_folders


class CosmicDownloadService(ParameterConfiguration):
    CONFIG_KEY_DATA_DOWNLOADER = 'cosmic_data'
    CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    CONFIG_COSMIC_SERVER = 'cosmic_server'
    CONFIG_COSMIC_FTP_URL = 'cosmic_ftp'
    CONFIG_API_ENDPOINT = 'api_endpoint'
    CONFIG_BUCKET = 'bucket'
    CONFIG_COSMIC_FTP_USER = "cosmic_user"
    CONFIG_COSMIC_FTP_PASSWORD = "cosmic_password"
    CONFIG_PRODUCTS = "products"

    # Built-in fallback list, used only if `products` is missing from config.
    # Real default lives in pgatk/config/cosmic_config.yaml and is the source of truth.
    _DEFAULT_PRODUCTS = (
        'grch38/cosmic/v103/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38.tar',
        'grch38/cosmic/v103/Cosmic_CompleteTargetedScreensMutant_Tsv_v103_GRCh38.tar',
        'grch38/cosmic/v103/Cosmic_Genes_Fasta_v103_GRCh38.tar',
        'grch38/cosmic/v103/Cosmic_Transcripts_Tsv_v103_GRCh38.tar',
        'grch38/cell_lines/v103/CellLinesProject_GenomeScreensMutant_Tsv_v103_GRCh38.tar',
        'grch38/cell_lines/v103/CellLinesProject_CompleteCNA_Tsv_v103_GRCh38.tar',
    )

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CosmicDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_file, pipeline_arguments)

        self._local_path_cosmic = self.get_configuration_default_params(
            variable=self.CONFIG_OUTPUT_DIRECTORY, default_value='./database_cosmic/')
        self._cosmic_ftp_url = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_FTP_URL, default_value='https://cancer.sanger.ac.uk')
        self._api_endpoint = self.get_configuration_default_params(
            variable=self.CONFIG_API_ENDPOINT, default_value='api/mono/products/v1/downloads/scripted')
        self._bucket = self.get_configuration_default_params(
            variable=self.CONFIG_BUCKET, default_value='downloads')
        self._cosmic_user = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_FTP_USER, default_value='')
        self._cosmic_password = self.get_configuration_default_params(
            variable=self.CONFIG_COSMIC_FTP_PASSWORD, default_value='')
        self._products = self.get_configuration_default_params(
            variable=self.CONFIG_PRODUCTS, default_value=list(self._DEFAULT_PRODUCTS))

        self._cosmic_token = base64.b64encode(
            "{}:{}".format(self._cosmic_user, self._cosmic_password).encode()
        ).decode('utf-8')

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

    def build_api_url(self, product_path):
        """Build the COSMIC scripted-download API URL for a single product."""
        return "{server}/{endpoint}?path={path}&bucket={bucket}".format(
            server=self._cosmic_ftp_url.rstrip('/'),
            endpoint=self._api_endpoint.strip('/'),
            path=product_path,
            bucket=self._bucket,
        )

    def download_mutation_file(self, url_file_name=None):
        """
        Download every product listed in `self._products` via COSMIC's v103+
        scripted-download API. Each product is a `path=` value such as
        `grch38/cosmic/v103/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38.tar`.

        If url_file_name is provided, write the API URLs to that file instead of downloading.
        """

        products = list(self._products)
        if not products:
            self.get_logger().warning("No COSMIC products configured; nothing to download.")
            return

        token = "Basic {}".format(self._cosmic_token)
        output_dir = self.get_local_path_root_cosmic_repo()

        if url_file_name is not None:
            with open(url_file_name, 'w', encoding='utf-8') as url_file:
                for path in products:
                    api_url = self.build_api_url(path)
                    output_file = "{}/{}".format(output_dir, os.path.basename(path))
                    url_file.write("{}\t{}\n".format(api_url, output_file))
            return

        for path in products:
            api_url = self.build_api_url(path)
            output_file = "{}/{}".format(output_dir, os.path.basename(path))
            self.download_file_cosmic(api_url, output_file, token)

    def download_file_cosmic(self, api_url, local_file, token):
        """
        Two-step download from COSMIC's scripted-download API:
        1) GET api_url with Basic auth -> returns JSON {"url": "<s3-presigned-url>"}
        2) GET the presigned URL (no auth, 1-hour TTL) -> stream to local_file

        If local_file ends in .tar, extract its contents to the same directory
        and delete the archive.
        """
        api_response = requests.get(api_url, headers={'Authorization': token}, timeout=30)
        if api_response.status_code != 200:
            msg = ("COSMIC API request failed: HTTP {} for {}. Body: {!r}"
                   .format(api_response.status_code, api_url, api_response.text[:200]))
            self.get_logger().error(msg)
            raise AppConfigException(msg)

        try:
            payload = api_response.json()
            download_url = payload['url']
        except (ValueError, KeyError) as exc:
            msg = ("COSMIC API returned 200 but the body is not a JSON object with a 'url' "
                   "key (got: {!r}). This usually means the `path=` value is wrong for the "
                   "release/assembly, or your account lacks access. URL was: {}. Error: {}"
                   .format(api_response.text[:200], api_url, exc))
            self.get_logger().error(msg)
            raise AppConfigException(msg)

        self.get_logger().debug("Downloading file from signed URL '{}'".format(download_url))
        data_response = requests.get(download_url, stream=True, timeout=30)
        if data_response.status_code != 200:
            msg = ("COSMIC S3 download failed: HTTP {} for {}"
                   .format(data_response.status_code, download_url))
            self.get_logger().error(msg)
            raise AppConfigException(msg)

        total_size = int(data_response.headers.get('content-length') or 0)
        chunk_size = 1024 * 1024
        progress = tqdm(
            total=total_size if total_size > 0 else None,
            unit='B', unit_scale=True, unit_divisor=1024,
            desc=os.path.basename(local_file),
            leave=True,
        )
        try:
            with open(local_file, 'wb') as f:
                for chunk in data_response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        progress.update(len(chunk))
        finally:
            progress.close()
        self.get_logger().debug("Download finished: '{}'".format(local_file))

        if local_file.endswith('.tar'):
            extract_dir = os.path.dirname(local_file) or '.'
            print("Extracting {}...".format(os.path.basename(local_file)), flush=True)
            with tarfile.open(local_file, 'r') as tar:
                tar.extractall(path=extract_dir)
            os.remove(local_file)
            self.get_logger().debug("Extracted archive into '{}'".format(extract_dir))
