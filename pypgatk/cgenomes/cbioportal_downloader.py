import csv
from concurrent.futures import as_completed
from concurrent.futures.thread import ThreadPoolExecutor

from requests import get

from pypgatk.toolbox.exceptions import AppException
from pypgatk.toolbox.general import ParameterConfiguration, check_create_folders, download_file, clear_cache
from pypgatk.toolbox.rest import call_api, call_api_raw


class CbioPortalDownloadService(ParameterConfiguration):
    CONFIG_KEY_DATA_DOWNLOADER = 'cbioportal_data_downloader'
    CONFIG_KEY_CBIOPORTAL_DOWNLOAD_URL = 'cbioportal_download_url'
    CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    CONFIG_CBIOPORTAL_API = 'cbioportal_api'
    CONFIG_CBIOPORTAL_API_SERVER = 'base_url'
    CONFIG_CBIOPORTAL_API_CANCER_STUDIES = "cancer_studies"
    CONFIG_LIST_STUDIES = "list_studies"
    CONFIG_MULTITHREADING = "multithreading"

    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CbioPortalDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_file,
                                                        pipeline_arguments)

        self._cbioportal_studies = []
        if self.CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_cbioportal = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_DIRECTORY]
        else:
            self._local_path_cbioportal = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_OUTPUT_DIRECTORY]

        self._list_studies = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_LIST_STUDIES]
        if self.CONFIG_LIST_STUDIES in self.get_pipeline_parameters():
            self._list_studies = self.get_pipeline_parameters()[self.CONFIG_LIST_STUDIES]

        self._multithreading = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
            self.CONFIG_MULTITHREADING]
        if self.CONFIG_MULTITHREADING in self.get_pipeline_parameters():
            self._multithreading = self.get_pipeline_parameters()[self.CONFIG_MULTITHREADING]

        self.prepare_local_cbioportal_repository()

    def prepare_local_cbioportal_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cbioportal_repo()))
        check_create_folders([self.get_local_path_root_cbioportal_repo()])
        self.get_logger().debug(
            "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cbioportal_repo()))

    def get_local_path_root_cbioportal_repo(self):
        return self._local_path_cbioportal

    def get_cancer_studies(self):
        """
        This method will print the list of all cancer studies for the user.
        :return:
        """
        server = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
            self.CONFIG_CBIOPORTAL_API_SERVER]
        endpoint = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
            self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES]
        self._cbioportal_studies = call_api_raw(server + "?" + endpoint).text
        return self._cbioportal_studies

    def download_study(self, download_study):
        """
        This function will download a study from cBioPortal using the study ID
        :param download_study: Study to be download, if the study is empty or None, all the studies will be
        downloaded.
        :return: None
        """

        clear_cache()

        if self._cbioportal_studies is None or len(self._cbioportal_studies) == 0:
            self.get_cancer_studies()

        if 'all' not in download_study:
            if not self.check_study_identifier(download_study):
                msg = "The following study accession '{}' is not present in cBioPortal Studies".format(download_study)
                self.get_logger().debug(msg)
                raise AppException(msg)
            else:
                self.download_one_study(download_study)
        else:
            csv_reader = csv.reader(self._cbioportal_studies.splitlines(), delimiter="\t")
            line_count = 0
            if self._multithreading:
                processes = []
                with ThreadPoolExecutor(max_workers=10, thread_name_prefix='Thread-Download') as executor:
                    for row in csv_reader:
                        if line_count != 0:
                            processes.append(executor.submit(self.download_one_study, row[0]))
                        line_count = line_count + 1
                for task in as_completed(processes):
                    print(task.result())
            else:
                for row in csv_reader:
                    if line_count != 0:
                        self.download_one_study(row[0])
                    line_count = line_count + 1

    def download_one_study(self, download_study):
        file_name = '{}.tar.gz'.format(download_study)
        file_url = '{}/{}'.format(
            self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_KEY_CBIOPORTAL_DOWNLOAD_URL],
            file_name)
        file_name = download_file(file_url, self.get_local_path_root_cbioportal_repo() + '/' + file_name, self.get_logger())
        if file_name is not None:
            msg = "The following study '{}' has been downloaded. ".format(download_study)
        else:
            msg = "The following study '{}' hasn't been downloaded. ".format(download_study)
        self.get_logger().debug(msg)
        return file_name

    def check_study_identifier(self, download_study):
        return download_study in self._cbioportal_studies
