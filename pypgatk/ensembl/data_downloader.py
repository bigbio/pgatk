"""
This module implements en Ensembl data grabber for a given Ensembl Service instance.

Some use cases for this module:
    1. Given a species ID, download its protein sequence data, with the option of decompressing it or not.
    2. Given a species ID, collect its GTF data, with the option of decompressing it or not.
"""

# App imports
from json import loads
from pypgatk.toolbox.general import ParameterConfiguration, check_create_folders, download_file
from pypgatk.toolbox.rest import call_api


class EnsemblDataDownloadService(ParameterConfiguration):
    """
  This Service is in charge of grabbing data (download) from Ensembl to a local repository
  """

    CONFIG_KEY_DATA_DOWNLOADER = 'ensembl_data_downloader'
    CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    CONFIG_PREFIX_RELEASE_FOLDER = 'folder_prefix_release'
    CONFIG_ENSEMBL_API = 'ensembl_api'
    CONFIG_ENSEMBL_API_SERVER = 'server'
    CONFIG_ENSEMBL_API_SPECIES = 'species'
    CONFIG_KEY_BASE_URL = 'ensembl_base_url'
    CONFIG_KEY_FOLDER_NAME_FASTA = 'folder_name_fasta'
    CONFIG_KEY_FOLDER_NAME_PROTEIN_SEQUENCES = 'folder_name_protein_sequences'
    CONFIG_KEY_FOLDER_NAME_GTF = 'folder_name_gtf'
    CONFIG_KEY_REWRITE_LOCAL_PATH_ENSEMBL_REPO = 'rewrite_local_path_ensembl_repo'
    CONFIG_KEY_ENSEMBL_FILE_NAMES = 'ensembl_file_names'
    CONFIG_KEY_PROTEIN_SEQUENCE_FILE = 'protein_sequence_file'
    CONFIG_KEY_FILE_TYPE = 'file_type'
    CONFIG_KEY_FILE_SUFFIXES = 'file_suffixes'
    CONFIG_KEY_FILE_EXTENSION = 'file_extension'
    CONFIG_KEY_GTF_FILE = 'gtf_file'
    CONFIG_REST_API_TAXON_ID = 'taxon_id'
    CONFIG_REST_API_NAME = 'name'
    CONFIG_TAXONOMY = 'taxonomy'
    CONFIG_ENSEMBL_NAME = 'name'
    CONFIG_KEY_SKIP_PROTEIN = 'skip_protein'
    CONFIG_KEY_SKIP_GTF = 'skip_gtf'
    CONFIG_KEY_SKIP_CDS = 'skip_cds'
    CONFIG_KEY_SKIP_CDNA = 'skip_cdna'
    CONFIG_KEY_SKIP_NCRNA = 'skip_ncrna'
    CONFIG_KEY_SKIP_DNA = 'skip_dna'
    CONFIG_KEY_SKIP_VCF = 'skip_vcf'
    CONFIG_KEY_GRCh37 = 'grch37'

    def __init__(self, config_file, pipeline_arguments):
        """
    Init the class with the specific parameters.
    :param config_file configuration file
    :param pipeline_arguments pipelines arguments
    """
        super(EnsemblDataDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_file,
                                                         pipeline_arguments)

        self._rest_api = 'http://rest.ensembl.org'
        self._rest_endpoint = '/info/species'

        self._skip_protein_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_PROTEIN,
                                                                        default_value=False)
        self._skip_gtf_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_GTF,
                                                                    default_value=False)
        self._skip_cds_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_CDS,
                                                                    default_value=False)
        self._skip_cdna_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_CDNA,
                                                                     default_value=False)
        self._skip_ncrna_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_NCRNA,
                                                                      default_value=False)
        self._skip_dna_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_DNA,
                                                                    default_value=False)
        self._skip_vcf_database = self.get_data_download_parameters(variable=self.CONFIG_KEY_SKIP_VCF,
                                                                    default_value=False)
        self._folder_prefix_release = self.get_data_download_parameters(variable=self.CONFIG_PREFIX_RELEASE_FOLDER,
                                                                    default_value='release-')
        self._local_path_ensembl = 'database_ensembl'
        if self.CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_ensembl = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_DIRECTORY]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters():
            self._local_path_ensembl = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_OUTPUT_DIRECTORY]

        self._ensembl_base_ftp = 'ftp://ftp.ensembl.org/pub'
        if self.CONFIG_KEY_BASE_URL in self.get_pipeline_parameters():
            self._ensembl_base_ftp = self.get_pipeline_parameters()[self.CONFIG_KEY_BASE_URL]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters():
            self._ensembl_base_ftp = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_KEY_BASE_URL]

        self._species_list = []
        if self.CONFIG_TAXONOMY in self.get_pipeline_parameters():
            species_parameters = self.get_pipeline_parameters()[self.CONFIG_TAXONOMY]
            self._species_list = species_parameters.split(",")

        self._ensembl_names = []
        if self.CONFIG_ENSEMBL_NAME in self.get_pipeline_parameters():
            species_name = self.get_pipeline_parameters()[self.CONFIG_ENSEMBL_NAME]
            self._ensembl_names = species_name.split(",")

        self._grch37 = self.get_data_download_parameters(variable=self.CONFIG_KEY_GRCh37, default_value=False)

        self.prepare_local_ensembl_repository()

    def get_data_download_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters():
            value_return = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][variable]
        return value_return

    def get_local_path_root_ensembl_repo(self):
        return self._local_path_ensembl

    def prepare_local_ensembl_repository(self):
        self.get_logger().debug(
            "Preparing local Ensembl repository, root folder - '{}'".format(self.get_local_path_root_ensembl_repo()))
        check_create_folders([self.get_local_path_root_ensembl_repo()])
        self.get_logger().debug(
            "Local path for Ensembl Release - '{}'".format(self.get_local_path_root_ensembl_repo()))

    def get_species_from_rest(self):
        """
    Get the list of species from ENSEMBL rest API.
    :return:
    """
        server = self._rest_api
        if self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and \
                self.CONFIG_ENSEMBL_API in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]:
            server = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_ENSEMBL_API][
                self.CONFIG_ENSEMBL_API_SERVER]

        endpoint = self._rest_endpoint
        if self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and \
                self.CONFIG_ENSEMBL_API_SPECIES in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]:
            endpoint = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_ENSEMBL_API][
                self.CONFIG_ENSEMBL_API_SPECIES]

        species_info = loads(call_api(server + endpoint).text)
        self._ensembl_species = species_info['species']

        return self._ensembl_species

    def download_database_by_species(self, url_file_name: str):
        """
    This method takes a list of Taxonomies from the commandline parameters
    and download the Protein fasta files and the gtf files.
    :return:
    """
        grch37 = self._grch37
        self.get_species_from_rest()

        if self._grch37:
            self.get_logger().info("gtch37 is only available for human not other species")

        url_file = None
        if url_file_name is not None:
            url_file = open(url_file_name, 'w')

        total_files = []
        files = []
        if self._species_list is not None and len(self._species_list) > 0:
            for species_id in self._species_list:
                for species in self._ensembl_species:
                    if species_id == species[self.CONFIG_REST_API_TAXON_ID]:
                        self.get_logger().debug(
                            "Downloading the data for the specie -- " + species[self.CONFIG_REST_API_TAXON_ID])
                        if not self._skip_protein_database:
                            prot_files = self.get_pep_files(species=species, grch37=grch37, url_file=url_file)
                            files.extend(prot_files)
                        if not self._skip_gtf_database:
                            gtf_files = self.get_gtf_files(species, grch37)
                            files.extend(gtf_files)
                        if not self._skip_cds_database:
                            cds_files = self.get_cds_files(species=species, grch37=grch37)
                            files.extend(cds_files)
                        if not self._skip_cdna_database:
                            cdna_files = self.get_cdna_files(species=species, grch37=grch37)
                            files.extend(cdna_files)
                        if not self._skip_ncrna_database:
                            ncrna_files = self.get_ncrna_files(species=species, grch37=grch37)
                            files.extend(ncrna_files)
                        if not self._skip_dna_database:
                            dna_files = self.get_genome_assembly_files(species=species, grch37=grch37)
                            files.extend(dna_files)
                        if not self._skip_vcf_database:
                            vcf_files = self.get_vcf_files(species=species)
                            files.extend(vcf_files)

                        self.get_logger().debug("Files downloaded -- " + ",".join(
                            [x for x in files if x]))
                        total_files.extend(files)
        elif self._ensembl_names is not None and len(self._ensembl_names) > 0:
            for species_name in self._ensembl_names:
                for species in self._ensembl_species:
                    if species_name == species[self.CONFIG_REST_API_NAME]:
                        self.get_logger().debug(
                            "Downloading the data for the specie --" + species[self.CONFIG_REST_API_NAME])
                        if not self._skip_protein_database:
                            prot_files = self.get_pep_files(species=species, url_file=url_file)
                            files.extend(prot_files)
                        if not self._skip_gtf_database:
                            gtf_files = self.get_gtf_files(species=species, url_file=url_file)
                            files.extend(gtf_files)
                        if not self._skip_cds_database:
                            cds_files = self.get_cds_files(species=species, url_file=url_file)
                            files.extend(cds_files)
                        if not self._skip_cdna_database:
                            cdna_files = self.get_cdna_files(species=species, url_file=url_file)
                            files.extend(cdna_files)
                        if not self._skip_ncrna_database:
                            ncrna_files = self.get_ncrna_files(species=species, url_file=url_file)
                            files.extend(ncrna_files)
                        if not self._skip_dna_database:
                            dna_files = self.get_genome_assembly_files(species=species, url_file=url_file)
                            files.extend(dna_files)
                        if not self._skip_vcf_database:
                            vcf_files = self.get_vcf_files(species=species, url_file=url_file)
                            files.extend(vcf_files)

                        self.get_logger().debug("Files downloaded -- " + ",".join(
                            [x for x in files if x]))
                        total_files.extend(files)
        else:
            for species in self._ensembl_species:
                self.get_logger().debug(
                    "Downloading the data for the specie --" + species[self.CONFIG_REST_API_TAXON_ID])
                if not self._skip_protein_database:
                    prot_files = self.get_pep_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(prot_files)
                if not self._skip_gtf_database:
                    gtf_files = self.get_gtf_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(gtf_files)
                if not self._skip_cds_database:
                    cds_files = self.get_cds_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(cds_files)
                if not self._skip_cdna_database:
                    cdna_files = self.get_cdna_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(cdna_files)
                if not self._skip_ncrna_database:
                    ncrna_files = self.get_ncrna_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(ncrna_files)
                if not self._skip_dna_database:
                    dna_files = self.get_genome_assembly_files(species=species, grch37=grch37, url_file=url_file)
                    files.extend(dna_files)
                if not self._skip_vcf_database:
                    vcf_files = self.get_vcf_files(species=species, url_file=url_file)
                    files.extend(vcf_files)

                self.get_logger().debug("Files downloaded -- " + ",".join(
                    [x for x in files if x]))
                total_files.extend(files)

        return total_files

    def get_release_path(self, species):
        """
        Get the release path in ensembl
        :param species: specie with to be downloaded
        :return: release path in ensembl
        """
        if self._folder_prefix_release == 'release-':
            return 'release-{}'.format(species['release'])
        return self._folder_prefix_release

    def get_cds_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    Get the cds files for a specific species object.
    :return: List of files names.
    """
        files = []
        release_version = self.get_release_path(species = species)
        try:
            if grch37:
                species['assembly'] = 'GRCh37'
            file_name = '{}.{}.cds.all.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                     species['assembly'])
            file_url = '{}/{}/fasta/{}/cds/{}'.format(self._ensembl_base_ftp, release_version,
                                                              species['name'], file_name)
            if grch37:
                file_url = '{}/grch37/{}/fasta/{}/cds/{}'.format(self._ensembl_base_ftp, release_version,
                                                                         species['name'], file_name)
            files.append(
                download_file(file_url=file_url, file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                              log=self.get_logger(), url_file=url_file))
        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def get_cdna_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    Get the cds files for a specific species object.
    :return: List of files names.
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            if grch37:
                species['assembly'] = 'GRCh37'
            file_name = '{}.{}.cdna.all.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                      species['assembly'])
            file_url = '{}/{}/fasta/{}/cdna/{}'.format(self._ensembl_base_ftp, release_version, species['name'], file_name)
            if grch37:
                file_url = '{}/grch37/{}/fasta/{}/cdna/{}'.format(self._ensembl_base_ftp,
                                                                          release_version, species['name'],
                                                                          file_name)
            files.append(
                download_file(file_url=file_url, file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                              log=self.get_logger(), url_file=url_file))
        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def get_ncrna_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    Get the cds files for a specific species object.
    :return: List of files names.
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            if grch37:
                species['assembly'] = 'GRCh37'

            file_name = '{}.{}.ncrna.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                   species['assembly'])
            file_url = '{}/{}/fasta/{}/ncrna/{}'.format(self._ensembl_base_ftp, release_version, species['name'], file_name)
            if grch37:
                file_url = '{}/grch37/{}/fasta/{}/ncrna/{}'.format(self._ensembl_base_ftp, release_version, species['name'],
                                                                           file_name)
            files.append(
                download_file(file_url=file_url, file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                              log=self.get_logger(), url_file=url_file))
        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def get_pep_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    Get the peptide files for a specific species object.
    :return: List of files names.
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            # TODO: Would be better to check by API the assembly version
            if grch37:
                species['assembly'] = 'GRCh37'
            file_name = '{}.{}.pep.all.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                     species['assembly'])
            file_url = '{}/{}/fasta/{}/pep/{}'.format(self._ensembl_base_ftp, release_version, species['name'], file_name)
            if grch37:
                file_url = '{}/grch37/{}/fasta/{}/pep/{}'.format(self._ensembl_base_ftp, release_version, species['name'], file_name)
            print('Downloading -- ', file_url)
            files.append(
                download_file(file_url=file_url, file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                              log=self.get_logger(), url_file=url_file))

        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def get_gtf_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    This method retrieve the gtf files for a specific specie object
    :param grch37: if the GrCh37 genome assembly is desired enable to true
    :param species: species to download the file.
    :param url_file: save the urls to a file for download
    :return:
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            if grch37:
                species['assembly'] = 'GRCh37'
                species['release'] = '87'
            file_name = '{}.{}.{}.gtf.gz'.format(species['name'][0].upper() + species['name'][1:], species['assembly'],
                                                 release_version, )
            file_url = '{}/{}/gtf/{}/{}'.format(self._ensembl_base_ftp, release_version, species['name'],
                                                        file_name)
            if grch37:
                file_url = '{}/grch37/{}/gtf/{}/{}'.format(self._ensembl_base_ftp, release_version,
                                                                   species['name'], file_name)
            if url_file is None:
                files.append(download_file(file_url=file_url,
                                           file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                           log=self.get_logger(), url_file=url_file))
            else:
                url_file.write("{}\t{}\n".format(file_url, self.get_local_path_root_ensembl_repo() + '/' + file_name))
                files.append(self.get_local_path_root_ensembl_repo() + '/' + file_name)

        except KeyError:
            self.get_logger().debug("No valid info is available species: ", species)

        return files

    def get_vcf_files(self, species: dict, url_file=None) -> list:
        """
    This method retrieve the vcf file for a specific specie object
    :param species:
    :param url_file: save the urls to a file for download
    :return:
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            file_name = '{}_incl_consequences.vcf.gz'.format(species['name'])
            file_url = '{}/{}/variation/vcf/{}/'.format(self._ensembl_base_ftp, release_version, species['name'])

            downloaded_file = download_file(file_url=file_url + file_name,
                                            file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                            log=self.get_logger(), url_file=url_file)
            if downloaded_file is not None:
                files.append(downloaded_file)

            elif species['name'] == 'homo_sapiens':
                # for humans the variants are stored per chromosome, so we need to download them all and combine them into one file here
                chrN = 1
                file_name = '{}_incl_consequences-chr{}.vcf.gz'.format(species['name'], chrN)
                downloaded_file = download_file(file_url=file_url + file_name,
                                                file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                                log=self.get_logger(), url_file=url_file)
                if downloaded_file is not None:
                    # if chr1 is downloaded then try all others
                    files.append(downloaded_file)
                    for chrN in range(2, 23):  # chr2-22
                        file_name = '{}_incl_consequences-chr{}.vcf.gz'.format(species['name'], chrN)
                        files.append(
                            download_file(file_url=file_url + file_name,
                                          file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                          log=self.get_logger(), url_file=url_file))
                    file_name = '{}_incl_consequences-chr{}.vcf.gz'.format(species['name'], 'X')
                    files.append(download_file(file_url=file_url + file_name,
                                               file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                               log=self.get_logger(), url_file=url_file))
                    file_name = '{}_incl_consequences-chr{}.vcf.gz'.format(species['name'], 'Y')
                    files.append(download_file(file_url=file_url + file_name,
                                               file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                               log=self.get_logger(), url_file=url_file))
                    file_name = '{}_incl_consequences-chr{}.vcf.gz'.format(species['name'], 'MT')
                    files.append(download_file(file_url=file_url + file_name,
                                               file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                                               log=self.get_logger(), url_file=url_file))

        except KeyError:
            self.get_logger().debug("No valid info is available species: ", species)

        return files

    def get_genome_assembly_files(self, species: dict, grch37=False, url_file=None) -> list:
        """
    This method retrieve the genome assembly files for a specific specie object
    :param grch37: if the GrCh37 genome assembly is desired enable to true
    :param species: species to download the file.
    :param url_file: save the urls to a file for download
    :return:
    """
        files = []
        release_version = self.get_release_path(species=species)
        try:
            if grch37:
                species['assembly'] = 'GRCh37'
            file_name = '{}.{}.dna_sm.toplevel.fa.gz'.format(species['name'][0].upper() + species['name'][1:],
                                                             species['assembly'])
            file_url = '{}/{}/fasta/{}/dna/{}'.format(self._ensembl_base_ftp,
                                                              release_version, species['name'], file_name)
            if grch37:
                file_url = '{}/grch37/release-{}/fasta/{}/dna/{}'.format(self._ensembl_base_ftp,
                                                                         release_version, species['name'], file_name)
            files.append(
                download_file(file_url=file_url, file_name=self.get_local_path_root_ensembl_repo() + '/' + file_name,
                              log=self.get_logger(), url_file=url_file))
        except KeyError:
            print("No valid info is available species: ", species)

        return files

    def validate_taxonomies(self):
        if self._species_list is not None:
            list_of_taxonomies = self.get_species_from_rest()
            taxonomies = list(map(lambda name: name['taxon_id'], list_of_taxonomies))
            for taxonomy in self._species_list:
                if taxonomy not in taxonomies:
                    print('Taxonomy -- {} not found in the ensembl taxonomies'.format(taxonomy))
                    print('Valid Taxonomies:')
                    for i in list_of_taxonomies:
                        print('Taxonomy {} --- {}'.format(i['taxon_id'], i['common_name']))
                    return False
        return True


