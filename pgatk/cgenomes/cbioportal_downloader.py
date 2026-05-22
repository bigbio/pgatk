import csv
import os
from concurrent.futures import as_completed
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Optional

import requests

from pgatk.toolbox.exceptions import AppException
from pgatk.toolbox.general import ParameterConfiguration, check_create_folders, clear_cache
from pgatk.toolbox.rest import call_api_raw

_CBIO_PAGE_SIZE = 10_000

MAF_HEADER = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Consequence",
    "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
    "Verification_Status", "Validation_Status", "Mutation_Status",
    "Sequencing_Phase", "Sequence_Source", "Validation_Method",
    "Score", "BAM_File", "Sequencer", "t_ref_count", "t_alt_count",
    "n_ref_count", "n_alt_count", "HGVSc", "HGVSp", "HGVSp_Short",
    "Transcript_ID", "RefSeq", "Protein_position", "Codons", "Hotspot",
]


def _json_or_raise(resp, label: str, log):
    """Parse JSON from a response, raising with a clear message if the body is empty or invalid."""
    if not resp.text:
        raise ValueError(
            f"{label}: server returned HTTP {resp.status_code} with an empty body"
        )
    try:
        return resp.json()
    except Exception:
        snippet = resp.text[:300]
        raise ValueError(f"{label}: could not parse JSON (HTTP {resp.status_code}). Body: {snippet!r}")


def _get_sample_ids(base_url: str, study_id: str, log) -> list:
    """Return all sequenced sample IDs for a study via the sample-list API."""
    # Prefer the sequenced list; fall back to the _all list.
    for list_id in (f"{study_id}_sequenced", f"{study_id}_all"):
        url = f"{base_url}/sample-lists/{list_id}/sample-ids"
        log.debug("GET %s", url)
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        if resp.status_code == 200 and resp.text:
            return _json_or_raise(resp, f"sample-ids ({list_id})", log)
    log.warning("Could not retrieve sample list for study '%s'", study_id)
    return []


def _fetch_study_mutations(base_url: str, study_id: str, log) -> list:
    """Fetch all mutations for a study from the cBioPortal REST API."""
    profiles_url = f"{base_url}/studies/{study_id}/molecular-profiles"
    log.debug("GET %s", profiles_url)
    resp = requests.get(profiles_url, headers={"Accept": "application/json"}, timeout=30)
    resp.raise_for_status()
    profiles = _json_or_raise(resp, "molecular-profiles", log)

    profile_id = None
    for p in profiles:
        if p.get("molecularAlterationType") == "MUTATION_EXTENDED":
            profile_id = p["molecularProfileId"]
            break
    if profile_id is None:
        log.warning("No MUTATION_EXTENDED profile found for study '%s'", study_id)
        return []
    log.info("Using molecular profile '%s'", profile_id)

    sample_ids = _get_sample_ids(base_url, study_id, log)
    if not sample_ids:
        log.warning("No samples found for study '%s'", study_id)
        return []
    log.info("Fetching mutations for %d samples ...", len(sample_ids))

    # sampleListId is silently broken in the current API; use explicit identifiers.
    sample_mol_ids = [
        {"molecularProfileId": profile_id, "sampleId": sid} for sid in sample_ids
    ]

    mutations = []
    page = 0
    while True:
        url = (
            f"{base_url}/mutations/fetch"
            f"?molecularProfileId={profile_id}&projection=DETAILED"
            f"&pageSize={_CBIO_PAGE_SIZE}&pageNumber={page}"
        )
        log.debug("POST %s (page %d)", url, page)
        resp = requests.post(
            url,
            json={"sampleMolecularIdentifiers": sample_mol_ids},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=120,
        )
        resp.raise_for_status()
        batch = _json_or_raise(resp, f"mutations page {page}", log)
        if not batch:
            break
        mutations.extend(batch)
        log.info("  Page %d: %d mutations (total: %d)", page, len(batch), len(mutations))
        if len(batch) < _CBIO_PAGE_SIZE:
            break
        page += 1
    return mutations


def _mutation_to_row(m: dict) -> list:
    """Map a cBioPortal API mutation object to a MAF-compatible row list."""
    gene = m.get("gene") or {}
    hugo = gene.get("hugoGeneSymbol", "")
    entrez = str(gene.get("entrezGeneId", ""))
    chrom = str(m.get("chr", ""))
    start = str(m.get("startPosition", ""))
    end = str(m.get("endPosition", ""))
    ref = m.get("referenceAllele", "")
    alt = m.get("variantAllele", "")
    varclass = m.get("mutationType", "")
    vartype = m.get("variantType", "")
    sample_id = m.get("sampleId", "")
    ncbi_build = m.get("ncbiBuild", "GRCh37")
    center = m.get("center", "")
    mut_status = m.get("mutationStatus", "Somatic")
    val_status = m.get("validationStatus", "")
    t_ref = str(m.get("tumorRefCount", ""))
    t_alt = str(m.get("tumorAltCount", ""))
    n_ref = str(m.get("normalRefCount", ""))
    n_alt = str(m.get("normalAltCount", ""))
    refseq = m.get("refseqMrnaId", "")
    hgvsc = m.get("hgvsc", "") or ""
    pc = m.get("proteinChange", "")
    hgvsp_short = pc if pc.startswith("p.") else (f"p.{pc}" if pc else "")
    # Use RefSeq transcript ID; Ensembl IDs are not returned by the public API.
    transcript_id = refseq.split(".")[0] if refseq else ""
    protein_pos = str(m.get("proteinPosStart", ""))

    return [
        hugo, entrez, center, ncbi_build, chrom, start, end, "+",
        "", varclass, vartype, ref, ref, alt, "", "", sample_id,
        "", "", "", "", "", "", "", "", val_status, mut_status, "", "", "", "", "", "",
        t_ref, t_alt, n_ref, n_alt, hgvsc, "", hgvsp_short, transcript_id, refseq,
        protein_pos, "", "",
    ]


def _fetch_clinical_data(base_url: str, study_id: str, output_dir: str, log) -> Optional[str]:
    """Download all sample-level clinical attributes for a study into a TSV file.

    Calls ``GET /studies/{studyId}/clinical-data?clinicalDataType=SAMPLE`` with
    pagination, then pivots the attribute-per-row API response into a standard
    SAMPLE_ID × attribute-column TSV that ``CancerGenomesService.get_value_per_sample``
    can read.
    """
    records: list[dict] = []
    page = 0
    while True:
        url = (
            f"{base_url}/studies/{study_id}/clinical-data"
            f"?clinicalDataType=SAMPLE&pageSize={_CBIO_PAGE_SIZE}&pageNumber={page}"
        )
        log.debug("GET %s", url)
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
        if resp.status_code != 200:
            log.warning("Clinical data request failed (HTTP %d) for study '%s'", resp.status_code, study_id)
            return None
        if not resp.text:
            break
        batch = _json_or_raise(resp, f"clinical-data page {page}", log)
        if not batch:
            break
        records.extend(batch)
        log.debug("  Clinical page %d: %d records (total %d)", page, len(batch), len(records))
        if len(batch) < _CBIO_PAGE_SIZE:
            break
        page += 1

    if not records:
        log.warning("No clinical data returned for study '%s'", study_id)
        return None

    # Pivot: collect all attributes in insertion order, then build one row per sample.
    sample_attrs: dict[str, dict[str, str]] = {}
    attr_order: list[str] = []
    for rec in records:
        sid = rec.get("sampleId", "")
        attr = rec.get("clinicalAttributeId", "")
        val = rec.get("value", "")
        if attr and attr not in attr_order:
            attr_order.append(attr)
        if sid not in sample_attrs:
            sample_attrs[sid] = {}
        sample_attrs[sid][attr] = val

    out_path = os.path.join(output_dir, "data_clinical_sample.txt")
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        header = ["SAMPLE_ID"] + attr_order
        fh.write("\t".join(header) + "\n")
        for sid, attrs in sample_attrs.items():
            row = [sid] + [attrs.get(a, "") for a in attr_order]
            fh.write("\t".join(row) + "\n")

    log.info("Wrote clinical data for %d samples to %s", len(sample_attrs), out_path)
    return out_path


class CbioPortalDownloadService(ParameterConfiguration):
    CONFIG_KEY_DATA_DOWNLOADER = 'cbioportal_data_downloader'
    CONFIG_KEY_CBIOPORTAL_DOWNLOAD_URL = 'cbioportal_download_url'
    CONFIG_OUTPUT_DIRECTORY = 'output_directory'
    CONFIG_CBIOPORTAL_API = 'cbioportal_api'
    CONFIG_CBIOPORTAL_API_SERVER = 'base_url'
    CONFIG_CBIOPORTAL_API_CANCER_STUDIES = "cancer_studies"
    CONFIG_LIST_STUDIES = "list_studies"
    CONFIG_MULTITHREADING = "multithreading"
    PROTEINDB = 'proteindb'
    FILTER_INFO = 'filter_info'
    FILTER_COLUMN = 'filter_column'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(CbioPortalDownloadService, self).__init__(self.CONFIG_KEY_DATA_DOWNLOADER, config_data,
                                                        pipeline_arguments)

        self._local_path_cbioportal = 'output_directory'
        self._list_studies = []
        self._multithreading = True

        self._cbioportal_base_url = 'https://www.cbioportal.org/api'
        self._cancer_studies_command = 'studies'

        if self.CONFIG_OUTPUT_DIRECTORY in self.get_pipeline_parameters():
            self._local_path_cbioportal = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_DIRECTORY]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and \
                self.CONFIG_OUTPUT_DIRECTORY in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]:
            self._local_path_cbioportal = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_OUTPUT_DIRECTORY]

        if self.CONFIG_LIST_STUDIES in self.get_pipeline_parameters():
            self._list_studies = self.get_pipeline_parameters()[self.CONFIG_LIST_STUDIES]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and \
                self.CONFIG_LIST_STUDIES in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]:
            self._list_studies = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_LIST_STUDIES]

        if self.CONFIG_MULTITHREADING in self.get_pipeline_parameters():
            self._multithreading = self.get_pipeline_parameters()[self.CONFIG_MULTITHREADING]
        elif self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and \
                self.CONFIG_MULTITHREADING in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]:
            self._multithreading = self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                self.CONFIG_MULTITHREADING]

        if self.CONFIG_CBIOPORTAL_API_SERVER in self.get_pipeline_parameters():
            self._cbioportal_base_url = self.get_pipeline_parameters()[self.CONFIG_CBIOPORTAL_API_SERVER]
        elif (self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and
              self.CONFIG_CBIOPORTAL_API in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]
              and self.CONFIG_CBIOPORTAL_API_SERVER in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                  self.CONFIG_CBIOPORTAL_API]):
            self._cbioportal_base_url = \
                self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
                    self.CONFIG_CBIOPORTAL_API_SERVER]

        if self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES in self.get_pipeline_parameters():
            self._cancer_studies_command = self.get_pipeline_parameters()[self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES]
        elif (self.CONFIG_KEY_DATA_DOWNLOADER in self.get_default_parameters() and
              self.CONFIG_CBIOPORTAL_API in self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER]
              and self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES in
              self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][
                  self.CONFIG_CBIOPORTAL_API]):
            self._cancer_studies_command = \
                self.get_default_parameters()[self.CONFIG_KEY_DATA_DOWNLOADER][self.CONFIG_CBIOPORTAL_API][
                    self.CONFIG_CBIOPORTAL_API_CANCER_STUDIES]

        self.prepare_local_cbioportal_repository()
        self.get_cancer_studies()

    def prepare_local_cbioportal_repository(self):
        self.get_logger().debug("Preparing local cbioportal repository, root folder - '{}'".format(
            self.get_local_path_root_cbioportal_repo()))
        check_create_folders([self.get_local_path_root_cbioportal_repo()])
        self.get_logger().debug(
            "Local path for cbioportal Release - '{}'".format(self.get_local_path_root_cbioportal_repo()))

    def get_local_path_root_cbioportal_repo(self):
        return self._local_path_cbioportal

    def get_filter_options(self, variable, default_value):
        return_value = default_value
        if variable in self.get_default_parameters():
            return_value = self.get_default_parameters()[variable]
        elif self.PROTEINDB in self.get_default_parameters() and \
                self.FILTER_INFO in self.get_default_parameters()[self.PROTEINDB] and \
                variable in self.get_default_parameters()[self.PROTEINDB][self.FILTER_INFO]:
            return_value = self.get_default_parameters()[self.PROTEINDB][self.FILTER_INFO][variable]
        return return_value

    def get_cancer_studies(self):
        """
        This method will print the list of all cancer studies for the user.
        :return:
        """
        server = self._cbioportal_base_url
        endpoint = self._cancer_studies_command
        self._cbioportal_studies = call_api_raw(server + "/" + endpoint).text
        return self._cbioportal_studies

    def download_study(self, download_study, url_file_name=None):
        """
        This function will download a study from cBioPortal using the study ID
        :param download_study: Study to be downloaded, if the study is empty or None, all the studies will be
        downloaded.
        :param url_file_name: file tsv containing the urls to be downloaded.
        :return: None
        """

        clear_cache()

        if self._cbioportal_studies is None or len(self._cbioportal_studies) == 0:
            self.get_cancer_studies()

        if url_file_name is not None:
            with open(url_file_name, 'w', encoding='utf-8') as url_file:
                if 'all' not in download_study:
                    if not self.check_study_identifier(download_study):
                        msg = "The following study accession '{}' is not present in cBioPortal Studies".format(download_study)
                        self.get_logger().debug(msg)
                        raise AppException(msg)
                    else:
                        self.download_one_study(download_study, url_file=url_file)
                else:
                    csv_reader = csv.reader(self._cbioportal_studies.splitlines(), delimiter="\t")
                    line_count = 0
                    if self._multithreading:
                        processes = []
                        with ThreadPoolExecutor(max_workers=10, thread_name_prefix='Thread-Download') as executor:
                            for row in csv_reader:
                                if line_count != 0:
                                    processes.append(executor.submit(self.download_one_study, row[0], url_file=url_file))
                                line_count = line_count + 1
                        for task in as_completed(processes):
                            print(task.result())
                    else:
                        for row in csv_reader:
                            if line_count != 0:
                                self.download_one_study(row[0], url_file=url_file)
                            line_count = line_count + 1
        else:
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

    def download_one_study(self, download_study: str, url_file=None) -> Optional[str]:
        log = self.get_logger()
        study_dir = os.path.join(self.get_local_path_root_cbioportal_repo(), download_study)
        check_create_folders([study_dir])
        out_path = os.path.join(study_dir, "data_mutations.txt")

        log.info("Fetching mutations for study '%s' via cBioPortal API ...", download_study)
        try:
            mutations = _fetch_study_mutations(self._cbioportal_base_url, download_study, log)
        except Exception as exc:
            log.error("Failed to fetch mutations for study '%s': %s", download_study, exc)
            return None

        if not mutations:
            log.warning("No mutations returned for study '%s'", download_study)
            return None

        with open(out_path, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(MAF_HEADER) + "\n")
            for m in mutations:
                row = _mutation_to_row(m)
                fh.write("\t".join("" if x is None else str(x) for x in row) + "\n")

        log.info("Wrote %d mutations to %s", len(mutations), out_path)

        _fetch_clinical_data(self._cbioportal_base_url, download_study, study_dir, log)

        if url_file is not None:
            url_file.write(out_path + "\n")
        return out_path

    def check_study_identifier(self, download_study):
        return download_study in self._cbioportal_studies
