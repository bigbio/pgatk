from pyopenms import *

from pypgatk.toolbox.general import ParameterConfiguration


class OpenmsDataService(ParameterConfiguration):
  CONFIG_KEY_OPENMS_ANALYSIS = 'openms_analysis'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_PEPTIDE_CLASS_PREFIX = "peptide_classes_prefix"
  CONFIG_PEPTIDE_FDR_CUTOFF = "psm_pep_fdr_cutoff"
  CONFIG_PEPTIDE_CLASS_FDR_CUTOFF = "psm_pep_class_fdr_cutoff"

  def __init__(self, config_file, pipeline_arguments):
    super(OpenmsDataService, self).__init__(self.CONFIG_KEY_OPENMS_ANALYSIS, config_file,
                                                pipeline_arguments)

    self._decoy_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_DECOY_PREFIX]
    if self.CONFIG_DECOY_PREFIX in self.get_pipeline_parameters():
      self._decoy_prefix = self.get_pipeline_parameters()[self.CONFIG_DECOY_PREFIX]

    self._peptide_class_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_CLASS_PREFIX]
    if self.CONFIG_PEPTIDE_CLASS_PREFIX in self.get_pipeline_parameters():
      self._peptide_class_prefix = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_CLASS_PREFIX]

    self._min_peptide_length = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_MIN_PEPTIDE_LENGTH]
    if self.CONFIG_MIN_PEPTIDE_LENGTH in self.get_pipeline_parameters():
      self._min_peptide_length = self.get_pipeline_parameters()[self.CONFIG_MIN_PEPTIDE_LENGTH]

    self._psm_pep_fdr_cutoff = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_FDR_CUTOFF]
    if self.CONFIG_PEPTIDE_FDR_CUTOFF in self.get_pipeline_parameters():
      self._psm_pep_fdr_cutoff = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_FDR_CUTOFF]

    self._psm_pep_class_fdr_cutoff = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_CLASS_FDR_CUTOFF]
    if self.CONFIG_PEPTIDE_CLASS_FDR_CUTOFF in self.get_pipeline_parameters():
      self._psm_pep_fdr_cutoff = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_CLASS_FDR_CUTOFF]

  def compute_global_fdr(self, peptide_ids):
    """
    Compute the global FDR and filter peptides
    :param peptide_ids: list of peptide identifications
    :return:  filtered peptides
    """

    target_count = 0
    decoy_count = 0

    # Map list of peptides with hits to -> map(hit, peptide)
    peptide_hit_dict = []

    # Here the Map contains Key=spectrum_reference (unique for each PeptideIdentification)
    # Value is the PeptideIdentification
    peptides_filtered = {}

    for peptide in peptide_ids:
      for peptide_hit in peptide.getHits():
        peptide_hit_dict.append((peptide_hit, peptide))

    # Todo: Here would be interesting to know the score order isHigherScoreBetter
    peptide_hit_dict.sort(key=lambda x: x[0].getScore())

    for peptide_tuple in peptide_hit_dict:
      (psm, peptide_identificaion) = peptide_tuple
      accessions = [ev.getProteinAccession() for ev in psm.getPeptideEvidences()]

      if any(self._decoy_prefix in s for s in accessions):
        decoy_count += 1
      else:
        target_count += 1

      FDR = float(decoy_count) / float(target_count)
      if FDR < self._psm_pep_fdr_cutoff:
        identifier = peptide_identificaion.getMetaValue("spectrum_reference")
        print(identifier)
        if identifier not in peptides_filtered:
          peptide_id = peptide_identificaion
          peptide_id.setHits([psm])
          peptides_filtered[identifier] = peptide_id
        else:
          peptide_id = peptides_filtered[identifier]
          hits = peptide_id.getHits()
          hits.append(psm)
          peptide_id.setHits(hits)
          peptides_filtered[identifier] = peptide_id

      print(FDR)

    return list(peptides_filtered.values())

  def compute_class_fdr(self, peptide_ids):

    # Here the Map contains Key=spectrum_reference (unique for each PeptideIdentification)
    # Value is the PeptideIdentification

    peptides_filtered = {}
    peptide_hit_dict = []

    for peptide in peptide_ids:
      for peptide_hit in peptide.getHits():
        peptide_hit_dict.append((peptide_hit, peptide))

    # Todo: Here would be interesting to know the score order isHigherScoreBetter
    peptide_hit_dict.sort(key=lambda x: x[0].getScore())

    for peptide_hit in peptide_hit_dict:
      accessions = [ev.getProteinAccession() for ev in peptide_hit[0].getPeptideEvidences()]

  def filter_peptide_class_fdr(self, input_idxml, output_idxml):
    """
    Filter peptides by Class FDR in idXML files
    :param input_idxml:  input idxml
    :param output_idxml: output idxml
    :return:
    """

    protein_ids = []
    peptide_ids = []

    idfilter = IDFilter()
    IdXMLFile().load(input_idxml, protein_ids, peptide_ids)
    print(peptide_ids[0])

    # Iterate over PeptideIdentification
    filtered_peptide_ids = list(
      filter(lambda peptide: (len(peptide.getHits()[0].getSequence().toUnmodifiedString()) > self._min_peptide_length),
             peptide_ids))

    filtered_peptide_ids = self.compute_global_fdr(filtered_peptide_ids)

    remove_peptides_without_reference = True
    idfilter.updateProteinReferences(filtered_peptide_ids, protein_ids, remove_peptides_without_reference)
    IdXMLFile().store(output_idxml, protein_ids, filtered_peptide_ids)






