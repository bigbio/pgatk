from pyopenms import IdXMLFile as idxml_parser
from pyopenms import IDFilter

from pypgatk.toolbox.general import ParameterConfiguration
import numpy.polynomial.polynomial as poly
import numpy as np


class OpenmsDataService(ParameterConfiguration):
  CONFIG_KEY_OPENMS_ANALYSIS = 'openms_analysis'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_PEPTIDE_CLASS_PREFIX = "peptide_classes_prefix"
  CONFIG_PEPTIDE_FDR_CUTOFF = "psm_pep_fdr_cutoff"
  CONFIG_PEPTIDE_CLASS_FDR_CUTOFF = "psm_pep_class_fdr_cutoff"
  CONFIG_PEPTIDE_GROUP_PREFIX = "peptide_groups_prefix"
  CONFIG_PEPTIDE_APPLY_CLASS_FDR = "enable_class_fdr"

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

    self._peptide_groups_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_GROUP_PREFIX]
    if self.CONFIG_PEPTIDE_GROUP_PREFIX in self.get_pipeline_parameters():
      self._peptide_groups_prefix = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_GROUP_PREFIX]

    self._peptide_class_fdr_enable = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_APPLY_CLASS_FDR]
    if self.CONFIG_PEPTIDE_APPLY_CLASS_FDR in self.get_pipeline_parameters():
      self._peptide_class_fdr_enable = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_APPLY_CLASS_FDR]

  @staticmethod
  def is_peptide_group(peptide_group_members, accessions):
    """
    Given a group of classes and a list of accessions of a peptide. Returns True if all accessions match to exactly one class in the group.
    :param peptide_group_members: all protein classes
    :param accessions:  all protein accessions associated with the peptide.
    :return: True if all protein accessions belows to one of these peptide_group_members.
    """

    accession_group = 0
    for accession in accessions:
      for class_peptide in peptide_group_members:
        if class_peptide in accession:
          accession_group += 1
    return len(accessions) == accession_group

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
      (psm, peptide_identification) = peptide_tuple
      accessions = [ev.getProteinAccession() for ev in psm.getPeptideEvidences()]

      if any(self._decoy_prefix in s for s in accessions):
        decoy_count += 1
      else:
        target_count += 1

      FDR = float(decoy_count) / float(target_count)
      if FDR < self._psm_pep_fdr_cutoff:
        identifier = peptide_identification.getMetaValue("spectrum_reference")
        if identifier not in peptides_filtered:
          peptide_id = peptide_identification
          peptide_id.setHits([psm])
          peptides_filtered[identifier] = peptide_id
        else:
          peptide_id = peptides_filtered[identifier]
          hits = peptide_id.getHits()
          hits.append(psm)
          peptide_id.setHits(hits)
          peptides_filtered[identifier] = peptide_id

    return list(peptides_filtered.values())

  def compute_class_fdr(self, peptide_ids):
    """
    Compute the global FDR and filter peptides
    :param peptide_ids: list of peptide identifications
    :return:  filtered peptides
    """

    # Map list of peptides with hits to -> map(hit, peptide)
    peptide_hit_dict = []

    # Here the Map contains Key=spectrum_reference (unique for each PeptideIdentification)
    # Value is the PeptideIdentification
    for peptide in peptide_ids:
      for peptide_hit in peptide.getHits():
        peptide_hit_dict.append((peptide_hit, peptide))

    # Todo: Here would be interesting to know the score order isHigherScoreBetter
    peptide_hit_dict.sort(key=lambda x: x[0].getScore())

    score_dic = {}
    decoy_dic = {}

    # A dictionary of peptide classes with number of decoys and targets per class
    novel_target_decoy_count = {}
    for peptide_class in self._peptide_groups_prefix:
      novel_target_decoy_count[peptide_class] = (0,0)
      decoy_dic[peptide_class] = {}
      score_dic[peptide_class] = {}

    global_decoy_count = 0
    global_target_count = 0

    for peptide_hit in peptide_hit_dict:
      (psm, peptide_identification) = peptide_hit
      accessions = [ev.getProteinAccession() for ev in psm.getPeptideEvidences()]
      score = -np.log10(float(psm.getScore()))
      if any(self._decoy_prefix in s for s in accessions):
        global_decoy_count += 1
        for peptide_group in self._peptide_groups_prefix:
          if self.is_peptide_group(self._peptide_groups_prefix[peptide_group], accessions):
            (target, decoy) = novel_target_decoy_count[peptide_group]
            decoy += 1
            novel_target_decoy_count[peptide_group] = (target, decoy)
            decoy_dic[peptide_group][score] = [decoy, global_decoy_count]
      else:
        global_target_count += 1
        for peptide_group in self._peptide_groups_prefix:
          if self.is_peptide_group(self._peptide_groups_prefix[peptide_group], accessions):
            (target, decoy) = novel_target_decoy_count[peptide_group]
            target += 1
            novel_target_decoy_count[peptide_group] = (target, decoy)

      for peptide_group in self._peptide_groups_prefix:
        (target, decoy) = novel_target_decoy_count[peptide_group]
        score_dic[peptide_group][score] = [global_target_count, global_decoy_count, target, decoy]

    peptide_dict_models = {}
    for peptide_class in score_dic:
      x = []
      y = []
      x_filter = []
      y_filter = []
      for score in decoy_dic[peptide_class]:
        x.append(score)
        frac = float(decoy_dic[peptide_class][score][0] / decoy_dic[peptide_class][score][1])
        y.append(frac)
        # Todo: In the original algorithm the authors used a filter for the score between 6-10
        #  https://github.com/yafeng/proteogenomics_python/blob/4b1638aa75903225e9ae45892af4cb9f078d7421/BayesClassSpecificFDR.py#L88
        # if 6 < score < 10:
        #   x_filter.append(score)
        #   y_filter.append(frac)
        x_filter.append(score)
        y_filter.append(frac)
      if len(x_filter) > 0 and len(y_filter) > 0:
        coefs = poly.polyfit(np.array(x_filter), np.array(y_filter), 1)
        fit = poly.polyval(x, coefs)
        peptide_dict_models[peptide_class] = (coefs, fit)

    peptides_filtered = {}
    global_decoy_count = 0
    global_target_count = 0
    for peptide_hit in peptide_hit_dict:
      (psm, peptide_identification) = peptide_hit
      accessions = [ev.getProteinAccession() for ev in psm.getPeptideEvidences()]
      score = -np.log10(float(psm.getScore()))
      pass_fdr = False
      canonical_peptide = True
      for peptide_class in peptide_dict_models:
        if self.is_peptide_group(self._peptide_groups_prefix[peptide_class], accessions):
          canonical_peptide = False
          counts = score_dic[peptide_class][score]
          target_count = float(counts[0])
          decoy_count = float(counts[1])
          FDR = decoy_count / target_count
          novel_targetcount = float(counts[2])
          gamma = poly.polyval(score, peptide_dict_models[peptide_class][0])
          try:
             novelFDR = FDR * gamma * (target_count / novel_targetcount)
          except ZeroDivisionError:
            novelFDR = 10000
            # If the model raise an error because novel_targetcount = 0
          if novelFDR < self._psm_pep_class_fdr_cutoff and FDR < self._psm_pep_fdr_cutoff:
            pass_fdr = True

      if any(self._decoy_prefix in s for s in accessions):
        global_decoy_count += 1
      else:
        global_target_count += 1

      FDR = float(global_decoy_count) / float(global_target_count)
      if FDR < self._psm_pep_fdr_cutoff and canonical_peptide:
        pass_fdr = True

      if pass_fdr:
        identifier = peptide_identification.getMetaValue("spectrum_reference")
        if identifier not in peptides_filtered:
          peptide_id = peptide_identification
          peptide_id.setHits([psm])
          peptides_filtered[identifier] = peptide_id
        else:
          peptide_id = peptides_filtered[identifier]
          hits = peptide_id.getHits()
          hits.append(psm)
          peptide_id.setHits(hits)
          peptides_filtered[identifier] = peptide_id

    return list(peptides_filtered.values())

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
    idxml_parser().load(input_idxml, protein_ids, peptide_ids)
    print(peptide_ids[0])

    # Iterate over PeptideIdentification
    filtered_peptide_ids = list(
      filter(lambda peptide: (len(peptide.getHits()[0].getSequence().toUnmodifiedString()) > self._min_peptide_length),
             peptide_ids))

    print(len(filtered_peptide_ids))
    if self._peptide_class_fdr_enable:
      filtered_peptide_ids = self.compute_class_fdr(filtered_peptide_ids)
    else:
      filtered_peptide_ids = self.compute_global_fdr(filtered_peptide_ids)

    print(len(filtered_peptide_ids))

    remove_peptides_without_reference = True
    idfilter.updateProteinReferences(filtered_peptide_ids, protein_ids, remove_peptides_without_reference)
    idxml_parser().store(output_idxml, protein_ids, filtered_peptide_ids)









