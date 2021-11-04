from pyopenms import *

from pypgatk.toolbox.general import ParameterConfiguration


class OpenmsDataService(ParameterConfiguration):
  CONFIG_KEY_OPENMS_ANALYSIS = 'openms_analysis'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'

  def __init__(self, config_file, pipeline_arguments):
    super(OpenmsDataService, self).__init__(self.CONFIG_KEY_OPENMS_ANALYSIS, config_file,
                                                pipeline_arguments)

    self._decoy_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_DECOY_PREFIX]
    if self.CONFIG_DECOY_PREFIX in self.get_pipeline_parameters():
      self._decoy_prefix = self.get_pipeline_parameters()[self.CONFIG_DECOY_PREFIX]


    self._min_peptide_length = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_MIN_PEPTIDE_LENGTH]
    if self.CONFIG_MIN_PEPTIDE_LENGTH in self.get_pipeline_parameters():
      self._min_peptide_length = self.get_pipeline_parameters()[self.CONFIG_MIN_PEPTIDE_LENGTH]

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
    filtered_peptide_ids = list(filter(lambda peptide: (len(peptide.getHits()[0].getSequence().toUnmodifiedString()) > self._min_peptide_length), peptide_ids))

    remove_peptides_without_reference = True
    idfilter.removeUnreferencedProteins(protein_ids, filtered_peptide_ids, remove_peptides_without_reference)

    IdXMLFile().store(output_idxml, protein_ids, peptide_ids)







