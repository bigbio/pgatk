from pandas import DataFrame
from pyopenms import IdXMLFile as idxml_parser
from pyopenms import IDFilter

from pypgatk.toolbox.general import ParameterConfiguration, is_peptide_group
import numpy.polynomial.polynomial as poly
import numpy as np
import pandas as pd


class OpenmsDataService(ParameterConfiguration):
  CONFIG_KEY_OPENMS_ANALYSIS = 'openms_analysis'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_PEPTIDE_CLASS_PREFIX = "peptide_classes_prefix"
  CONFIG_PEPTIDE_FDR_CUTOFF = "psm_pep_fdr_cutoff"
  CONFIG_PEPTIDE_CLASS_FDR_CUTOFF = "psm_pep_class_fdr_cutoff"
  CONFIG_PEPTIDE_GROUP_PREFIX = "peptide_groups_prefix"
  CONFIG_PEPTIDE_DISABLE_CLASS_FDR = "disable_class_fdr"
  CONFIG_PEPTIDE_DISABLE_BAYESIAN_FDR = "disable_bayesian_class_fdr"

  def __init__(self, config_file, pipeline_arguments):

    self._psm_df_index = "openms_psm_index"
    self._psm_spectrum_reference = "spectrum_reference"
    self._openms_exclude_columns = ['SpecId', 'Label', 'ScanNr', 'Peptide',
            'Proteins', 'FDR', 'q-val', 'class-specific-q-value',
            'Rank', "protein_references",
            'CountSequenceIsTop', 'CountSequenceCharges', 'CountSequenceIsXL', 'CountSequenceIsPeptide', # tend to overfit
            'NuXL:total_Morph', 'NuXL:total_HS', 'NuXL:total_MIC', # redundant
            'A_136.062309999999997', 'A_330.060330000000022',
            'C_112.051079999999999', 'C_306.04910000000001',
            'G_152.057230000000004', 'G_346.055250000000001',
            'U_113.035089999999997', 'U_307.033110000000022',
            'NuXLScore_score',
            'NuXL:z1 mass', 'NuXL:z2 mass', 'NuXL:z3 mass', 'NuXL:z4 mass',
            "NuXL:NA", "NuXL:NT", "NuXL:localization_scores", "NuXL:best_localization",
            'NuXL:best_localization_score', "CalcMass", "NuXL:Da difference",
             'NuXL:XL_U', 'NuXL:XL_C', 'NuXL:XL_G','NuXL:XL_A'
            ]

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

    self._peptide_class_fdr_disable = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_DISABLE_CLASS_FDR]
    if self.CONFIG_PEPTIDE_DISABLE_CLASS_FDR in self.get_pipeline_parameters():
      self._peptide_class_fdr_disable = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_DISABLE_CLASS_FDR]

    self._bayesian_class_fdr_disable = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_DISABLE_BAYESIAN_FDR]
    if self.CONFIG_PEPTIDE_DISABLE_BAYESIAN_FDR in self.get_pipeline_parameters():
      self._bayesian_class_fdr_disable = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_DISABLE_BAYESIAN_FDR]

  @staticmethod
  def _filter_by_group(accessions, peptide_classes):
    return is_peptide_group(peptide_classes, accessions)

  def _compute_class_fdr(self, df_psms: DataFrame):

    # Get the order of the score
    orders = df_psms["is_higher_score_better"].unique()
    if len(orders) != 1:
      raise ValueError("The Global FDR error do not support multiple orders for scores")
    ascending = (orders[0] == False)

    ls = []
    for c in self._peptide_groups_prefix:
      # split the dataframe and save the subset
      currClass = df_psms[
        df_psms['accessions'].apply(lambda x: self._filter_by_group(x, self._peptide_groups_prefix[c]))]
      ls.append(currClass)

      # calculate class-specific q-value
      currClass.sort_values("score", ascending=ascending, inplace=True)
      FDR = (range(1, len(currClass["target"]) + 1) / currClass["target"].cumsum()) - 1
      currClass['class-specific-q-value'] = FDR[::-1].cummin()[::-1]
    df = pd.concat(ls)

    df_psms['class-specific-q-value'] = df['class-specific-q-value']
    df_psms['class-specific-q-value'][df_psms['class-specific-q-value'].isnull()] = df_psms['q-value']
    df_psms.sort_values("score", ascending=ascending, inplace=True)

    return df_psms

  @staticmethod
  def _compute_global_fdr(df_psms: DataFrame):
    """
    Compute the global FDR and filter peptides
    :param df_psms: list of peptide identifications
    :return:  filtered peptides
    """

    # Get the order of the score
    orders = df_psms["is_higher_score_better"].unique()
    if len(orders) != 1:
      raise ValueError("The Global FDR error do not support multiple orders for scores")
    ascending = (orders[0] == False)

    df_psms.sort_values("score", ascending=ascending, inplace=True)
    df_psms['FDR'] = (range(1, len(df_psms) + 1) / df_psms['target'].cumsum()) - 1
    df_psms['q-value'] = df_psms['FDR'][::-1].cummin()[::-1]

    df_psms.sort_values("score", ascending=ascending, inplace=True)

    return df_psms

  def _compute_bayesian_class_fdr(self, df_psms: DataFrame):
    """
    Compute the global FDR and filter peptides
    :param df_psms: list of peptide identifications
    :return:  filtered peptides
    """

    score_dic = {}
    decoy_dic = {}

    # A dictionary of peptide classes with number of decoys and targets per class
    novel_target_decoy_count = {}
    for peptide_class in self._peptide_groups_prefix:
      novel_target_decoy_count[peptide_class] = (0, 0)
      decoy_dic[peptide_class] = {}
      score_dic[peptide_class] = {}

    global_decoy_count = 0
    global_target_count = 0

    # Get the order of the score
    orders = df_psms["is_higher_score_better"].unique()
    if len(orders) != 1:
      raise ValueError("The Global FDR error do not support multiple orders for scores")
    ascending = (orders[0] == False)

    df_psms.sort_values("score", ascending=ascending, inplace=True)

    for psm_index, score_raw, accessions  in zip(df_psms.index, df_psms["score"], df_psms["accessions"]):
      score = -np.log10(float(score_raw))
      if any(self._decoy_prefix in s for s in accessions):
        global_decoy_count += 1
        for peptide_group in self._peptide_groups_prefix:
          if is_peptide_group(self._peptide_groups_prefix[peptide_group], accessions):
            (target, decoy) = novel_target_decoy_count[peptide_group]
            decoy += 1
            novel_target_decoy_count[peptide_group] = (target, decoy)
            decoy_dic[peptide_group][score] = [decoy, global_decoy_count]
      else:
        global_target_count += 1
        for peptide_group in self._peptide_groups_prefix:
          if is_peptide_group(self._peptide_groups_prefix[peptide_group], accessions):
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

    global_decoy_count = 0
    global_target_count = 0

    psm_index_qvalue = []
    class_qvalue = []
    for psm_index, score_raw, accessions, qvalue  in zip(df_psms.index, df_psms["score"], df_psms["accessions"], df_psms["q-value"]):
      score = -np.log10(float(score_raw))
      canonical_peptide = True
      for peptide_class in peptide_dict_models:
        if is_peptide_group(self._peptide_groups_prefix[peptide_class], accessions):
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
          psm_index_qvalue.append(psm_index)
          class_qvalue.append(novelFDR)

      if any(self._decoy_prefix in s for s in accessions):
        global_decoy_count += 1
      else:
        global_target_count += 1

      if canonical_peptide:
        psm_index_qvalue.append(psm_index)
        class_qvalue.append(qvalue)

    d = {self._psm_df_index:psm_index_qvalue,'class-specific-q-value':class_qvalue}
    df = pd.DataFrame(d)
    df.set_index(self._psm_df_index)

    df_psms['class-specific-q-value'] = df['class-specific-q-value']
    df_psms['class-specific-q-value'][df_psms['class-specific-q-value'].isnull()] = df_psms['q-value']
    df_psms.sort_values("score", ascending=ascending, inplace=True)

    return df_psms

  def filter_peptide_class_fdr(self, input_idxml: str, output_idxml: str):
    """
    Filter peptides by Class FDR in idXML files
    :param input_idxml:  input idxml
    :param output_idxml: output idxml
    :return:
    """

    self._new_columns = []

    df_psms = self._psm_idxml_todf(input_idxml)
    print(df_psms.size)

    df_psms = self._compute_global_fdr(df_psms)

    if self._peptide_class_fdr_disable:
      df_psms = df_psms[df_psms['q-value'] < self._psm_pep_fdr_cutoff]
    elif(self._bayesian_class_fdr_disable):
      df_psms = self._compute_class_fdr(df_psms)
      df_psms = df_psms[((df_psms['q-value'] < self._psm_pep_fdr_cutoff) & (df_psms['class-specific-q-value'] < self._psm_pep_class_fdr_cutoff))]
    else:
      df_psms = self._compute_bayesian_class_fdr(df_psms)
      df_psms = df_psms[((df_psms['q-value'] < self._psm_pep_fdr_cutoff) & (df_psms[
        'class-specific-q-value'] < self._psm_pep_class_fdr_cutoff))]

    print(df_psms.size)

    self._filter_write_idxml_with_df(df_psms, self._new_columns, input_idxml, output_idxml)

  @staticmethod
  def _get_psm_index(ms_run: str, spectrum_reference: str, psm_index: int):
    """
    Get the psm index for each spectrum from idXML
    :param ms_run: name of the file containing the spectrum
    :param spectrum_reference: spectrum reference from the PeptideIdentification in the idXML
    :param psm_index: Peptide Hits PSMs for each PeptideIdentification
    :return:
    """
    return spectrum_reference + "_" + str(psm_index)

  def _filter_write_idxml_with_df(self, df: DataFrame, new_columns: list, input_file: str, output_file: str):
    """
    This function allows to filter an input idXML with a set of PMSs in a panda daframe (df).
    :param input_file: original idXML file
    :param df: the PSMs after filtering or processing
    :param output_file: output idXML
    :param new_columns: columns to be added to psms peptide hits
    :return:
    """

    prot_ids = []
    pep_ids = []
    idxml_parser().load(input_file, prot_ids, pep_ids)

    new_pep_ids = []

    protein_dic = {}
    for pro in prot_ids:
      pro_acc = pro.getIdentifier()
      ms_run = []
      pro.getPrimaryMSRunPath(ms_run)
      ms_run = [n.decode() for n in ms_run]
      protein_dic[pro_acc] = "_".join(ms_run)

    for peptide_id in pep_ids:
      hits = peptide_id.getHits()
      psmid = 1
      specid = peptide_id.getMetaValue(self._psm_spectrum_reference)

      pep_acc = peptide_id.getIdentifier()
      if pep_acc not in protein_dic:
        raise ValueError(
          "The reference peptide in the PeptideIdentification {}--{} can't be found in the ProteinIdentifications".format(
            pep_acc, specid))

      ms_run_acc = protein_dic[pep_acc]

      new_hits = []
      for h in hits:
        key = self._get_psm_index(ms_run, specid, psmid)
        if key in df.index:
          for col in new_columns:
            value = df.at[key, col]
            if value is not None:
              h.setMetaValue(col, value)
          new_hits.append(h)
        psmid += 1
      peptide_id.setHits(new_hits)
      new_pep_ids.append(peptide_id)

    remove_peptides_without_reference = True
    idfilter = IDFilter()
    idfilter.removeEmptyIdentifications(new_pep_ids)
    idfilter.updateProteinReferences(new_pep_ids, prot_ids, remove_peptides_without_reference)
    idxml_parser().store(output_file, prot_ids, pep_ids)

  @staticmethod
  def _str_to_int(df: DataFrame):
    """
    Convert string in dataframe to Flot values
    :param df: DataFrame
    :return: Dataframe with all values in Flot
    """
    for col in df:
      try:
        if all([i.is_integer() for i in df[col]]):
          df[col] = [int(x) for x in df[col]]
      except AttributeError:
        continue
    return df

  def _psm_idxml_todf(self, input_file: str):

    prot_ids = []
    pep_ids = []
    idxml_parser().load(input_file, prot_ids, pep_ids)

    protein_dic = {}
    for pro in prot_ids:
      pro_acc = pro.getIdentifier()
      ms_run = []
      pro.getPrimaryMSRunPath(ms_run)
      ms_run = [n.decode() for n in ms_run]
      protein_dic[pro_acc] = "_".join(ms_run)
      print("{} -- {}".format(pro_acc,"_".join(ms_run)))

    meta_value_keys = []
    rows = []
    all_columns = []

    for peptide_id in pep_ids:
      spectrum_id = peptide_id.getMetaValue(self._psm_spectrum_reference)
      scan_nr = spectrum_id[spectrum_id.rfind('=') + 1:]
      psm_index = 1
      pep_acc = peptide_id.getIdentifier()
      if pep_acc not in protein_dic:
        raise ValueError("The reference peptide in the PeptideIdentification {}--{} can't be found in the ProteinIdentifications".format(pep_acc, spectrum_id))

      ms_run_acc = protein_dic[pep_acc]
      hits = peptide_id.getHits()
      order = peptide_id.isHigherScoreBetter()

      for h in hits:
        charge = h.getCharge()
        if "target" in h.getMetaValue("target_decoy"):
          label = 1
        else:
          label = 0
        sequence = h.getSequence().toString()
        unmodified_sequence = h.getSequence().toUnmodifiedString()
        accessions = [ev.getProteinAccession() for ev in h.getPeptideEvidences()]
        score = h.getScore()

        if len(meta_value_keys) == 0:
          h.getKeys(meta_value_keys)
          meta_value_keys = [x.decode() for x in meta_value_keys if not ("target_decoy" in x.decode() or "spectrum_reference" in x.decode() or "rank" in x.decode() or x.decode() in self._openms_exclude_columns)]
          all_columns = [self._psm_df_index, "target", "scanNr", "charge", "mz", "peptide", "unmodified_peptide", "peptide_length", "accessions", "score", "is_higher_score_better"] + meta_value_keys

        df_psm_index = self._get_psm_index(ms_run_acc, spectrum_id, psm_index)
        row = [df_psm_index, label, scan_nr, charge, peptide_id.getMZ(), sequence, unmodified_sequence, str(len(unmodified_sequence)), accessions, score, order]
        # scores in meta values
        for k in meta_value_keys:
          if not (
            "target_decoy" in k or "spectrum_reference" in k or "rank" in k or k in self._openms_exclude_columns):  # don't add them twice
            s = h.getMetaValue(k)
            if isinstance(s,bytes):
              s = s.decode()
            row.append(s)
        rows.append(row)
        psm_index += 1
    df = pd.DataFrame(rows, columns=all_columns)

    print(df.head())

    df = self._str_to_int(df)
    df.set_index(self._psm_df_index, inplace=True)

    return df









