import multiprocessing

from pandas import DataFrame
from pyopenms import IdXMLFile as idxml_parser
from pyopenms import IDFilter
from pyopenms import ModificationsDB
import json

from pypgatk.toolbox.general import ParameterConfiguration, is_peptide_group
import numpy.polynomial.polynomial as poly
import numpy as np
import pandas as pd

num_processes = multiprocessing.cpu_count()


class OpenmsDataService(ParameterConfiguration):
  CONFIG_KEY_OPENMS_ANALYSIS = 'openms_analysis'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_PEPTIDE_CLASS_PREFIX = "peptide_classes_prefix"
  CONFIG_PEPTIDE_FDR_CUTOFF = "psm_pep_fdr_cutoff"
  CONFIG_PEPTIDE_CLASS_FDR_CUTOFF = "psm_pep_class_fdr_cutoff"
  CONFIG_PEPTIDE_GROUP_PREFIX = "peptide_groups_prefix"
  CONFIG_PEPTIDE_DISABLE_CLASS_FDR = "disable_class_fdr"
  CONFIG_FILE_TYPE = "file_type"

  def __init__(self, config_file, pipeline_arguments):

    self._psm_df_index = "openms_psm_index"
    self._psm_spectrum_reference = "spectrum_reference"
    self._openms_exclude_columns = ['SpecId', 'Label', 'ScanNr', 'Peptide',
                                    'Proteins', 'FDR', 'q-val', 'class-specific-q-value',
                                    'Rank', "protein_references",
                                    'CountSequenceIsTop', 'CountSequenceCharges', 'CountSequenceIsXL',
                                    'CountSequenceIsPeptide',  # tend to overfit
                                    'NuXL:total_Morph', 'NuXL:total_HS', 'NuXL:total_MIC',  # redundant
                                    'A_136.062309999999997', 'A_330.060330000000022',
                                    'C_112.051079999999999', 'C_306.04910000000001',
                                    'G_152.057230000000004', 'G_346.055250000000001',
                                    'U_113.035089999999997', 'U_307.033110000000022',
                                    'NuXLScore_score',
                                    'NuXL:z1 mass', 'NuXL:z2 mass', 'NuXL:z3 mass', 'NuXL:z4 mass',
                                    "NuXL:NA", "NuXL:NT", "NuXL:localization_scores", "NuXL:best_localization",
                                    'NuXL:best_localization_score', "CalcMass", "NuXL:Da difference",
                                    'NuXL:XL_U', 'NuXL:XL_C', 'NuXL:XL_G', 'NuXL:XL_A'
                                    ]

    super(OpenmsDataService, self).__init__(self.CONFIG_KEY_OPENMS_ANALYSIS, config_file,
                                            pipeline_arguments)

    self._decoy_prefix = 'DECOY_'
    self._peptide_class_prefix = 'altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs'
    self._file_type = 'idxml'
    self._min_peptide_length = 5
    self._psm_pep_fdr_cutoff = 0.01
    self._psm_pep_class_fdr_cutoff = 0.01
    self._peptide_groups_prefix = {'non_canonical':['altorf','pseudo','ncRNA'],'mutations':['COSMIC','cbiomut'],'variants':['var_mut','var_rs']}
    self._peptide_class_fdr_disable = False

    self._decoy_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_DECOY_PREFIX]
    if self.CONFIG_DECOY_PREFIX in self.get_pipeline_parameters():
      self._decoy_prefix = self.get_pipeline_parameters()[self.CONFIG_DECOY_PREFIX]

    self._peptide_class_prefix = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_PEPTIDE_CLASS_PREFIX]
    if self.CONFIG_PEPTIDE_CLASS_PREFIX in self.get_pipeline_parameters():
      self._peptide_class_prefix = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_CLASS_PREFIX]

    self._file_type = self.get_default_parameters()[self.CONFIG_KEY_OPENMS_ANALYSIS][
      self.CONFIG_FILE_TYPE]
    if self.CONFIG_FILE_TYPE in self.get_pipeline_parameters():
      self._file_type = self.get_pipeline_parameters()[self.CONFIG_FILE_TYPE]

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

  @staticmethod
  def _filter_by_group(accessions, peptide_classes):
    return is_peptide_group(peptide_classes, accessions)

  @staticmethod
  def _series_get_qvalue(arr_fdr: list):
    s = pd.Series(arr_fdr)
    return s[::-1].cummin()[::-1]

  @staticmethod
  def _get_msrescore_modification(modification):
    specificity_index = modification.getTermSpecificity()
    specificity_str = modification.getTermSpecificityName(specificity_index)

    c_term = True
    n_term = True
    aa_mod = modification.getOrigin()
    if specificity_str == 'none':
      c_term = False
      n_term = False
    elif specificity_str == "Protein C-term" or  specificity_str == "C-term":
      n_term = False
      aa_mod = None
    else:
      c_term = False
      aa_mod = None

    formula = str(modification.getDiffFormula())
    return {"name": modification.getId(), "unimod_accession": modification.getUniModRecordId(), "mass_shift": modification.getDiffMonoMass(),
            "atomic_composition": formula, "amino_acid": aa_mod, "n_term": n_term, "c_term": c_term}

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

    # df_psms['class-specific-q-value'] = df['class-specific-q-value']
    df_psms = df_psms.merge(df['class-specific-q-value'], left_index=True, right_index=True, how='left')
    df_psms.loc[df_psms['class-specific-q-value'].isnull(), 'class-specific-q-value'] = df_psms['q-value']
    df_psms.sort_values("score", ascending=ascending, inplace=True)

    return df_psms

  def _generate_msrescore_file(self, input_xml: str, quant_method: str, decoy_pattern: str, output_json: str):
    """
     This function generates an msRescore configuration file for an idXML file.
    :param input_xml: Input IdXML with the peptides and PTMs
    :param output_json: Output json file
    :return:
    """

    prot_ids = []
    pep_ids = []
    idxml_parser().load(input_xml, prot_ids, pep_ids)

    # Get the Modification parameters
    modifications = []
    for protein_hit in prot_ids:
      search_params = protein_hit.getSearchParameters()
      print(" - Search params:", search_params)

      if (search_params is not None and search_params.fixed_modifications is not None):
        for mod in search_params.fixed_modifications:
          ox = ModificationsDB().getModification(mod.decode())
          modifications.append(self._get_msrescore_modification(ox))
      if (search_params is not None and search_params.variable_modifications is not None):
        for mod in search_params.variable_modifications:
          ox = ModificationsDB().getModification(mod.decode())
          modifications.append(self._get_msrescore_modification(ox))

      fragment_error = search_params.fragment_mass_tolerance

      model = 'HCD2021'
      if quant_method == 'TMT':
        model = 'TMT'

    print(modifications)

    mappings = []
    for mod in modifications:
      amino_acid = mod["amino_acid"]
      if amino_acid is None:
        amino_acid = "."
      mappings.append({"amino_acid": amino_acid, "unimod_accession": mod["unimod_accession"], "name": mod["name"]})

    config_msrescore = {}
    config_msrescore["$schema"] = "https://raw.githubusercontent.com/compomics/ms2rescore/master/ms2rescore/package_data/config_schema.json"
    config_msrescore["general"] = {"pipeline":"infer", "run_percolator":False, "id_decoy_pattern": decoy_pattern,  "log_level": "info"}
    config_msrescore["ms2pip"]  = {"model": model, "frag_error": fragment_error, "modifications": modifications}
    config_msrescore["idxml_to_rescore"] = {"modification_mapping": mappings}

    # Write the json output
    with open(output_json, "w") as write_file:
      json.dump(config_msrescore, write_file, indent=4)

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
    Compute the bayesian class FDR from manuscript (https://pubmed.ncbi.nlm.nih.gov/24200586/). From previous discussions
    with @yanfeng is clear that the Bayesian FDR should be tested for different score systems. The original paper only use the
    algorithm with MSGF+ and SpecEval score.
    TODO: This algorithm should be only use for research purpose but not for production. The method to be use in production should be _compute_class_fdr

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

    for psm_index, score_raw, accessions in zip(df_psms.index, df_psms["score"], df_psms["accessions"]):
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
    for psm_index, score_raw, accessions, qvalue in zip(df_psms.index, df_psms["score"], df_psms["accessions"],
                                                        df_psms["q-value"]):
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

    d = {self._psm_df_index: psm_index_qvalue, 'class-specific-q-value': class_qvalue}
    df = pd.DataFrame(d)
    df.set_index(self._psm_df_index)

    df_psms = df_psms.merge(df['class-specific-q-value'], left_index=True, right_index=True, how='left')
    df_psms.loc[df_psms['class-specific-q-value'].isnull(), 'class-specific-q-value'] = df_psms['q-value']
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
    df_psms = None
    if self._file_type == 'idxml':
      df_psms = self._psm_idxml_todf(input_idxml)
    # elif (self._file_type == "triqler"):
    #   df_psms = self._psms_triqler_todf(input_idxml)

    self.get_logger().info("Number of PSM in the file {} : {}".format(input_idxml, len(df_psms.index)))

    df_psms = self._compute_global_fdr(df_psms)

    if self._peptide_class_fdr_disable:
      df_psms = df_psms[df_psms['q-value'] < self._psm_pep_fdr_cutoff]
      self.get_logger().info("Number of PSM after Global FDR filtering: {}".format(len(df_psms.index)))
    else:
      df_psms = self._compute_class_fdr(df_psms)
      df_psms = df_psms[((df_psms['q-value'] < self._psm_pep_fdr_cutoff) & (
      df_psms['class-specific-q-value'] < self._psm_pep_class_fdr_cutoff))]
      self.get_logger().info("Number of PSM after Class FDR filtering: {}".format(len(df_psms.index)))

    if self._file_type == 'idxml':
      self._filter_write_idxml_with_df(df_psms, self._new_columns, input_idxml, output_idxml)
    elif self._file_type == "triqler":
      self._export_df_triqler(df_psms, output_idxml)

  @staticmethod
  def _get_psm_index(ms_run: str, spectrum_reference: str, psm_index: int):
    """
    Get the psm index for each spectrum from idXML
    :param ms_run: name of the file containing the spectrum
    :param spectrum_reference: spectrum reference from the PeptideIdentification in the idXML
    :param psm_index: Peptide Hits PSMs for each PeptideIdentification
    :return:
    """
    return ms_run + "_" + spectrum_reference + "_" + str(psm_index)

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
      pro_run_id = pro.getIdentifier()
      ms_run = []
      pro.getPrimaryMSRunPath(ms_run)
      ms_run = [n.decode() for n in ms_run]
      protein_dic[pro_run_id] = "_".join(ms_run)

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
        key = self._get_psm_index(ms_run_acc, specid, psmid)
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

  @staticmethod
  def _get_ptm_str(psm):
    """
    This function converts a Peptide Hit in idXML into a PTM modification format as requested by DeepLC
    see documentation of the PTMs (https://github.com/compomics/DeepLC#input-files)
    :param psm: Peptide Hit
    :return: modification string position and name of the PTM
    """
    sequence = psm.getSequence()
    mod_str = ""
    if(sequence.hasNTerminalModification()):
       mod_str = mod_str + "0|" + sequence.getNTerminalModificationName()
    if (sequence.hasCTerminalModification()):
      if len(mod_str) > 0:
        mod_str = mod_str + "|"
        mod_str = mod_str + "1|" + sequence.getCTerminalModificationName()

    i = 1
    for aa in sequence:
      mod = aa.getModificationName()
      if(len(mod) > 0):
        if len(mod_str) > 0:
          mod_str = mod_str + "|"
        mod_str = mod_str + str(i) + "|" + mod
      i = i + 1

    return mod_str

  def _psm_idxml_todf(self, input_file: str):
    """
    This function converts an idXML file into a pandas dataframe
    :param input_file: input idXML file
    :return:
    """

    prot_ids = []
    pep_ids = []
    idxml_parser().load(input_file, prot_ids, pep_ids)

    protein_dic = {}
    for pro in prot_ids:
      pro_run_id = pro.getIdentifier()
      ms_run = []
      pro.getPrimaryMSRunPath(ms_run)
      ms_run = [n.decode() for n in ms_run]
      protein_dic[pro_run_id] = "_".join(ms_run)

    meta_value_keys = []
    rows = []
    all_columns = []

    for peptide_id in pep_ids:
      spectrum_id = peptide_id.getMetaValue(self._psm_spectrum_reference)
      scan_nr = spectrum_id[spectrum_id.rfind('=') + 1:]
      psm_index = 1
      pep_acc = peptide_id.getIdentifier()
      if pep_acc not in protein_dic:
        raise ValueError(
          "The reference peptide in the PeptideIdentification {}--{} can't be found in the ProteinIdentifications".format(
            pep_acc, spectrum_id))

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
        rt = peptide_id.getRT()

        ptm_str = self._get_ptm_str(psm = h)

        if len(meta_value_keys) == 0:
          h.getKeys(meta_value_keys)
          meta_value_keys = [x.decode() for x in meta_value_keys if not (
            "target_decoy" in x.decode() or "spectrum_reference" in x.decode() or "rank" in x.decode() or x.decode() in self._openms_exclude_columns)]
          all_columns = [self._psm_df_index, "target", "scanNr", "charge", "mz", "rt", "peptide", "unmodified_peptide", "mod_str",
                         "peptide_length", "accessions", "score", "is_higher_score_better"] + meta_value_keys

        df_psm_index = self._get_psm_index(ms_run_acc, spectrum_id, psm_index)
        row = [df_psm_index, label, scan_nr, charge, peptide_id.getMZ(), rt, sequence, unmodified_sequence, ptm_str,
               str(len(unmodified_sequence)), accessions, score, order]
        # scores in meta values
        for k in meta_value_keys:
          if not (
            "target_decoy" in k or "spectrum_reference" in k or "rank" in k or k in self._openms_exclude_columns):  # don't add them twice
            s = h.getMetaValue(k)
            if isinstance(s, bytes):
              s = s.decode()
            row.append(s)
        rows.append(row)
        psm_index += 1
    df = pd.DataFrame(rows, columns=all_columns)

    df = self._str_to_int(df)
    df.set_index(self._psm_df_index, inplace=True)

    return df

  def _psms_triqler_todf(self, input_file: str):
    """
    Parse triqler file .
    :param input_file:
    :return:
    """
    all_columns = [self._psm_df_index, "run", "condition", "charge", "score", "intensity", "peptide", "accessions",
                   "is_higher_score_better"]
    rows = []
    first_line = 0
    with open(input_file, 'r') as source:
      for line in source:
        line = line.rstrip()
        if first_line != 0:
          a = line.split('\t')
          a = ["_".join(a), *a, True]
          rows.append(a)
        first_line += 1

    df = pd.DataFrame(rows, columns=all_columns)

    print(len(df.index))
    # with concurrent.futures.ProcessPoolExecutor(num_processes) as pool:
    #   df['target'] = pd.DataFrame(pool.map(df['accessions'].map(lambda x: is_peptide_decoy(x, decoy_prefix))))  # With a progressbar
    # print(len(df.index))
    df["target"] = df['accessions'].map(lambda x: not (self._decoy_prefix in x))

    df = self._str_to_int(df)
    df.set_index(self._psm_df_index, inplace=True)
    return df

  @staticmethod
  def _export_df_triqler(df_psms: DataFrame, output_file: str):
    """
    Export dataframe to triqler format
    :param df_psms:  dataframe containing the psms
    :param output_file: output triqler file
    :return:
    """
    result_df = df_psms[["run", "condition", "charge", "score", "intensity", "peptide", "accessions"]]
    result_df.rename(columns={"score": "searchScore", "accessions": "proteins"}, errors="raise")
    result_df.to_csv(output_file, sep='\t', index=False, header=True)

  def _generate_deepLC_file(self, input_xml: str, output_deepLC: str, decoy_pattern: str,
                                        peptide_class_prefix: str, novel_peptides: bool):
    peptides = self._psm_idxml_todf(input_file = input_xml)
    print(peptides.head())

    # remove the decoy peptides
    peptides = peptides.loc[peptides['target'] == 1]

    if peptide_class_prefix is None:
      peptide_class_prefix = self._peptide_class_prefix.split(",")
    else:
      peptide_class_prefix = peptide_class_prefix.split(",")

    peptides = peptides[peptides['accessions'].apply(lambda x: not (any(self._decoy_prefix in s for s in x)))]

    if novel_peptides:
      peptides = peptides[peptides['accessions'].apply(lambda x: self._filter_by_group(x, peptide_class_prefix))]
      peptides = peptides[["unmodified_peptide", "mod_str"]]
      peptides = peptides.rename(columns={'unmodified_peptide': 'seq', 'mod_str': 'modifications'})
    else:
      peptides = peptides[peptides['accessions'].apply(lambda x: not self._filter_by_group(x, peptide_class_prefix))]
      peptides = peptides[["unmodified_peptide", "mod_str", "rt"]]
      peptides = peptides.rename(columns={'unmodified_peptide': 'seq', 'mod_str': 'modifications', 'rt': 'tr'})

    peptides.to_csv(output_deepLC, sep=',', index=False, header=True)

    print(peptides.head())
