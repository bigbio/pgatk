import gffutils
import vcf
from Bio import SeqIO
from Bio.Seq import Seq

from pypgatk.toolbox.general import ParameterConfiguration


class EnsemblDataService(ParameterConfiguration):
  CONFIG_KEY_VCF = "ensembl_translation"
  INPUT_FASTA = "input_fasta"
  TRANSLATION_TABLE = "translation_table"
  MITO_TRANSLATION_TABLE = "mito_translation_table"
  HEADER_VAR_PREFIX = "var_prefix"
  REPORT_REFERENCE_SEQ = "report_ref_seq"
  PROTEIN_DB_OUTPUT = "proteindb_output_file"
  ANNOTATION_FIELD_NAME = "annotation_field_name"
  AF_FIELD = "af_field"
  AF_THRESHOLD = "af_threshold"
  TRANSCRIPT_INDEX = "transcript_index"
  CONSEQUENCE_INDEX = "consequence_index"
  EXCLUDE_BIOTYPES = "exclude_biotypes"
  EXCLUDE_CONSEQUENCES = "exclude_consequences"
  SKIP_INCLUDING_ALL_CDS = "skip_including_all_cds"
  INCLUDE_BIOTYPES = "include_biotypes"
  INCLUDE_CONSEQUENCES = "include_consequences"
  BIOTYPE_STR = "biotype_str"
  SKIP_INCLUDING_ALL_CDSS = "skip_including_all_CDSs"
  CONFIG_KEY_DATA = "ensembl_translation"
  NUM_ORFS = "num_orfs"
  NUM_ORFS_COMPLEMENT = "num_orfs_complement"
  EXPRESSION_STR = "expression_str"
  EXPRESSION_THRESH = "expression_thresh"
  IGNORE_FILTERS = "ignore_filters"
  ACCEPTED_FILTERS = "accepted_filters"

  def __init__(self, config_file, pipeline_arguments):
    """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
    super(EnsemblDataService, self).__init__(self.CONFIG_KEY_DATA, config_file,
                                             pipeline_arguments)

    self._proteindb_output = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.PROTEIN_DB_OUTPUT]
    if self.PROTEIN_DB_OUTPUT in self.get_pipeline_parameters():
      self._proteindb_output = self.get_pipeline_parameters()[self.PROTEIN_DB_OUTPUT]

    self._translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.TRANSLATION_TABLE]
    if self.TRANSLATION_TABLE in self.get_pipeline_parameters():
      self._translation_table = self.get_pipeline_parameters()[self.TRANSLATION_TABLE]

    self._mito_translation_table = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.MITO_TRANSLATION_TABLE]
    if self.MITO_TRANSLATION_TABLE in self.get_pipeline_parameters():
      self._mito_translation_table = self.get_pipeline_parameters()[self.MITO_TRANSLATION_TABLE]

    self._header_var_prefix = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.HEADER_VAR_PREFIX]
    if self.HEADER_VAR_PREFIX in self.get_pipeline_parameters():
      self._header_var_prefix = self.get_pipeline_parameters()[self.HEADER_VAR_PREFIX]

    self._report_reference_seq = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.REPORT_REFERENCE_SEQ]
    if self.REPORT_REFERENCE_SEQ in self.get_pipeline_parameters():
      self._report_reference_seq = self.get_pipeline_parameters()[self.REPORT_REFERENCE_SEQ]

    self._annotation_field_name = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.ANNOTATION_FIELD_NAME]
    if self.ANNOTATION_FIELD_NAME in self.get_pipeline_parameters():
      self._annotation_field_name = self.get_pipeline_parameters()[self.ANNOTATION_FIELD_NAME]

    self._af_field = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.AF_FIELD]
    if self.AF_FIELD in self.get_pipeline_parameters():
      self._af_field = self.get_pipeline_parameters()[self.AF_FIELD]

    self._af_threshold = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.AF_THRESHOLD]
    if self.AF_THRESHOLD in self.get_pipeline_parameters():
      self._af_threshold = self.get_pipeline_parameters()[self.AF_THRESHOLD]

    self._transcript_index = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.TRANSCRIPT_INDEX]
    if self.TRANSCRIPT_INDEX in self.get_pipeline_parameters():
      self._transcript_index = self.get_pipeline_parameters()[self.TRANSCRIPT_INDEX]

    self._consequence_index = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.CONSEQUENCE_INDEX]
    if self.CONSEQUENCE_INDEX in self.get_pipeline_parameters():
      self._consequence_index = self.get_pipeline_parameters()[self.CONSEQUENCE_INDEX]

    self._exclude_biotypes = self.get_multiple_options(
      self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.EXCLUDE_BIOTYPES])
    if self.EXCLUDE_BIOTYPES in self.get_pipeline_parameters():
      self._exclude_biotypes = self.get_multiple_options(self.get_pipeline_parameters()[self.EXCLUDE_BIOTYPES])

    self._exclude_consequences = self.get_multiple_options(
      self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.EXCLUDE_CONSEQUENCES])
    if self.EXCLUDE_CONSEQUENCES in self.get_pipeline_parameters():
      self._exclude_consequences = self.get_multiple_options(
        self.get_pipeline_parameters()[self.EXCLUDE_CONSEQUENCES])

    self._skip_including_all_cds = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.SKIP_INCLUDING_ALL_CDS]
    if self.SKIP_INCLUDING_ALL_CDS in self.get_pipeline_parameters():
      self._skip_including_all_cds = self.get_pipeline_parameters()[self.SKIP_INCLUDING_ALL_CDS]

    self._include_biotypes = self.get_multiple_options(
      self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.INCLUDE_BIOTYPES])
    if self.INCLUDE_BIOTYPES in self.get_pipeline_parameters():
      self._include_biotypes = self.get_multiple_options(self.get_pipeline_parameters()[self.INCLUDE_BIOTYPES])

    self._include_consequences = self.get_multiple_options(
      self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.INCLUDE_CONSEQUENCES])
    if self.INCLUDE_CONSEQUENCES in self.get_pipeline_parameters():
      self._include_consequences = self.get_multiple_options(
        self.get_pipeline_parameters()[self.INCLUDE_CONSEQUENCES])

    self._biotype_str = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.BIOTYPE_STR]
    if self.BIOTYPE_STR in self.get_pipeline_parameters():
      self._biotype_str = self.get_pipeline_parameters()[self.BIOTYPE_STR]

    self._num_orfs = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.NUM_ORFS]
    if self.NUM_ORFS in self.get_pipeline_parameters():
      self._num_orfs = self.get_pipeline_parameters()[self.NUM_ORFS]

    self._num_orfs_complement = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.NUM_ORFS_COMPLEMENT]
    if self.NUM_ORFS_COMPLEMENT in self.get_pipeline_parameters():
      self._num_orfs_complement = self.get_pipeline_parameters()[self.NUM_ORFS_COMPLEMENT]

    self._expression_str = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.EXPRESSION_STR]
    if self.EXPRESSION_STR in self.get_pipeline_parameters():
      self._expression_str = self.get_pipeline_parameters()[self.EXPRESSION_STR]

    self._expression_thresh = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.EXPRESSION_THRESH]
    if self.EXPRESSION_THRESH in self.get_pipeline_parameters():
      self._expression_thresh = self.get_pipeline_parameters()[self.EXPRESSION_THRESH]

    self._ignore_filters = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][
      self.IGNORE_FILTERS]
    if self.IGNORE_FILTERS in self.get_pipeline_parameters():
      self._ignore_filters = self.get_pipeline_parameters()[self.IGNORE_FILTERS]

    self._accepted_filters = self.get_multiple_options(
      self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_KEY_VCF][self.ACCEPTED_FILTERS])
    if self.ACCEPTED_FILTERS in self.get_pipeline_parameters():
      self._accepted_filters = self.get_multiple_options(
        self.get_pipeline_parameters()[self.ACCEPTED_FILTERS])

  def three_frame_translation(self, input_fasta):
    """
        This function translate a transcriptome into a 3'frame translation protein sequence database
        :param input_fasta: fasta input file
        :return:
        """

    input_handle = open(input_fasta, 'r')
    output_handle = open(self._proteindb_output, 'w')

    for record in SeqIO.parse(input_handle, 'fasta'):
      seq = record.seq
      RF1 = seq.translate(table=self._translation_table)
      RF2 = seq[1::].translate(table=self._translation_table)
      RF3 = seq[2::].translate(table=self._translation_table)

      if record.id == "":
        print("skip entries without id", record.description)
        continue
      output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF1', RF1))
      output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF2', RF2))
      output_handle.write("%s\n%s\n" % ('>' + record.id + '_RF3', RF3))

  @staticmethod
  def get_multiple_options(options_str: str):
    """
        This method takes an String like option1, option2, ... and produce and array [option1, option2,... ]
        :param options_str:
        :return: Array
        """
    return list(map(lambda x: x.strip(), options_str.split(",")))

  @staticmethod
  def check_overlap(var_start, var_end, features_info=None):
    """
        This function returns true when the variant overlaps any of the features
        :param var_start: Start location
        :param var_end: End location
        :param features_info: Feature information (default = [[0, 1, 'type']])
        :return:
        """
    if features_info is None:
      features_info = [[0, 1, 'type']]
    if var_start == -1:
      return True
    # check if the var overlaps any of the features
    for feature_pos in features_info:
      pep_start = feature_pos[0]
      pep_end = feature_pos[1]
      if var_start <= pep_start <= var_end:  # fully contained or partial overlap from the end
        return True
      elif var_start <= pep_end <= var_end:  # partial overlap in the begining
        return True
      elif pep_start <= var_start and pep_end >= var_end:  # fully covered
        return True
    return False

  @staticmethod
  def get_altseq(ref_seq, ref_allele, var_allele, var_pos, strand, features_info, cds_info=None):
    """
    The given sequence in the fasta file represents all exons of the transcript combined.
    for protein coding genes, the CDS is specified therefore the sequence position has to be
    calculated based on the CDS's positions
    However, for non-protein coding genes, the whole sequence is used
    :param ref_seq:
    :param ref_allele:
    :param var_allele:
    :param var_pos:
    :param strand:
    :param features_info:
    :param cds_info:
    :return:
    """
    if cds_info is None:
      cds_info = []
    alt_seq = ""
    if len(cds_info) == 2:
      start_coding_index = cds_info[0] - 1  # it should be index not pos
      stop_coding_index = cds_info[1]  # get end position of the  last cds
    else:
      start_coding_index = 0
      total_len = 0
      for x in features_info:
        total_len += x[1] - x[0] + 1
      stop_coding_index = total_len  # the features are sorted by end therefroe the end pos of the last item is the last coding nc

    if strand == '-':  # ge the correct orientation, because exons are oredered based on their position
      ref_seq = ref_seq[
                ::-1]  # in order to calculate from the first base of the first feature (sorted by genomic coordinates)
      ref_allele = ref_allele.complement()  # the reverse will be done on return
      var_allele = var_allele.complement()  # the reverse will be done on return

    ref_seq = ref_seq[
              start_coding_index:stop_coding_index]  # just keep the coding regions (mostly effective in case of protein-coding genes)
    nc_index = 0
    if len(ref_allele) == len(var_allele) or ref_allele[0] == var_allele[0]:
      for feature in features_info:  # for every exon, cds or stop codon
        if var_pos in range(feature[0], feature[
                                          1] + 1):  # get index of the var relative to the position of the overlapping feature in the coding region
          var_index_in_cds = nc_index + (var_pos - feature[0])
          # modify the coding reference sequence accoding to the var_allele
          c = len(ref_allele)
          alt_seq = ref_seq[0:var_index_in_cds] + var_allele + ref_seq[
                                                               var_index_in_cds + c::]  # variant and ref strand??
          if strand == '-':
            return ref_seq[::-1], alt_seq[::-1]
          else:
            return ref_seq, alt_seq

        nc_index += (feature[1] - feature[0] + 1)

    return ref_seq, alt_seq

  @staticmethod
  def parse_gtf(gene_annotations_gtf, gtf_db_file):
    """
        Convert GTF file into a FeatureDB
        :param gene_annotations_gtf:
        :param gtf_db_file:
        :return:
        """
    try:
      gffutils.create_db(gene_annotations_gtf, gtf_db_file, merge_strategy="create_unique",
                         keep_order=True, disable_infer_transcripts=True, disable_infer_genes=True,
                         verbose=True,
                         force=False)
    except Exception as e:  # already exists
      print("Databae already exists" + str(e), gtf_db_file)

    db = gffutils.FeatureDB(gtf_db_file)
    return db

  @staticmethod
  def get_features(db, feature_id, biotype_str, feature_types=None):
    """
        Get chr, genomic positions, strand and biotype for feature_id
        also genomic positions for all its elements (exons/cds&start_codon)
        :param db:
        :param feature_id:
        :param biotype_str:
        :param feature_types:
        :return:
        """
    if feature_types is None:
      feature_types = ['exon']
    try:
      feature = db[feature_id]
    except gffutils.exceptions.FeatureNotFoundError:  # remove version number from the ID
      try:
        feature = db[feature_id.split('.')[0]]
      except gffutils.exceptions.FeatureNotFoundError:
        print("""Feature {} found in fasta file but not in gtf file. Check that the fasta file and the gtf files match.
                        A common issue is when the fasta file have chromosome patches but not the gtf""".format(
          feature_id))
        return None, None, None, None
    coding_features = []
    features = db.children(feature, featuretype=feature_types, order_by='end')
    for f in features:
      f_type = f.featuretype
      coding_features.append([f.start, f.end, f_type])
    return feature.chrom, feature.strand, coding_features, feature.attributes[biotype_str][0]

  @staticmethod
  def get_orfs_vcf(ref_seq: str, alt_seq: str, translation_table: int, num_orfs=1):
    """
        Translate the coding_ref and the coding_alt into ORFs
        :param ref_seq:
        :param alt_seq:
        :param translation_table:
        :param num_orfs:
        :return:
        """

    ref_orfs = []
    alt_orfs = []
    for n in range(0, num_orfs):
      ref_orfs.append(ref_seq[n::].translate(translation_table))
      alt_orfs.append(alt_seq[n::].translate(translation_table))

    return ref_orfs, alt_orfs

  @staticmethod
  def get_orfs_dna(ref_seq: str, translation_table: int, num_orfs: int, num_orfs_complement: int, to_stop: bool):
    """translate the coding_ref into ORFs"""

    ref_orfs = []
    for n in range(0, num_orfs):
      ref_orfs.append(ref_seq[n::].translate(translation_table, to_stop=to_stop))

    rev_ref_seq = ref_seq.reverse_complement()
    for n in range(0, num_orfs_complement):
      ref_orfs.append(rev_ref_seq[n::].translate(translation_table, to_stop=to_stop))

    return ref_orfs

  def dnaseq_to_proteindb(self, input_fasta):
    """
        translates DNA sequences to protein sequences
        :param input_fasta:
        :return:
        """

    seq_dict = SeqIO.index(input_fasta, "fasta")

    with open(self._proteindb_output, 'w') as prots_fn:
      for record_id in seq_dict.keys():

        ref_seq = seq_dict[record_id].seq  # get the seq and desc for the record from the fasta of the gtf
        desc = str(seq_dict[record_id].description)

        key_values = {}  # extract key=value in the desc into a dict
        sep = ' '
        if '|' in desc:
          sep = '|'
        for value in desc.split(sep):
          if value.split('=')[0] == 'cds' or value.split(':')[0] == 'cds':
            value.replace('cds', 'CDS')
          if '=' in value:
            key_values[value.split('=')[0]] = value.split('=')[1]
          elif ':' in value:
            key_values[value.split(':')[0]] = value.split(':')[1]
          elif value.split('=')[
            0] == 'CDS':  # when only it is specified to be a CDS, it means the whole sequence to be used
            key_values[value.split('=')[0]] = '{}-{}'.format(1, len(ref_seq))

        feature_biotype = ""
        if self._biotype_str:
          try:
            feature_biotype = key_values[self._biotype_str]
          except KeyError:
            msg = "Biotype info was not found in the header using {} for record {} {}".format(self._biotype_str,
                                                                                              record_id, desc)
            self.get_logger().debug(msg)

        # only include features that have the specified biotypes or they have CDSs info
        if 'CDS' in key_values.keys() and (
          not self._skip_including_all_cds or 'altORFs' in self._include_biotypes):
          pass
        elif self._biotype_str and (feature_biotype == "" or (feature_biotype in self._exclude_biotypes or
                                                              (
                                                                feature_biotype not in self._include_biotypes and self._include_biotypes != [
                                                                'all']))):
          continue

        # check wether to filter on expression and if it passes
        if self._expression_str:
          try:
            if float(key_values[self._expression_str]) < self._expression_thresh:
              continue
          except KeyError:
            msg = "Expression information not found in the fasta header with expression_str: {} for record {} {}".format(
              self._expression_str, record_id, desc)
            self.get_logger().debug(msg)
            continue
          except TypeError:
            msg = "Expression value is not of valid type (float) at record: {} {}".format(record_id, key_values[
              self._expression_str])
            self.get_logger().debug(msg)
            continue

        # translate the whole sequences (3 ORFs) for non CDS sequences and not take alt_ORFs for CDSs
        if 'CDS' not in key_values.keys() or ('CDS' in key_values.keys() and
                                              ('altORFs' in self._include_biotypes or
                                               self._include_biotypes == ['all'])):
          ref_orfs = self.get_orfs_dna(ref_seq, self._translation_table, self._num_orfs,
                                       self._num_orfs_complement, to_stop=False)
          print(self._header_var_prefix)
          self.write_output(seq_id=self._header_var_prefix + record_id, desc=desc, seqs=ref_orfs, prots_fn=prots_fn)

        # also allow for direct translation of the CDS, when the cds info exists in the fasta header skip_including_all_cds is false
        if 'CDS' in key_values.keys() and not self._skip_including_all_cds:
          try:
            cds_info = [int(x) for x in key_values['CDS'].split('-')]
            ref_seq = ref_seq[cds_info[0] - 1:cds_info[1]]
            ref_orfs = self.get_orfs_dna(ref_seq, self._translation_table, 1, 0, to_stop=True)
            self.write_output(seq_id=record_id, desc=desc, seqs=ref_orfs, prots_fn=prots_fn)
          except (ValueError, IndexError, KeyError):
            print("Could not extra cds position from fasta header for: ", record_id, desc)

    return self._proteindb_output

  @staticmethod
  def get_key(fasta_header):
    return fasta_header.split('|')[0].split(' ')[0]

  def vcf_to_proteindb(self, vcf_file, input_fasta, gene_annotations_gtf):
    """
        Generate peps for variants by modifying sequences of affected transcripts (VCF - VEP annotated).
        It only considers variants within potential coding regions of the transcript
        (CDSs & stop codons for protein-coding genes, exons for non-protein coding genes).
        :param vcf_file:
        :param input_fasta:
        :param gene_annotations_gtf:
        :return:
        """

    db = self.parse_gtf(gene_annotations_gtf, gene_annotations_gtf.replace('.gtf', '.db'))

    transcripts_dict = SeqIO.index(input_fasta, "fasta", key_function=self.get_key)
    # handle cases where the transcript has version in the GTF but not in the VCF
    transcript_id_mapping = {k.split('.')[0]: k for k in transcripts_dict.keys()}
    with open(self._proteindb_output, 'w') as prots_fn:
      vcf_reader = vcf.Reader(open(vcf_file, 'r'))

      for record in vcf_reader:
        if record.ALT == [None] or record.REF == [None]:
          msg = "Invalid VCF record, skipping: {}".format(record)
          self.get_logger().debug(msg)
          continue
        if not self._ignore_filters:
          if record.FILTER:  # if not PASS: None and empty means PASS
            if not (set(record.FILTER[0].split(',')) <= set(self._accepted_filters)):
              continue

        # only process variants above a given allele frequency threshold if the AF string is not empty
        if self._af_field:
          # get AF from the INFO field
          try:
            af = float(record.INFO[self._af_field])
          except TypeError:
            af = float(record.INFO[self._af_field][0])
          except KeyError:
            continue

          # check if the AF passed the threshold
          if af < self._af_threshold:
            continue

        trans_table = self._translation_table
        consequences = []
        if str(record.CHROM).lstrip('chr').upper() in ['M', 'MT']:
          trans_table = self._mito_translation_table

        processed_transcript_allele = []

        for transcript_record in record.INFO[self._annotation_field_name]:
          transcript_info = transcript_record.split('|')
          try:
            consequence = transcript_info[self._consequence_index]
          except IndexError:
            msg = "Give a valid index for the consequence in the INFO field for: {}".format(transcript_record)
            self.get_logger().debug(msg)
            continue
          consequences.append(consequence)
          try:
            transcript_id = transcript_info[self._transcript_index]
          except IndexError:
            msg = "Give a valid index for the Transcript ID in the INFO field for: {}".format(transcript_record)
            self.get_logger().debug(msg)
            continue
          if transcript_id == "":
            continue

          try:
            transcript_id_v = transcript_id_mapping[transcript_id]
          except KeyError:
            transcript_id_v = transcript_id

          try:
            row = transcripts_dict[transcript_id_v]
            ref_seq = row.seq  # get the seq and desc for the transcript from the fasta of the gtf
            desc = str(row.description)
          except KeyError:
            msg = "Transcript {} not found in fasta of the GTF file {}".format(transcript_id_v, record)
            self.get_logger().debug(msg)
            continue

          feature_types = ['exon']
          # check if cds info exists in the fasta header otherwise translate all exons
          cds_info = []
          num_orfs = 3
          if 'CDS=' in desc:
            try:
              cds_info = [int(x) for x in desc.split(' ')[1].split('=')[1].split('-')]
              feature_types = ['CDS', 'stop_codon']
              num_orfs = 1
            except (ValueError, IndexError):
              msg = "Could not extra cds position from fasta header for: {}".format(desc)
              self.get_logger().debug(msg)

          chrom, strand, features_info, feature_biotype = self.get_features(db, transcript_id_v,
                                                                            self._biotype_str,
                                                                            feature_types)
          if chrom is None:  # the record info was not found
            continue
          # skip transcripts with unwanted consequences
          if (consequence in self._exclude_consequences or
            (consequence not in self._include_consequences and
             self._include_consequences != ['all'])):
            continue

          # only include features that have the specified biotypes or they have CDSs info
          if 'CDS' in feature_types and not self._skip_including_all_cds:
            pass
          elif (feature_biotype in self._exclude_biotypes or
                (feature_biotype not in self._include_biotypes and
                 self._include_biotypes != ['all'])):
            continue
          for alt in record.ALT:  # in cases of multiple alternative alleles consider all
            if alt is None:
              continue
            if transcript_id + str(record.REF) + str(
              alt) in processed_transcript_allele:  # because VEP reports affected transcripts per alt allele
              continue

            processed_transcript_allele.append(transcript_id + str(record.REF) + str(alt))
            # for non-CDSs, only consider the exon that actually overlaps the variant

            try:
              overlap_flag = self.check_overlap(record.POS, record.POS + len(alt), features_info)
            except TypeError:
              msg = "Wrong VCF record in {}".format(record)
              self.get_logger().debug(msg)
              continue

            if (chrom.lstrip("chr") == str(record.CHROM).lstrip("chr") and
              overlap_flag):
              coding_ref_seq, coding_alt_seq = self.get_altseq(ref_seq, Seq(str(record.REF)),
                                                               Seq(str(alt)), int(record.POS), strand,
                                                               features_info, cds_info)
              if coding_alt_seq != "":
                ref_orfs, alt_orfs = self.get_orfs_vcf(coding_ref_seq, coding_alt_seq, trans_table,
                                                       num_orfs)
                record_id = ""
                if record.ID:
                  record_id = '_' + str(record.ID)
                self.write_output(seq_id='_'.join([self._header_var_prefix + str(record_id),
                                                   '.'.join([str(record.CHROM), str(record.POS),
                                                             str(record.REF), str(alt)]),
                                                   transcript_id_v]),
                                  desc=feature_biotype + ":" + consequence,
                                  seqs=alt_orfs,
                                  prots_fn=prots_fn)

                if self._report_reference_seq:
                  self.write_output(seq_id=transcript_id_v,
                                    desc=feature_biotype,
                                    seqs=ref_orfs,
                                    prots_fn=prots_fn)

    return self._proteindb_output

  @staticmethod
  def add_protein_to_map(seq: str, new_desc_string: str, protein_id: str, proteins, output_handle):
    protein = {'description': new_desc_string, 'sequence': seq, 'accession': protein_id}
    proteins.append(protein)
    output_handle.write(">{}\t{}\n{}\n".format(protein_id, new_desc_string, seq))
    return proteins

  def check_proteindb(self, input_fasta: str = None, add_stop_codon: bool = False, num_aa: int = 6):

    input_handle = open(input_fasta, 'r')
    output_handle = open(self._proteindb_output, 'w')
    proteins = []
    pcount = 0
    stop_count = 0
    gap_count = 0
    no_met = 0
    less = 0

    for record in SeqIO.parse(input_handle, 'fasta'):

      seq = str(record.seq)
      pcount += 1

      # parse the description string into a dictionary
      new_desc_string = record.description
      new_desc_string = new_desc_string[new_desc_string.find(' ') + 1:]
      # test for odd amino acids, stop codons, gaps
      if not seq.startswith('M'):
        no_met += 1
      if seq.endswith('*'):
        seq = seq[:-1]
      if '-' in seq:
        gap_count += 1
        new_desc_string = new_desc_string + ' (Contains gaps)'
      if '*' in seq:
        stop_count += 1
        if add_stop_codon:
          seq_list = seq.split("*")
          codon_index = 1
          for codon in seq_list:
            new_desc_string = new_desc_string + ' codon ' + str(codon_index)
            protein_id = record.id + '_codon_' + str(codon_index)
            seq = codon
            if len(seq) > num_aa:
              proteins = self.add_protein_to_map(seq, new_desc_string, protein_id, proteins, output_handle)
            codon_index = codon_index + 1
        else:
          cut = seq.index('*')
          string = ' (Premature stop %s/%s)' % (cut, len(seq))
          new_desc_string = new_desc_string + string
          seq = seq[:cut]
          # save the protein in list
          if len(seq) > num_aa:
            protein_id = record.id
            proteins = self.add_protein_to_map(seq, new_desc_string, protein_id, proteins, output_handle)
          else:
            less += 1
      else:
        if len(seq) > num_aa:
          protein_id = record.id
          proteins = self.add_protein_to_map(seq, new_desc_string, protein_id, proteins, output_handle)
        else:
          less += 1

    print("   translations that do not start with Met:", no_met)
    print("   translations that have premature stop codons:", stop_count)
    print("   translations that contain gaps:", gap_count)
    print("   total number of input sequences was:", pcount)
    print("   total number of sequences written was:", len(proteins))
    print("   total number of proteins less than {} aminoacids: {}".format(num_aa, less))

  @staticmethod
  def write_output(seq_id, desc, seqs, prots_fn):
    """
    write the orfs to the output file
    :param seq_id: Sequence Accession
    :param desc: Sequence Description
    :param seqs: Sequence
    :param prots_fn:
    :return:
    """
    write_i = False
    if len(seqs) > 1:  # only add _num when multiple ORFs are generated (e.g in 3 ORF)
      write_i = True

    for i, orf in enumerate(seqs):
      if write_i:  # only add _num when multiple ORFs are generated (e.g in 3 ORF)
        prots_fn.write('>{} {}\n{}\n'.format(seq_id + "_" + str(i + 1), desc, orf))
      else:
        prots_fn.write('>{} {}\n{}\n'.format(seq_id, desc, orf))


if __name__ == '__main__':
  print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
