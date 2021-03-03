import random
import os
from Bio import SeqIO
from pyteomics.fasta import decoy_sequence
from pyteomics.parser import cleave

from pypgatk.proteomics.models import PYGPATK_ENZYMES
from pypgatk.toolbox.exceptions import AppException
from pypgatk.toolbox.general import ParameterConfiguration


class ProteinDBDecoyService(ParameterConfiguration):
  CONFIG_KEY_PROTEINDB_DECOY = 'proteindb_decoy'
  CONFIG_PROTEINDB_OUTPUT = 'output'
  CONFIG_INPUT_FILE = 'input'
  CONFIG_DECOY_METHOD='method'
  CONFIG_ENZYME = 'enzyme'
  CONFIG_CLEAVAGE_POSITION = 'cleavage_position'
  CONFIG_MAX_MISSED_CLEAVAGES = 'max_missed_cleavages'
  CONFIG_ANTI_CLEAVAGE_SITES = 'anti_cleavage_sites'
  CONFIG_MIN_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_PEPTIDE_LENGTH_MAX = 'max_peptide_length'
  CONFIG_MAX_ITERATIONS = 'max_iterations'
  CONFIG_DO_NOT_SUFFLE = 'do_not_shuffle'
  CONFIG_DO_NOT_SWITCH = 'do_not_switch'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_TEMP_FILE = 'temp_file'
  CONFIG_NO_ISOBARIC = 'no_isobaric'
  CONFIG_MEMORY_SAVE = 'memory_save'
  CONFIG_KEEP_TARGET_HITS = 'keep_target_hits'

  def __init__(self, config_file, pipeline_arguments):
    """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """

    super(ProteinDBDecoyService, self).__init__(self.CONFIG_KEY_PROTEINDB_DECOY, config_file,
                                                pipeline_arguments)

    self._temp_file = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_TEMP_FILE]
    if self.CONFIG_TEMP_FILE in self.get_pipeline_parameters():
      self._temp_file = self.get_pipeline_parameters()[self.CONFIG_TEMP_FILE]

    self._input_fasta = self.get_pipeline_parameters()[self.CONFIG_INPUT_FILE]

    self._isobaric = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_NO_ISOBARIC]
    if self.CONFIG_NO_ISOBARIC in self.get_pipeline_parameters():
      self._isobaric = self.get_pipeline_parameters()[self.CONFIG_NO_ISOBARIC]

    self._memory_save = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_MEMORY_SAVE]
    if self.CONFIG_MEMORY_SAVE in self.get_pipeline_parameters():
      self._memory_save = self.get_pipeline_parameters()[self.CONFIG_MEMORY_SAVE]

    self._enzyme = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_ENZYME]
    if self.CONFIG_ENZYME in self.get_pipeline_parameters():
      self._enzyme = self.get_pipeline_parameters()[self.CONFIG_ENZYME]

    self._decoy_prefix = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_DECOY_PREFIX]
    if self.CONFIG_DECOY_PREFIX in self.get_pipeline_parameters():
      self._decoy_prefix = self.get_pipeline_parameters()[self.CONFIG_DECOY_PREFIX]

    self._cleavage_position = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_CLEAVAGE_POSITION]
    if self.CONFIG_CLEAVAGE_POSITION in self.get_pipeline_parameters():
      self._cleavage_position = self.get_pipeline_parameters()[self.CONFIG_CLEAVAGE_POSITION]

    self._anti_cleavage_sites = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_ANTI_CLEAVAGE_SITES]
    if self.CONFIG_ANTI_CLEAVAGE_SITES in self.get_pipeline_parameters():
      self._anti_cleavage_sites = self.get_pipeline_parameters()[self.CONFIG_ANTI_CLEAVAGE_SITES]

    self._min_peptide_length = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_MIN_PEPTIDE_LENGTH]
    if self.CONFIG_MIN_PEPTIDE_LENGTH in self.get_pipeline_parameters():
      self._min_peptide_length = self.get_pipeline_parameters()[self.CONFIG_MIN_PEPTIDE_LENGTH]

    self._max_iterations = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_MAX_ITERATIONS]
    if self.CONFIG_MAX_ITERATIONS in self.get_pipeline_parameters():
      self._max_iterations = self.get_pipeline_parameters()[self.CONFIG_MAX_ITERATIONS]

    self._no_switch = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_DO_NOT_SWITCH]
    if self.CONFIG_DO_NOT_SWITCH in self.get_pipeline_parameters():
      self._no_switch = self.get_pipeline_parameters()[self.CONFIG_DO_NOT_SWITCH]

    self._output_file = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_PROTEINDB_OUTPUT]
    if self.CONFIG_PROTEINDB_OUTPUT in self.get_pipeline_parameters():
      self._output_file = self.get_pipeline_parameters()[self.CONFIG_PROTEINDB_OUTPUT]

    self._no_suffle = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_DO_NOT_SUFFLE]
    if self.CONFIG_DO_NOT_SUFFLE in self.get_pipeline_parameters():
      self._no_suffle = self.get_pipeline_parameters()[self.CONFIG_DO_NOT_SUFFLE]

    self._method = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_DECOY_METHOD]
    if self.CONFIG_DECOY_METHOD in self.get_pipeline_parameters():
      self._method = self.get_pipeline_parameters()[self.CONFIG_DECOY_METHOD]

    self._max_peptide_length  = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_PEPTIDE_LENGTH_MAX]
    if self.CONFIG_PEPTIDE_LENGTH_MAX in self.get_pipeline_parameters():
      self._max_peptide_length = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_LENGTH_MAX]

    self._max_missed_cleavages = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_MAX_MISSED_CLEAVAGES]
    if self.CONFIG_MAX_MISSED_CLEAVAGES in self.get_pipeline_parameters():
      self._max_missed_cleavages = self.get_pipeline_parameters()[self.CONFIG_MAX_MISSED_CLEAVAGES]

    self._max_missed_cleavages = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_KEEP_TARGET_HITS]
    if self.CONFIG_KEEP_TARGET_HITS in self.get_pipeline_parameters():
      self._keep_target_hits = self.get_pipeline_parameters()[self.CONFIG_KEEP_TARGET_HITS]


  @staticmethod
  def revswitch(protein, noswitch, sites):
    """
    Return a reversed protein sequence with cleavage residues switched with preceding residue. This method is used by
    DecoyPyrat.
    :return protein
    """
    # reverse protein sequence with a reverse splice convert to list
    revseq = list(protein[::-1])

    if not noswitch:
      # loop sequence list
      for i, c in enumerate(revseq):
        # if value is cleavage site switch with previous amino acid
        for s in sites:
          if c == s:
            aa = revseq[i - 1]
            revseq[i - 1] = revseq[i]
            revseq[i] = aa

    # return reversed with/without switched proteins as string
    return ''.join(revseq)

  @staticmethod
  def shuffle(peptide):
    """
    shuffle peptide without moving c-terminal amino acid cleavage site.
    :param peptide sequence
    :return new shuffle peptide.
    """

    # extract terminal aa
    s = peptide[-1]
    # convert peptide to list (remove K/R) and shuffle the list
    l = list(peptide[:-1])
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s

  def protein_database_decoy(self, method ='reverse'):
    """
    Reverse protein sequences and attach them to the output fasta file
    :param method the method used to compute the decoy
    :return:
    """
    fasta = SeqIO.parse(self._input_fasta, 'fasta')
    with open(self._output_file, "wt") as output_file:
      for record in fasta:
        seq = str(record.seq)
        decoy_seq = decoy_sequence(seq, mode=method)

        output_file.write('>' + record.id + "\t" + record.description + '\n')
        output_file.write(str(record.seq) + '\n')
        output_file.write('>' + self._decoy_prefix + record.id + '\n')
        output_file.write(decoy_seq + '\n')

    output_file.close()

  def print_target_decoy_composition(self, min_length: int = 0, max_length: int = 100):
    """
    Print the number of target peptides vs decoy peptides in a Fasta database
    :param min_length the minimum length of a peptide
    :param max_length the maximum length of a peptide
    :return:
    """

    fasta = SeqIO.parse(self._output_file, 'fasta')
    target_peptides = {}
    decoy_peptides = {}
    pep_count_in_both = 0
    for record in fasta:
      peptides = cleave(sequence = str(record.seq), rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=self._max_missed_cleavages, min_length = self._min_peptide_length)
      for peptide in peptides:
         if self._decoy_prefix in record.id:
           decoy_peptides[peptide] = 'decoy'
           if peptide in target_peptides:
             target_peptides.pop(peptide)
             pep_count_in_both += 1
         else:
           if peptide not in decoy_peptides:
             target_peptides[peptide] = 'target'

    print('Number of target peptides: {} and Decoy Peptides: {}'.format(len(target_peptides), len(decoy_peptides)))
    target_percentage = (len(target_peptides)/(len(target_peptides)+len(decoy_peptides)))*100
    print('% Target peptides {:.1f}'.format(target_percentage))
    decoy_percentage = (len(decoy_peptides) / (len(target_peptides) + len(decoy_peptides))) * 100
    print('% Decoy peptides {:.1f}'.format(decoy_percentage))
    duplicate_percentage = (pep_count_in_both/(len(target_peptides) + len(decoy_peptides))) * 100
    print('Number of peptides in Target and Decoy {}, Percentage {:.1f}'.format(pep_count_in_both, duplicate_percentage))

  def generate_decoypyrat_database(self):
    """
    Create a decoy database from a proteomics database this method is presented in manuscript:
    J Proteomics Bioinform. 2016 Jun 27; 9(6): 176–180. PMCID: PMC4941923
    DecoyPyrat: Fast Non-redundant Hybrid Decoy Sequence Generation for Large Scale Proteomics
    :return:
    """

    # Create empty sets to add all target and decoy peptides
    upeps = set()
    dpeps = set()

    # Counter for number of decoy sequences
    dcount = 0

    # Open FASTA file using first cmd line argument
    fasta = SeqIO.parse(self._input_fasta, 'fasta')

    # open temporary decoy FASTA file
    with open(self._temp_file, 'w') as outfa:

      # loop each seq in the file
      for record in fasta:
        seq = str(record.seq)
        dcount += 1
        # make sequence isobaric (check args for switch off)
        if not self._isobaric:
          seq = seq.replace('I', 'L')

        # digest sequence add peptides to set
        upeps.update(cleave(sequence = seq, rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=0, min_length=self._min_peptide_length))

        # reverse and switch protein sequence
        decoyseq = self.revswitch(seq, self._no_switch, PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage sites'])

        # do not store decoy peptide set in reduced memory mode
        if not self._memory_save:
          # update decoy peptide set
          dpeps.update(cleave(sequence = decoyseq, rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=0, min_length=self._min_peptide_length))

        # write decoy protein accession and sequence to file
        outfa.write('>' + self._decoy_prefix + record.id + '\n')
        outfa.write(decoyseq + '\n')

    # Summarise the numbers of target and decoy peptides and their intersection
    nonDecoys = set()
    print("proteins:" + str(dcount))
    print("target peptides:" + str(len(upeps)))

    # Reloop decoy file in reduced memory mode to store only intersecting decoys
    if self._memory_save:
      # open temp decoys
      with open(self._temp_file, "rt") as fin:
        for line in fin:
          # if line is not accession
          if line[0] != '>':
            # digest protein
            for p in cleave(sequence = line.rstrip(), rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=0, min_length=self._min_peptide_length):
              # check if in target peptides if true then add to nonDecoys
              if p in upeps:
                nonDecoys.add(p)
      fin.close()
      print("decoy peptides: !Memory Saving Made!")
    else:
      # can only report total number in normal memory mode
      print("decoy peptides:" + str(len(dpeps)))
      # find intersecting peptides
      nonDecoys = upeps.intersection(dpeps)

    print("#intersection:" + str(len(nonDecoys)))

    # if there are decoy peptides that are in the target peptide set
    if len(nonDecoys) > 0 and self._no_suffle == False:

      # create empty dictionary with bad decoys as keys
      dAlternative = dict.fromkeys(nonDecoys, '')
      noAlternative = list()

      # loop bad decoys / dictionary keys
      for dPep in dAlternative:
        i = 0
        aPep = dPep

        # shuffle until aPep is not in target set (maximum of 10 iterations)
        while aPep in upeps and i < self._max_iterations:

          # increment iteration counter
          i += 1

          # shuffle peptide
          aPep = self.shuffle(dPep)

          # check if shuffling has an effect if not end iterations
          if aPep == dPep:
            i = self._max_iterations

        # update dictionary with alternative shuffled peptide
        dAlternative[dPep] = aPep

        # warn if peptide has no suitable alternative, add to removal list
        if i == self._max_iterations:
          noAlternative.append(dPep)

      print(str(len(noAlternative)) + ' have no alternative peptide')
      # remove peptides with no alternative
      for p in noAlternative:
        del dAlternative[p]

      # Free up memory by clearing large sets of peptides
      upeps.clear()
      dpeps.clear()

      # open second decoy file
      with open(self._output_file, "wt") as fout:

        # Attach the target sequences to the database
        fasta = SeqIO.parse(self._input_fasta, 'fasta')
        for record in fasta:
          seq = str(record.seq)
          id  = record.id
          description = record.description
          fout.write('>'+ id + '\t' + description + '\n')
          fout.write(seq + '\n')

        # open original decoy file
        with open(self._temp_file, "rt") as fin:
          # loop each line of original decoy fasta
          for line in fin:
            # if line is not accession replace peptides in dictionary with alternatives
            if line[0] != '>':
              # digest decoy sequence
              for p in cleave(sequence = line.rstrip(), rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=0, min_length=self._min_peptide_length):
                # store decoy peptide for final count
                dpeps.add(p)

                # if decoy peptide is in dictionary replace with alternative
                if p in dAlternative:
                  line = line.replace(p, dAlternative[p])

            fout.write(line)
        fin.close()
      fout.close()

      # delete temporary file
      os.remove(self._temp_file)
    else:
      os.rename(self._temp_file, self._output_file)

    print("final decoy peptides:" + str(len(dpeps)))

  def pypgatk_decoy_database(self):
    """
    Create a decoy database from a proteomics database
    target db is digested and only digested peptides > _peptide_length are kept
    next, each target protein is reversed and digested, all peptides are kept
    regardless of their length.
    The list of digested peptides from the reversed protein are iterated:
     - small peptides are kept (len < _peptide_length)
     - peptides not found in target are kep
     - peptides with a match are shuffled for max_iterations, if a non-target
       peptide was found then written otherwise the peptide is skipped unless
      the _keep_target_hits option is true.
     :return:
    """

    # Create empty sets to add all target and decoy peptides
    upeps = set()
    noAlternative = set()
    # Open FASTA file using first cmd line argument
    fasta = SeqIO.parse(self._input_fasta, 'fasta')
    # loop each seq in the file
    for record in fasta:
      seq = str(record.seq)
      if not self._isobaric:
        seq = seq.replace('I', 'L')

        # digest sequence add peptides to set
        upeps.update(cleave(sequence = seq, rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=self._max_missed_cleavages, min_length=self._min_peptide_length))

    # open orary decoy FASTA file
    with open(self._output_file, 'w') as outfa:
      fasta = SeqIO.parse(self._input_fasta, 'fasta')
      targets = []
      decoys = []
      for i, record in enumerate(fasta):
        protseq = str(record.seq)
        targets.append(protseq)
        revprotseq = []

        # output target protein
        seq = str(record.seq)
        id = record.id
        description = record.description
        outfa.write('>' + id + '\t' + description + '\n')
        outfa.write(seq + '\n')

        for seq in protseq.split('*'):
          if not seq:
            continue
          if not self._isobaric:
            seq = seq.replace('I', 'L')

          # reverse and switch protein sequence
          decoyseq = self.revswitch(seq, self._no_switch, PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage sites'])

          decoy_peps = cleave(sequence = decoyseq, rule = PYGPATK_ENZYMES.enzymes[self._enzyme]['cleavage rule'], missed_cleavages=0, min_length=0)

          # if any of the digested peptides are found in the targets (upeps) then shuffle
          checked_decoy_peps = []
          for decoy_pep in decoy_peps:
            if len(decoy_pep) < self._min_peptide_length:
              checked_decoy_peps.append(decoy_pep)
              continue

            found_in_target = False
            aPep = ''
            if decoy_pep in upeps:
              found_in_target = True
            else:
              checked_decoy_peps.append(decoy_pep)
              continue

            if found_in_target and not self._no_suffle and decoy_pep not in noAlternative:
              aPep = decoy_pep
              # shuffle until aPep is not in target set (maximum of 10 iterations)
              i = 0
              while aPep in upeps and i < self._max_iterations:
                # increment iteration counter
                i += 1
                # shuffle peptide
                aPep = self.shuffle(aPep)

                # check if shuffling has an effect if not end iterations
                if aPep == decoy_pep:
                  i = self._max_iterations

                # warn if peptide has no suitable alternative, add to removal list
                if i == self._max_iterations:
                  noAlternative.add(decoy_pep)
                  aPep = ''
            # if decoy is generated then add to the list of peptides
            if aPep:
              checked_decoy_peps.append(aPep)
            else:
              if self._keep_target_hits:
                checked_decoy_peps.append(decoy_pep)
          # finally join the peptides to generate protein decoy
          if checked_decoy_peps:
            revprotseq.append(''.join(checked_decoy_peps))

        outfa.write('>{}\n{}\n'.format(self._decoy_prefix + str(record.id), '*'.join(revprotseq)))
        decoys.append('*'.join(revprotseq))

      with open(self._output_file.replace('.fa', '') + '_noAlternative.fa', 'w') as noAlternative_outfa:
        noAlternative_outfa.write('\n'.join(noAlternative) + '\n')
      print('Number of skipped tryptic peptides in decoy db (no alternatives): {}'.
            format(len(noAlternative)))
      print('Total number of amino acids in target and decoy databases: ',
            len(''.join(targets)), len(''.join(decoys)))



  def decoy_database(self):
    """
    This method is used to pick the rigth decoy method to generate the decoys.
    :return:
    """
    if self._method == 'protein-reverse':
      self.protein_database_decoy(method='reverse')
    elif self._method == 'protein-shuffle':
      self.protein_database_decoy(method='shuffle')
    elif self._method == 'decoypyrat':
      self.generate_decoypyrat_database()
    elif self._method == 'pgdbdeep':
      self.pypgatk_decoy_database()
    else:
      raise AppException("The following method {} is not supported by the library".format(self._method))

