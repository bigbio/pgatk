#
# DecoyPYrat - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses
#
# MIT License
#
# Copyright (c) 2016 James Christopher Wright - Wellcome Trust Sanger Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import random
import os
from Bio import SeqIO

from pypgatk.toolbox.general import ParameterConfiguration


class ProteinDBService(ParameterConfiguration):
  CONFIG_KEY_PROTEINDB_DECOY = 'proteindb_decoy'
  CONFIG_PROTEINDB_OUTPUT = 'output'
  CONFIG_INPUT_FILE = 'input'
  CONFIG_CLEAVAGE_SITES = 'cleavage_sites'
  CONFIG_CLEAVAGE_POSITION = 'cleavage_position'
  CONFIG_ANTI_CLEAVAGE_SITES = 'anti_cleavage_sites'
  CONFIG_PEPTIDE_LENGTH = 'min_peptide_length'
  CONFIG_MAX_ITERATIONS = 'max_iterations'
  CONFIG_DO_NOT_SUFFLE = 'do_not_shuffle'
  CONFIG_KEEP_TARGET_HITS = 'keep_target_hits'
  CONFIG_DO_NOT_SWITCH = 'do_not_switch'
  CONFIG_DECOY_PREFIX = 'decoy_prefix'
  CONFIG_TEMP_FILE = 'temp_file'
  CONFIG_NO_ISOBARIC = 'no_isobaric'
  
  def __init__(self, config_file, pipeline_arguments):
    """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """

    super(ProteinDBService, self).__init__(self.CONFIG_KEY_PROTEINDB_DECOY, config_file,
                                           pipeline_arguments)

    self._temp_file = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_TEMP_FILE]
    if self.CONFIG_TEMP_FILE in self.get_pipeline_parameters():
      self._temp_file = self.get_pipeline_parameters()[self.CONFIG_TEMP_FILE]

    self._input_fasta = self.get_pipeline_parameters()[self.CONFIG_INPUT_FILE]

    self._isobaric = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_NO_ISOBARIC]
    if self.CONFIG_NO_ISOBARIC in self.get_pipeline_parameters():
      self._isobaric = self.get_pipeline_parameters()[self.CONFIG_NO_ISOBARIC]

    self._cleavage_sites = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_CLEAVAGE_SITES]
    if self.CONFIG_CLEAVAGE_SITES in self.get_pipeline_parameters():
      self._cleavage_sites = self.get_pipeline_parameters()[self.CONFIG_CLEAVAGE_SITES]

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

    self._peptide_length = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_PEPTIDE_LENGTH]
    if self.CONFIG_PEPTIDE_LENGTH in self.get_pipeline_parameters():
      self._peptide_length = self.get_pipeline_parameters()[self.CONFIG_PEPTIDE_LENGTH]

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

    self._keep_target_hits = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][
      self.CONFIG_KEEP_TARGET_HITS]
    if self.CONFIG_KEEP_TARGET_HITS in self.get_pipeline_parameters():
      self._keep_target_hits = self.get_pipeline_parameters()[self.CONFIG_KEEP_TARGET_HITS]
    
    
  @staticmethod
  def digest(protein, sites, pos, no, min_peptide_length):
    """
        Return a list of cleaved peptides with minimum length in protein sequence.
        :param protein = sequence
        :param sites = string of amino acid cleavage sites
        :param pos = n or c for n-terminal or c-terminal cleavage
        :param no = amino acids following site that would prevent cleavage ie proline
        :param min_peptide_length = minimum length of peptides returned
        """

    # for each possible cleavage site insert a comma with before or after depending on pos
    for s in sites:
      r = s + ','
      if pos == 'n':
        r = ',' + s
      protein = protein.replace(s, r)

    # for each possible cleavage and all none cleavage remove comma
    for s in sites:
      for n in no:
        a = s + ',' + n
        if pos == 'n':
          a = ',' + s + n
        r = s + n
        protein = protein.replace(a, r)

    # filter peptides into list by minimum size
    return list(filter(lambda x: len(x) >= min_peptide_length, (protein.split(','))))

  @staticmethod
  def revswitch(protein, noswitch, sites):
    """
        Return a reversed protein sequence with cleavage residues switched with preceding residue
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
        shuffle peptide without moving c-terminal amino acid cleavage site
        """

    # extract terminal aa
    s = peptide[-1]
    # convert peptide to list (remove K/R) and shuffle the list
    l = list(peptide[:-1])
    random.shuffle(l)
    # return new peptide
    return ''.join(l) + s

  def decoy_database_full_protein(self):
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
        upeps.update(ProteinDBService.digest(seq, self._cleavage_sites, self._cleavage_position,
                                             self._anti_cleavage_sites, self._peptide_length))

    # open orary decoy FASTA file
    with open(self._output_file, 'w') as outfa:
      fasta = SeqIO.parse(self._input_fasta, 'fasta')
      targets = []
      decoys = []
      for i,record in enumerate(fasta):
        protseq = str(record.seq)
        targets.append(protseq)
        revprotseq = []
        for seq in protseq.split('*'):
          if not seq:
              continue
          if not self._isobaric:
            seq = seq.replace('I', 'L')

          # reverse and switch protein sequence
          decoyseq = ProteinDBService.revswitch(seq, self._no_switch, self._cleavage_sites)

          decoy_peps = ProteinDBService.digest(decoyseq, self._cleavage_sites, self._cleavage_position,
                                               self._anti_cleavage_sites, 0)
          #if any of the digested peptides are found in the targets (upeps) then shuffle
          checked_decoy_peps = []
          for decoy_pep in decoy_peps:
            if len(decoy_pep) < self._peptide_length:
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
                aPep = ProteinDBService.shuffle(aPep)

                # check if shuffling has an effect if not end iterations
                if aPep == decoy_pep:
                  i = self._max_iterations

                # warn if peptide has no suitable alternative, add to removal list
                if i == self._max_iterations:
                  noAlternative.add(decoy_pep)
                  aPep = ''
            #if decoy is generated then add to the list of peptides
            if aPep:
              checked_decoy_peps.append(aPep)
            else:
              if self._keep_target_hits:
                checked_decoy_peps.append(decoy_pep)
          #finally join the peptides to generate protein decoy
          if checked_decoy_peps:
            revprotseq.append(''.join(checked_decoy_peps))

        outfa.write('>{}\n{}\n'.format(self._decoy_prefix+str(record.id), '*'.join(revprotseq)))
        decoys.append('*'.join(revprotseq))

      with open(self._output_file.replace('.fa','')+'_noAlternative.fa', 'w') as noAlternative_outfa:
        noAlternative_outfa.write('\n'.join(noAlternative) + '\n')
      print('Number of skipped tryptic peptides in decoy db (no alternatives): {}'.
            format(len(noAlternative)))
      print('Total number of amino acids in target and decoy databases: ', 
            len(''.join(targets)), len(''.join(decoys)))

