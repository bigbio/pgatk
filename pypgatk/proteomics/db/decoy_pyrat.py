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
    CONFIG_DO_NOT_SWITCH = 'do_not_switch'
    CONFIG_DECOY_PREFIX = 'decoy_prefix'
    CONFIG_TEMP_FILE = 'temp_file'
    CONFIG_NO_ISOBARIC = 'no_isobaric'
    CONFIG_MEMORY_SAVE = 'memory_save'

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

        self._memory_save = self.get_default_parameters()[self.CONFIG_KEY_PROTEINDB_DECOY][self.CONFIG_MEMORY_SAVE]
        if self.CONFIG_MEMORY_SAVE in self.get_pipeline_parameters():
            self._memory_save = self.get_pipeline_parameters()[self.CONFIG_MEMORY_SAVE]

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

    def decoy_database(self):
        """
        Create a decoy database from a proteomics database
        :return:
        """

        # Create empty sets to add all target and decoy peptides
        upeps = set()
        dpeps = set()

        # Counter for number of decoy sequences
        dcount = 0

        # empty protein sequence
        seq = ''

        # open temporary decoy FASTA file
        outfa = open(self._temp_file, 'w')

        # Open FASTA file using first cmd line argument
        fasta = SeqIO.parse(self._input_fasta, 'fasta')
        # loop each seq in the file
        for record in fasta:
                seq = str(record.seq)
                dcount+=1
                #make sequence isobaric (check args for switch off)
                if not self._isobaric:
                        seq = seq.replace('I', 'L')

                #digest sequence add peptides to set
                upeps.update(ProteinDBService.digest(seq, self._cleavage_sites, self._cleavage_position,
                                                         self._anti_cleavage_sites, self._peptide_length))

                #reverse and switch protein sequence
                decoyseq = ProteinDBService.revswitch(seq, self._no_switch, self._cleavage_sites)

                #do not store decoy peptide set in reduced memory mode
                if not self._memory_save:
                        #update decoy peptide set
                        dpeps.update( ProteinDBService.digest(decoyseq, self._cleavage_sites, self._cleavage_position,
                                                             self._anti_cleavage_sites, self._peptide_length))

                #write decoy protein accession and sequence to file
                outfa.write('>' + self._decoy_prefix + record.id + '\n')
                outfa.write(decoyseq + '\n')

        # Close files
        fasta.close()
        outfa.close()

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
                        for p in ProteinDBService.digest(line.rstrip(), self._cleavage_sites, self._cleavage_position,
                                                         self._anti_cleavage_sites, self._peptide_length):
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
                    aPep = ProteinDBService.shuffle(dPep)

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
                # open original decoy file
                with open(self._temp_file, "rt") as fin:
                    # loop each line of original decoy fasta
                    for line in fin:
                        # if line is not accession replace peptides in dictionary with alternatives
                        if line[0] != '>':
                            # digest decoy sequence
                            for p in ProteinDBService.digest(line.rstrip(), self._cleavage_sites,
                                                             self._cleavage_position, self._anti_cleavage_sites,
                                                             self._peptide_length):
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
