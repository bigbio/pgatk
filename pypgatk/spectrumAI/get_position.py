import pandas as pd
import re
from Bio import SeqIO

from pypgatk.toolbox.general import ParameterConfiguration

class GetPosition(ParameterConfiguration):
    CONFIG_KEY_GET_POSITION = 'get_position'
    CONFIG_INPUT_PSM_TABLE = 'input_psm_table'
    CONFIG_INPUT_FASTA = 'input_fasta'
    CONFIG_OUTPUT_PSM_TABLE = 'output_psm_table'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(GetPosition, self).__init__(self.CONFIG_KEY_GET_POSITION, config_data, pipeline_arguments)

        if self.CONFIG_INPUT_PSM_TABLE in self.get_pipeline_parameters():
            self._input_psm_table = self.get_pipeline_parameters()[self.CONFIG_INPUT_PSM_TABLE]
        elif self.CONFIG_KEY_GET_POSITION in self.get_default_parameters() and \
                self.CONFIG_INPUT_PSM_TABLE in self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION]:
            self._input_psm_table = self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION][
                self.CONFIG_INPUT_PSM_TABLE]
        
        if self.CONFIG_INPUT_FASTA in self.get_pipeline_parameters():
            self._input_fasta = self.get_pipeline_parameters()[self.CONFIG_INPUT_FASTA]
        elif self.CONFIG_KEY_GET_POSITION in self.get_default_parameters() and \
                self.CONFIG_INPUT_FASTA in self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION]:
            self._input_fasta = self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION][
                self.CONFIG_INPUT_FASTA]

        if self.CONFIG_OUTPUT_PSM_TABLE in self.get_pipeline_parameters():
            self._output_psm_table = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_PSM_TABLE]
        elif self.CONFIG_KEY_GET_POSITION in self.get_default_parameters() and \
                self.CONFIG_OUTPUT_PSM_TABLE in self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION]:
            self._output_psm_table = self.get_default_parameters()[self.CONFIG_KEY_GET_POSITION][
                self.CONFIG_OUTPUT_PSM_TABLE]
    
    def get_pep_pos(self, protein, sequence, fasta):
        for seq in SeqIO.parse(fasta, "fasta"):
            if seq.id == protein:
                sequence_canonical = str(seq.seq)
                rr = re.compile(sequence, re.I)
                match = re.finditer(rr, sequence_canonical)
                seq_position = int(re.sub("\D", "", protein.split(":")[2]))
                for i in match:
                    start = i.start() + 1
                    end = i.end()
                    if seq_position >= start and seq_position <= end:
                        return seq_position - start + 1
                    else:
                        continue

    def get_canonical_aa(self, protein):
        return protein.split(":")[2].split(".")[1][0]

    def get_variant_aa(self, protein):
        return protein.split(":")[2].split(".")[1][-1]

    def get_canonical_peptide(self, variant_peptide,canonical_aa,position):
        if position==0:
            return "Unmutated"
        else:
            str1 = variant_peptide[0:position - 1]
            str2 = variant_peptide[position:]
            return str1 + canonical_aa + str2
    
    def get_position(self, input_psm_table, input_fasta, output_psm_table):
        PSM = pd.read_table(input_psm_table, header = 0, dtype = "str", sep = "\t")
        PSM.loc[:, "Variant Peptide"] = PSM.apply(lambda x: re.sub("[^A-Z]","",x["Peptide"]), axis = 1)
        PSM.loc[:, "Canonical AA"] = PSM.apply(lambda x: self.get_canonical_aa(x["Protein"]), axis = 1)
        PSM.loc[:, "Variant AA"] = PSM.apply(lambda x: self.get_variant_aa(x["Protein"]), axis = 1)
        PSM.loc[:, "position"] = PSM.apply(lambda x: self.get_pep_pos(x["Protein"], x["Variant Peptide"], input_fasta), axis = 1)
        PSM["position"].fillna(0, inplace = True)

        PSM.loc[:, "Canonical Peptide"] = PSM.apply(lambda x: self.get_canonical_peptide(x["Variant Peptide"],x["Canonical AA"],int(x["position"])), axis = 1)
        psm = list(PSM)
        psm.insert(16,psm.pop(psm.index('Canonical Peptide')))
        PSM = PSM.loc[:,psm]

        PSM.to_csv(output_psm_table, sep = "\t", index = 0)
