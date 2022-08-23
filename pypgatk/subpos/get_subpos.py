import pandas as pd
import re
from Bio import SeqIO

from pypgatk.toolbox.general import ParameterConfiguration

class GetSubPos(ParameterConfiguration):
    CONFIG_KEY_GET_SUBPOS = 'get_subpos'
    CONFIG_INPUT_PSM_TABLE = 'input_psm_table'
    CONFIG_INPUT_FASTA = 'input_fasta'
    CONFIG_OUTPUT_PSM_TABLE = 'output_psm_table'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(GetSubPos, self).__init__(self.CONFIG_KEY_GET_SUBPOS, config_data, pipeline_arguments)

        if self.CONFIG_INPUT_PSM_TABLE in self.get_pipeline_parameters():
            self._input_psm_table = self.get_pipeline_parameters()[self.CONFIG_INPUT_PSM_TABLE]
        elif self.CONFIG_KEY_GET_SUBPOS in self.get_default_parameters() and \
                self.CONFIG_INPUT_PSM_TABLE in self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS]:
            self._input_psm_table = self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS][
                self.CONFIG_INPUT_PSM_TABLE]
        
        if self.CONFIG_INPUT_FASTA in self.get_pipeline_parameters():
            self._input_fasta = self.get_pipeline_parameters()[self.CONFIG_INPUT_FASTA]
        elif self.CONFIG_KEY_GET_SUBPOS in self.get_default_parameters() and \
                self.CONFIG_INPUT_FASTA in self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS]:
            self._input_fasta = self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS][
                self.CONFIG_INPUT_FASTA]

        if self.CONFIG_OUTPUT_PSM_TABLE in self.get_pipeline_parameters():
            self._output_psm_table = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_PSM_TABLE]
        elif self.CONFIG_KEY_GET_SUBPOS in self.get_default_parameters() and \
                self.CONFIG_OUTPUT_PSM_TABLE in self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS]:
            self._output_psm_table = self.get_default_parameters()[self.CONFIG_KEY_GET_SUBPOS][
                self.CONFIG_OUTPUT_PSM_TABLE]
    
    def find_position(self, target):
        result = target.split(":")[2].split(".")[1][1:]
        num = ""
        for i in range(len(result)):
            if result[i].isdigit():
                num += result[i]
            else:
                break
        return num
    
    def sub_pos(self, target, target_pep, pro_pos, fasta):
        for seq_record in SeqIO.parse(fasta, "fasta"):
            if seq_record.id == target:
                seq = re.sub('\(|\)|\,', '', str(seq_record.seq))

                rr = re.compile(target_pep, re.I)
                match = re.finditer(rr, seq)
                pro_pos = int(pro_pos)
                for i in match:
                    start = i.start() + 1
                    end = i.end() + 1
                    if pro_pos >= start and pro_pos <= end:
                        return pro_pos - start + 1
                    else: return 0
    
    def get_subpos(self, input_psm_table, input_fasta, output_psm_table):
        psm_table = pd.read_csv(input_psm_table, header = 0, dtype = "str", sep = "\t")
        psm_table.loc[:, "Sequence"] = psm_table.apply(lambda x: re.sub("[^A-Z]","",x["Peptide"]), axis = 1)
        psm_table.loc[:, "pro_pos"] = psm_table.apply(lambda x: self.find_position(x["Protein"]), axis = 1, result_type = "expand")
        psm_table.loc[:, "sub_pos"] = psm_table.apply(lambda x: self.sub_pos(x["Protein"], x["Sequence"], x["pro_pos"], input_fasta), axis = 1, result_type = "expand")
        psm_table["sub_pos"].fillna(0, inplace = True)
        psm_table.to_csv(output_psm_table, sep = "\t", index = 0)
