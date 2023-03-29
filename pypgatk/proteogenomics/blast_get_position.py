import pandas as pd
from Bio import pairwise2,SeqIO
import datetime
from pathos.multiprocessing import ProcessingPool as Pool

from pypgatk.toolbox.general import ParameterConfiguration

class BlastGetPositionService(ParameterConfiguration):
    CONFIG_KEY_BlastGetPosition = 'blast_get_position'
    CONFIG_CANONICAL_PEPTIDE_PREFIX = 'canonical_peptide_prefix'
    CONFIG_INPUT_REFERENCE_DATABASE = 'input_reference_database'
    CONFIG_NUMBER_OF_PROCESSES = 'number_of_processes'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(BlastGetPositionService, self).__init__(self.CONFIG_KEY_BlastGetPosition, config_data, pipeline_arguments)
        self._canonical_peptide_prefix = self.get_blast_parameters(variable=self.CONFIG_CANONICAL_PEPTIDE_PREFIX, default_value='sp,NP,ENSP')
        self._input_reference_database = self.get_blast_parameters(variable=self.CONFIG_INPUT_REFERENCE_DATABASE, default_value='')
        self._number_of_processes = self.get_blast_parameters(variable=self.CONFIG_NUMBER_OF_PROCESSES, default_value='40')

        self.fa_set =set()
        for j in SeqIO.parse(self._input_reference_database, "fasta"):
            self.fa_set.add(str(j.seq))

    def get_blast_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_BlastGetPosition in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_BlastGetPosition]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_BlastGetPosition][variable]
        return value_return

    def _is_canonical(self, accession):
        prefix_list = self._canonical_peptide_prefix.split(",")
        accession_list = accession.split(",")
        for i in accession_list:
            for j in prefix_list:
                if j in i:
                    return "canonical"
        return "waiting for blast"
    
    def _blast_set(self, fasta_set, peptide):
        length = len(peptide)
        position_set = set()
        for fasta in fasta_set:
            if len(fasta)>=length:
                alignments_score  = pairwise2.align.localms(sequenceA = fasta,sequenceB = peptide,match = 1,mismatch = 0,open = -1,extend = 0,score_only = True)
                if alignments_score == length:
                    return "canonical"
                elif alignments_score == length-1:
                    alignments_local  = pairwise2.align.localms(sequenceA = fasta,sequenceB = peptide,match = 1,mismatch = 0,open = -1,extend = 0)
                    for alignment in alignments_local:
                        #insertion e.g., ABCDMEFGH<----ABCDEFGH
                        if alignment.end - alignment.start == length+1:
                            s = fasta[alignment.start:alignment.end]
                            for i in range(length):
                                if peptide[i] != s[i]:
                                    position_set.add(i+1)
                                    break
                        #substitution e.g., ABCDMFGH<----ABCDEFGH
                        elif alignment.end - alignment.start == length:
                            s = fasta[alignment.start:alignment.end]
                            for i in range(length):
                                if peptide[i] != s[i]:
                                    position_set.add(i+1)
                                    break
                        # substitution e.g., ABCDEFGM<----ABCDEFGH
                        elif alignment.end - alignment.start == length-1:
                            s = fasta[alignment.start:alignment.end]
                            if peptide[0] != s[0]:
                                position_set.add(1)
                            elif peptide[-1] != s[-1]:
                                position_set.add(length)
                elif alignments_score == length-2:
                    alignments_local = pairwise2.align.localms(sequenceA=fasta, sequenceB=peptide, match=1, mismatch=-1,
                                                            open=-1, extend=0)
                    for alignment in alignments_local:
                        # deletion e.g., ABCEFGH<----ABCDEFGH
                        if alignment.end - alignment.start == length and alignment.score == length - 2:
                            s = fasta[alignment.start:alignment.end - 1]
                            if pairwise2.align.localms(sequenceA=s, sequenceB=peptide, match=1, mismatch=0, open=0,
                                                    extend=0, score_only=True) == length - 1:
                                for i in range(length - 1):
                                    if peptide[i] != s[i]:
                                        position_set.add(i+1)
                                        break
        if position_set:
            return position_set
        else:
            return "non-canonical"

    def _result(self, sequence):
        return self._blast_set(self.fa_set ,sequence)

    def blast(self, input_psm_to_blast, output_psm):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)

        PSM = pd.read_table(input_psm_to_blast, header=0,sep="\t")
        PSM.loc[:,"position"] = PSM.apply(lambda x: self._is_canonical(x["accession"]),axis = 1)

        first_filter = PSM[PSM.position=="canonical"]
        psm_to_blast = PSM[PSM.position=="waiting for blast"]
        psm_to_blast = psm_to_blast.copy()

        # seq = psm_to_blast["sequence"].to_list()
        # pool = Pool()
        # res = pool.map(self._result, seq)
        # psm_to_blast["position"] = res
        # pool.close()

        #Remove duplicate sequences
        seq_set = set(psm_to_blast["sequence"].to_list())
        seq_list = list(seq_set)

        pool = Pool(int(self._number_of_processes))
        res = pool.map(self._result, seq_list)
        seq_dict = dict()
        for i in range(len(seq_list)):
            seq_dict[seq_list[i]] = res[i]
        psm_to_blast["position"] = PSM.apply(lambda x: seq_dict.get(x["sequence"]),axis=1)
        pool.close()

        second_filter = psm_to_blast[psm_to_blast.position=="canonical"]
        non_filter = psm_to_blast[psm_to_blast.position=="non-canonical"]

        psm_to_findpos = psm_to_blast[psm_to_blast.position!="canonical"]
        psm_to_findpos = psm_to_findpos[psm_to_findpos.position!="non-canonical"]

        if len(psm_to_findpos)>0:
            psm_to_findpos["var_num"] = psm_to_findpos.apply(lambda x: len(x["position"]),axis=1)
            psm_to_findpos = psm_to_findpos.loc[psm_to_findpos.index.repeat(psm_to_findpos["var_num"])]
            psm_to_findpos["var_num"].iloc[0] = 0
            psm_id = psm_to_findpos["PSM_ID"].iloc[0]
            for i in range(1,psm_to_findpos.shape[0]):
                if psm_to_findpos["PSM_ID"].iloc[i] == psm_id:
                    psm_to_findpos["var_num"].iloc[i] = psm_to_findpos["var_num"].iloc[i-1]+1
                else:
                    psm_id = psm_to_findpos["PSM_ID"].iloc[i]
                    psm_to_findpos["var_num"].iloc[i] = 0
            psm_to_findpos["position"] = psm_to_findpos.apply(lambda x: str(x["position"])[1:-1].split(",")[x["var_num"]],axis=1)
            psm_to_findpos.drop(columns="var_num", axis = 1, inplace=True)
            psm_to_findpos["position"] = psm_to_findpos.apply(lambda x: x["position"].replace(' ', ''), axis=1)

        all_psm_out = pd.concat([first_filter, second_filter, non_filter, psm_to_findpos], axis=0, join='outer')
        all_psm_out = all_psm_out.sort_values("PSM_ID")
        all_psm_out.to_csv(output_psm, header=1,sep="\t", index=None)

        end_time = datetime.datetime.now()
        print("End time :", end_time)
        set_time_taken = end_time - start_time
        print("Time consumption :", set_time_taken)



