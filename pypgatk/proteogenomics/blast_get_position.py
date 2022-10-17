import pandas as pd
from Bio import pairwise2,SeqIO
import datetime

from pypgatk.toolbox.general import ParameterConfiguration

class BlastGetPositonService(ParameterConfiguration):
    CONFIG_KEY_BlastGetPositon = 'blast_get_positon'
    CONFIG_CANONICAL_PEPTIDE_PREFIX = 'canonical_peptide_prefix'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(BlastGetPositonService, self).__init__(self.CONFIG_KEY_BlastGetPositon, config_data, pipeline_arguments)

        self._canonical_peptide_prefix = self.get_blast_parameters(variable=self.CONFIG_CANONICAL_PEPTIDE_PREFIX, default_value='sp,NP,ENSP')

    def get_blast_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_BlastGetPositon in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_BlastGetPositon]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_BlastGetPositon][variable]
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

    def blast(self, input_psm_to_blast, output_psm, input_refence_proteome):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)

        PSM = pd.read_table(input_psm_to_blast, header=0,sep="\t")
        PSM.loc[:,"position"] = PSM.apply(lambda x: self._is_canonical(x["accession"]),axis = 1)
        first_filter = PSM[PSM.position=="canonical"]
        psm_to_blast = PSM[PSM.position=="waiting for blast"]

        fa_set =set()
        for j in SeqIO.parse(input_refence_proteome, "fasta"):
            fa_set.add(str(j.seq))
        psm_to_blast.loc[:, "position"] = psm_to_blast.apply(lambda x: self._blast_set(fa_set, x["sequence"]), axis=1)

        second_filter = psm_to_blast[psm_to_blast.position=="canonical"]
        non_filter = psm_to_blast[psm_to_blast.position=="non-canonical"]

        psm_to_findpos = psm_to_blast[psm_to_blast.position!="canonical"]
        psm_to_findpos = psm_to_findpos[psm_to_findpos.position!="non-canonical"]

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
        print("Set Time consumption :", set_time_taken)



