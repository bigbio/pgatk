import pandas as pd
from Bio import pairwise2, SeqIO
import datetime
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Manager
from tqdm import tqdm
import ahocorasick

from pypgatk.toolbox.general import ParameterConfiguration

def get_details(fasta, peptide):
    res = []
    i = 0
    j = 0
    for AA1, AA2 in zip(fasta, peptide):
        i += 1
        j += 1
        if AA1 == AA2:
            continue
        else:
            res.append(str(i) + "|" + AA1 + ">" + AA2)
    return res

def peptide_blast_protein(fasta, peptide):
    length = len(peptide)
    mismatch = []
    if len(fasta) >= length:
        score = pairwise2.align.localms(sequenceA=fasta, sequenceB=peptide,
                                                       match=1, mismatch=0, open=-2, extend=-2, score_only=True)
        if score == length-1:
            alignment = pairwise2.align.localms(sequenceA=fasta, sequenceB=peptide,
                                                       match=1, mismatch=0, open=-2, extend=-2)[0]
            if alignment.end - alignment.start == length:
                mismatch = get_details(alignment.seqA[alignment.start:alignment.end], alignment.seqB[alignment.start:alignment.end])
            elif alignment.end - alignment.start == length-1:
                res = get_details(alignment.seqA[alignment.start:alignment.end+1], alignment.seqB[alignment.start:alignment.end+1])
                if len(res) == 1:
                    if res[0].split(">")[1]!="-":
                        mismatch = res
                    else:
                        mismatch = get_details(alignment.seqA[alignment.start-1:alignment.end], alignment.seqB[alignment.start-1:alignment.end])
                elif len(res) == 0:
                    mismatch = get_details(alignment.seqA[alignment.start-1:alignment.end], alignment.seqB[alignment.start-1:alignment.end])
                else:
                    print("Number of mismatch Error")
    return mismatch

def _blast_set(fasta_dict, peptide):
    positions = dict()
    for fasta in fasta_dict.keys():
        mismatch = peptide_blast_protein(fasta, peptide)
        if len(mismatch) == 1:
            if positions.get(mismatch[0]):
                positions[mismatch[0]].update(fasta_dict[fasta])
            else:
                positions[mismatch[0]] = fasta_dict[fasta]
        elif len(mismatch) > 1:
            print("Number of mismatch > 1")
            print(peptide)
            print(fasta)
            print(mismatch)
    if positions:
        res = []
        for key,value in positions.items():
            splits = key.split("|")
            splits.append(",".join(value))
            res.append(splits)
        return res
    else:
        return "non-canonical"

class BlastGetPositionService(ParameterConfiguration):
    CONFIG_KEY_BlastGetPosition = 'blast_get_position'
    CONFIG_INPUT_REFERENCE_DATABASE = 'input_reference_database'
    CONFIG_NUMBER_OF_PROCESSES = 'number_of_processes'

    def __init__(self, config_data, pipeline_arguments):
        """
        init the class with the specific parameters.
        :param config_data configuration file
        :param pipeline_arguments pipelines arguments
        """

        super(BlastGetPositionService, self).__init__(self.CONFIG_KEY_BlastGetPosition, config_data, pipeline_arguments)
        self._input_reference_database = self.get_blast_parameters(variable=self.CONFIG_INPUT_REFERENCE_DATABASE,
                                                                   default_value='')
        self._number_of_processes = self.get_blast_parameters(variable=self.CONFIG_NUMBER_OF_PROCESSES,
                                                              default_value='40')


        self.fasta_dict = dict()
        for j in SeqIO.parse(self._input_reference_database, "fasta"):
            if self.fasta_dict.get(str(j.seq)):
                self.fasta_dict[str(j.seq)].add(j.id)
            else:
                self.fasta_dict[str(j.seq)] = {j.id}
        self.blast_dict = Manager().dict()

    def get_blast_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_BlastGetPosition in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_BlastGetPosition]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_BlastGetPosition][variable]
        return value_return

    def _blast_canonical(self, df):
        seq_set = set(df["sequence"].to_list())

        auto = ahocorasick.Automaton()
        seq_dict = dict()
        for seq_peptide in seq_set:
            auto.add_word(seq_peptide, seq_peptide)
            seq_dict[seq_peptide] = "waiting for blast"

        auto.make_automaton()

        for protein_seq in self.fasta_dict.keys():
            for end_ind, found in auto.iter(protein_seq):
                seq_dict[found] = "canonical"
                print("Found", found, "at position", end_ind, "in protein sequence")

        df["position"] = df["sequence"].map(seq_dict)
        return df

    def _result(self, sequence):
        self.blast_dict[sequence] = _blast_set(self.fasta_dict, sequence)

    def blast(self, input_psm_to_blast, output_psm):
        """
        Blast peptide and reference protein database to find variation sites.
        :param input_psm_to_blast: input PSM table to blast
        :param output_psm: output PSM table
        :return:
        """

        start_time = datetime.datetime.now()
        print("Start time :", start_time)

        if input_psm_to_blast.endswith(".csv.gz"):
            psm = pd.read_csv(input_psm_to_blast, header=0, sep=",", compression="gzip")
        elif input_psm_to_blast.endswith(".csv"):
            psm = pd.read_csv(input_psm_to_blast, header=0, sep=",")
        elif input_psm_to_blast.endswith(".tsv.gz"):
            psm = pd.read_table(input_psm_to_blast, header=0, sep="\t", compression="gzip")
        elif input_psm_to_blast.endswith(".tsv"):
            psm = pd.read_table(input_psm_to_blast, header=0, sep="\t")
        else:
            raise ValueError("The input file format is not supported.")


        psm = self._blast_canonical(psm)

        first_filter = psm[psm.position == "canonical"]
        psm_to_blast = psm[psm.position == "waiting for blast"]
        psm_to_blast = psm_to_blast.copy()

        # Remove duplicate sequences
        seq_set = set(psm_to_blast["sequence"].to_list())
        seq_list = list(seq_set)

        pool = Pool(int(self._number_of_processes))
        list(tqdm(pool.imap(self._result, seq_list), total=len(seq_list), desc="Blast", unit="peptide"))

        pool.close()
        pool.join()

        psm_to_blast["position"] = psm.apply(lambda x: self.blast_dict.get(x["sequence"]), axis=1)

        second_filter = psm_to_blast[psm_to_blast.position == "canonical"]
        non_filter = psm_to_blast[psm_to_blast.position == "non-canonical"]

        psm_to_findpos = psm_to_blast[psm_to_blast.position != "canonical"]
        psm_to_findpos = psm_to_findpos[psm_to_findpos.position != "non-canonical"]

        if len(psm_to_findpos) > 0:
            psm_to_findpos = psm_to_findpos.explode("position", ignore_index=True)
            psm_to_findpos["variant"] = psm_to_findpos["position"].apply(lambda x: x[1])
            psm_to_findpos["protein"] = psm_to_findpos["position"].apply(lambda x: x[2])
            psm_to_findpos["position"] = psm_to_findpos["position"].apply(lambda x: x[0])

        all_psm_out = pd.concat([first_filter, second_filter, non_filter, psm_to_findpos], axis=0, join='outer')
        all_psm_out = all_psm_out.sort_values("usi")

        if output_psm.endswith(".csv.gz"):
            all_psm_out.to_csv(output_psm, header=True, sep=",", index=None, compression="gzip")
        elif output_psm.endswith(".csv"):
            all_psm_out.to_csv(output_psm, header=True, sep=",", index=None)
        elif output_psm.endswith(".tsv.gz"):
            all_psm_out.to_csv(output_psm, header=True, sep="\t", index=None, compression="gzip")
        elif output_psm.endswith(".tsv"):
            all_psm_out.to_csv(output_psm, header=True, sep="\t", index=None)
        else:
            all_psm_out.to_csv(output_psm, header=True, sep="\t", index=None)

        end_time = datetime.datetime.now()
        print("End time :", end_time)
        set_time_taken = end_time - start_time
        print("Time consumption :", set_time_taken)
