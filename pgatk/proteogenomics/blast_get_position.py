import datetime
import functools
import logging

import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import ahocorasick

from pgatk.toolbox.general import ParameterConfiguration

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
        for i in range(len(fasta) - length + 1):
            window = fasta[i:i + length]
            details = get_details(window, peptide)
            if len(details) == 1:
                mismatch = details
                break
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
            logging.getLogger(__name__).warning(
                "Number of mismatch > 1: peptide=%s, fasta=%s, mismatch=%s", peptide, fasta, mismatch)
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
        self.blast_dict = {}

    def get_blast_parameters(self, variable: str, default_value):
        return self.get_config_value(variable, default_value)

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
                self.get_logger().info("Found %s at position %s in protein sequence", found, end_ind)

        df["position"] = df["sequence"].map(seq_dict)
        return df

    @staticmethod
    def _result_worker(fasta_dict, sequence):
        """Process a single peptide and return (sequence, blast_result)."""
        return (sequence, _blast_set(fasta_dict, sequence))

    def blast(self, input_psm_to_blast, output_psm):
        """
        Blast peptide and reference protein database to find variation sites.
        :param input_psm_to_blast: input PSM table to blast
        :param output_psm: output PSM table
        :return:
        """

        start_time = datetime.datetime.now()
        self.get_logger().info("Start time: %s", start_time)

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

        worker = functools.partial(BlastGetPositionService._result_worker, self.fasta_dict)
        with ProcessPoolExecutor(max_workers=int(self._number_of_processes)) as executor:
            results = list(tqdm(executor.map(worker, seq_list, chunksize=100),
                                total=len(seq_list), desc="Blast", unit="peptide"))

        for sequence, blast_result in results:
            self.blast_dict[sequence] = blast_result

        psm_to_blast["position"] = psm_to_blast["sequence"].map(self.blast_dict.get)

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
        self.get_logger().info("End time: %s", end_time)
        set_time_taken = end_time - start_time
        self.get_logger().info("Time consumption: %s", set_time_taken)
