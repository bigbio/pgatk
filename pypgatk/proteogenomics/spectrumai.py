import datetime
import os.path
import re
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Manager
from pyopenms import (TheoreticalSpectrumGenerator, MSSpectrum,
                      AASequence, Param, MzMLFile, MSExperiment, SpectrumLookup)
from tqdm import tqdm
from pypgatk.toolbox.general import ParameterConfiguration


class SpectrumAIService(ParameterConfiguration):
    CONFIG_KEY_VALIDATE_PEPTIDES = 'validate_peptides'
    CONFIG_MZML_PATH = 'mzml_path'
    CONFIG_MZML_FILES = 'mzml_files'
    CONFIG_INFILE_NAME = 'infile_name'
    CONFIG_OUTFILE_NAME = 'outfile_name'
    CONFIG_IONS_TOLERANCE = 'ions_tolerance'
    CONFIG_NUMBER_OF_PROCESSES = 'number_of_processes'
    CONFIG_RELATIVE = 'relative'
    CONFIG_MZTAB = 'mztab'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(SpectrumAIService, self).__init__(self.CONFIG_KEY_VALIDATE_PEPTIDES, config_data,
                                                pipeline_arguments)

        self._mzml_path = self.get_validate_parameters(variable=self.CONFIG_MZML_PATH, default_value=False)
        self._mzml_files = self.get_validate_parameters(variable=self.CONFIG_MZML_FILES, default_value=False)
        self._ions_tolerance = self.get_validate_parameters(variable=self.CONFIG_IONS_TOLERANCE, default_value=0.02)

        ## check if ions_tolerance is string, convert to float
        if isinstance(self._ions_tolerance, str):
            self._ions_tolerance = float(self._ions_tolerance)

        self._relative = self.get_validate_parameters(variable=self.CONFIG_RELATIVE, default_value=False)
        self._mztab = self.get_validate_parameters(variable=self.CONFIG_MZTAB, default_value=False)
        self._number_of_processes = self.get_validate_parameters(variable=self.CONFIG_NUMBER_OF_PROCESSES,
                                                                 default_value=40)

        self.df_list = Manager().list()

    def get_validate_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_VALIDATE_PEPTIDES in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_VALIDATE_PEPTIDES]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_VALIDATE_PEPTIDES][variable]
        return value_return

    def _predict_MS2_spectrum(self, peptide, size, product_ion_charge=1):

        tsg = TheoreticalSpectrumGenerator()
        spec = MSSpectrum()
        peptide = AASequence.fromString(peptide)

        p = Param()
        p.setValue("add_metainfo", "true")
        p.setValue("add_first_prefix_ion", "true")
        p.setValue("add_precursor_peaks", "true")
        tsg.setParameters(p)
        tsg.getSpectrum(spec, peptide, 1, 1)  # charge range 1:1

        b_y_ions = []
        for i in spec.getStringDataArrays()[0]:
            b_y_ions.append(i.decode())
        mz = []
        for i in spec:
            mz.append(i.getMZ())

        ions = pd.DataFrame({"mz": mz, "ion": b_y_ions, "z": 1})

        ions.loc[2 * size - 2, "ion"] = "b" + str(size)
        ions = ions.drop(2 * size - 1)
        ions.loc[2 * size, "ion"] = "y" + str(size)

        ions.loc[:, "ion"] = ions.apply(lambda x: re.sub("[+]", "", x["ion"]), axis=1)
        ions.loc[:, "pos"] = ions.apply(lambda x: re.sub(r"[^\d]", "", x["ion"]), axis=1)
        ions.loc[:, "type"] = ions.apply(lambda x: re.sub("[^a-z]", "", x["ion"]), axis=1)

        proton_mono_mass = 1.007276
        if product_ion_charge > 1:
            ions2 = ions.copy()
            ions2.loc[:, "mz"] = ions2.apply(lambda x: (x["mz"] + proton_mono_mass) / 2, axis=1)
            ions2.loc[:, "z"] = 2

            ions = ions.merge(ions2, how='outer')

        ions = ions.reset_index(drop=True)

        return ions

    @staticmethod
    def _get_intensity(exp_peak, ion_mz):
        exp_peak.loc[:, "mz_difference"] = exp_peak.apply(lambda x: abs(float(ion_mz) - x["mz"]), axis=1)
        min_index = exp_peak["mz_difference"].idxmin()
        return exp_peak.loc[exp_peak["mz_difference"] == exp_peak["mz_difference"].min()].loc[min_index, "intensity"]

    def _match_exp2predicted(self, exp_peak, pred_peak):
        pred_peak.loc[:, "error"] = pred_peak.apply(lambda x: min(abs(float(x["mz"]) - exp_peak["mz"])), axis=1)
        pred_peak.loc[:, "intensity"] = pred_peak.apply(lambda x: self._get_intensity(exp_peak, x["mz"]), axis=1)
        pred_peak.loc[:, "ppm"] = pred_peak.apply(lambda x: round(x["error"] / x["mz"] * 1000000, 2), axis=1)

        if self._relative:
            match_ions = pred_peak[pred_peak["ppm"] < self._ions_tolerance]
        else:
            match_ions = pred_peak[pred_peak["error"] < self._ions_tolerance]

        match_ions = match_ions.reset_index(drop=True)

        return match_ions

    def _inspect_spectrum(self, df, mzml_path, mzml_files):
        df.loc[:, "peptide_length"] = df.apply(lambda x: len(x["sequence"]), axis=1)

        df["status"] = "skiped"

        df["ions_support"] = "NO"
        df["support_ions"] = ""
        df["sum.supportions.intensity"] = float(0)

        df["flanking_ions_support"] = "NO"
        df["flanking_ions"] = ""
        df["sum.flanking.ions.intensity"] = float(0)

        df["matched_ions"] = ""
        df["sum.matchedions.intensity"] = float(0)
        df["sum.fragmentions.intensity"] = float(0)
        df["maxintensity"] = float(0)
        df["average_intensity"] = float(0)
        df["median_intensity"] = float(0)
        mzml_file = None

        if self._mztab:
            spectra_file = str(df.loc[0, "SpecFile"])
        else:
            spectra_file = str(df.loc[0, "reference_file_name"]) + ".mzML"

        if mzml_files and not mzml_path:
            mzml_list = mzml_files.split(",")
            for file in mzml_list:
                if spectra_file in file:
                    mzml_file = file
                    break
        elif not mzml_files and mzml_path:
            mzml_file = os.path.join(mzml_path, spectra_file)
        else:
            raise ValueError(
                "You only need to use either '--mzml_path' or '--mzml_files'.")

        exp = MSExperiment()
        try:
            MzMLFile().load(mzml_file, exp)
            look = SpectrumLookup()
            look.readSpectra(exp, r"((?<SCAN>)\d+$)")
        except Exception as e:
            print(mzml_file + " has ERROR!")
            print(e)
            df["ions_support"] = "mzML ERROR"
            return df

        for i in range(df.shape[0]):
            if self._mztab:
                scan_num = int(df.loc[i, "ScanNum"])
            else:
                scan_num = int(df.loc[i, "scan_number"])

            seq = df.loc[i, "sequence"]
            length = df.loc[i, "peptide_length"]

            # get peaks through ScanNum
            try:
                index = look.findByScanNumber(scan_num)
            except Exception as e:
                print("ERROR: " + str(e) + "; file:" + str(mzml_file) + "; scan_num:" + str(scan_num))
                continue

            exp_peaks = exp.getSpectrum(index).get_peaks()

            exp_peaks = pd.DataFrame({"mz": exp_peaks[0], "intensity": exp_peaks[1]})

            if self._mztab:
                predicted_peaks = self._predict_MS2_spectrum(
                    str(df.loc[i, "opt_global_cv_MS:1000889_peptidoform_sequence"]), length, 1)
            else:
                predicted_peaks = self._predict_MS2_spectrum(
                    str(df.loc[i, "peptidoform"]).replace("[","(").replace("]",")").replace("-",""), length, 1)
            match_ions = self._match_exp2predicted(exp_peaks, predicted_peaks)

            max_intensity = exp_peaks["intensity"].max()
            average_intensity = exp_peaks["intensity"].mean()
            median_intensity = exp_peaks["intensity"].median()

            df.loc[i, "sum.fragmentions.intensity"] = exp_peaks["intensity"].sum()
            df.loc[i, "maxintensity"] = max_intensity
            df.loc[i, "average_intensity"] = average_intensity
            df.loc[i, "median_intensity"] = median_intensity

            if match_ions.shape[0] == 0:
                continue
            df.loc[i, "matched_ions"] = ','.join(match_ions["ion"].unique().tolist())
            df.loc[i, "sum.matchedions.intensity"] = match_ions["intensity"].sum()

            if df.loc[i, "position"] == "canonical":
                continue
            if df.loc[i, "position"] == "non-canonical":
                continue
            position = int(df.loc[i, "position"])
            if position == 0:
                continue
            if position > length:
                continue

            df.loc[i, "status"] = "checked"
            supportions_intensity = 0
            ions_support = "NO"
            supportions = ""

            for j in range(match_ions.shape[0]):
                ion_type = match_ions.loc[j, "type"]
                pos = int(match_ions.loc[j, "pos"])
                ion = match_ions.loc[j, "ion"]

                if ion_type == "b" and pos >= position:
                    ions_support = "YES"
                    supportions_intensity = supportions_intensity + match_ions.loc[j, "intensity"]
                    supportions = supportions + ',' + ion
                elif ion_type == "y" and pos > length - position:
                    ions_support = "YES"
                    supportions_intensity = supportions_intensity + match_ions.loc[j, "intensity"]
                    supportions = supportions + ',' + ion

            df.loc[i, "ions_support"] = ions_support
            df.loc[i, "support_ions"] = supportions
            df.loc[i, "sum.supportions.intensity"] = supportions_intensity

            # check if it is a noise peak or isotope peak supporting mutant ions
            if df.loc[i, "sum.supportions.intensity"] < df.loc[i, "median_intensity"]:
                df.loc[i, "ions_support"] = "NO"

            flanking_ions_support = "NO"
            n1 = length
            n2 = position
            match_ions_set = set(match_ions["ion"].tolist())

            if n2 == 1:
                flanking_ions = {"b1", "y" + str(n1 - 1)}
                flanking_ions = flanking_ions.intersection(match_ions_set)
                if len(flanking_ions) > 0:
                    flanking_ions_support = "YES"
            elif n2 == n1:
                flanking_ions = {"y1", "b" + str(n1 - 1)}
                flanking_ions = flanking_ions.intersection(match_ions_set)
                if len(flanking_ions) > 0:
                    flanking_ions_support = "YES"
            else:
                flanking_ions_left = {"b" + str(n2 - 1), "y" + str(n1 - n2 + 1)}
                flanking_ions_right = {"b" + str(n2), "y" + str(n1 - n2)}

                flanking_ions_left = flanking_ions_left.intersection(match_ions_set)
                flanking_ions_right = flanking_ions_right.intersection(match_ions_set)

                flanking_ions = flanking_ions_left.union(flanking_ions_right)
                if len(flanking_ions_left) > 0 and len(flanking_ions_right) > 0:
                    flanking_ions_support = "YES"

            df.loc[i, "flanking_ions_support"] = flanking_ions_support
            df.loc[i, "flanking_ions"] = ",".join(flanking_ions)
            if flanking_ions:
                df.loc[i, "sum.flanking.ions.intensity"] = \
                    match_ions[match_ions['ion'].str.contains("|".join(flanking_ions))]["intensity"].sum()

            if df.loc[i, "sum.flanking.ions.intensity"] < df.loc[i, "median_intensity"]:
                df.loc[i, "flanking_ions_support"] = "NO"

            # fragmentation is not preferable at Cterm side of proline, so only require supporting ions
            if re.search("P", seq[position - 1:position]):
                df.loc[i, "flanking_ions_support"] = df.loc[i, "ions_support"]

        return df

    def _multiprocess_inspect_spectrum(self, df):
        self.df_list.append(self._inspect_spectrum(df, self._mzml_path, self._mzml_files))

    def validate(self, infile_name, outfile_name: str):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)

        if infile_name.endswith(".csv.gz"):
            df_psm = pd.read_csv(infile_name, header=0, sep=",", compression="gzip")
        elif infile_name.endswith(".csv"):
            df_psm = pd.read_csv(infile_name, header=0, sep=",")
        elif infile_name.endswith(".tsv.gz"):
            df_psm = pd.read_table(infile_name, header=0, sep="\t", compression="gzip")
        elif infile_name.endswith(".tsv"):
            df_psm = pd.read_table(infile_name, header=0, sep="\t")
        else:
            raise ValueError("The input file format is not supported.")

        if self._mztab:
            grouped_dfs = df_psm.groupby("SpecFile")
        else:
            grouped_dfs = df_psm.groupby("reference_file_name")

        list_of_dfs = [group_df.reset_index(drop=True) for name, group_df in grouped_dfs]

        pool = Pool(int(self._number_of_processes))
        list(tqdm(pool.imap(self._multiprocess_inspect_spectrum, list_of_dfs), total=len(list_of_dfs),
                  desc="Validate By Each mzMl", unit="mzML"))
        pool.close()
        pool.join()

        df_output = pd.concat(self.df_list, axis=0, ignore_index=True)

        if outfile_name.endswith(".csv.gz"):
            df_output.to_csv(outfile_name, header=True, sep=",", index=None, compression="gzip")
        elif outfile_name.endswith(".csv"):
            df_output.to_csv(outfile_name, header=True, sep=",", index=None)
        elif outfile_name.endswith(".tsv.gz"):
            df_output.to_csv(outfile_name, header=True, sep="\t", index=None, compression="gzip")
        else:  # .tsv or other formats
            df_output.to_csv(outfile_name, header=True, sep="\t", index=None)

        end_time = datetime.datetime.now()
        print("End time :", end_time)
        time_taken = end_time - start_time
        print("Time consumption :", time_taken)
