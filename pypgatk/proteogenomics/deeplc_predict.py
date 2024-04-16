import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from deeplc import DeepLC
import datetime
# from deeplcretrainer import deeplcretrainer


from pypgatk.toolbox.general import ParameterConfiguration

class DeepLCPredictService(ParameterConfiguration):
    CONFIG_KEY_DeepLCPredict = 'deeplc_predict'
    CONFIG_TRANSFER_LEARNING = 'transfer_learning'
    CONFIG_REDUNDANCY_REMOVAL_STRATEGY = 'redundancy_removal_strategy'
    CONFIG_NUMBER_OF_PROCESSES = 'number_of_processes'
    CONFIG_FILTRATION_RATIO = 'filtration_ratio'

    def __init__(self, config_data, pipeline_arguments):
        """
      Init the class with the specific parameters.
      :param config_data configuration file
      :param pipeline_arguments pipelines arguments
      """

        super(DeepLCPredictService, self).__init__(self.CONFIG_KEY_DeepLCPredict, config_data, pipeline_arguments)
        self._transfer_learning = self.get_deeplc_parameters(variable=self.CONFIG_TRANSFER_LEARNING,default_value=0)
        self._redundancy_removal_strategy = self.get_deeplc_parameters(variable=self.CONFIG_REDUNDANCY_REMOVAL_STRATEGY, default_value='median')
        self._number_of_processes = int(self.get_deeplc_parameters(variable=self.CONFIG_NUMBER_OF_PROCESSES, default_value=40))
        self._filtration_ratio = self.get_deeplc_parameters(variable=self.CONFIG_FILTRATION_RATIO, default_value=0.01)

        self.mod_rep = {"UNIMOD:4": "Carbamidomethyl", "UNIMOD:737": "TMT6plex", "UNIMOD:35": "Oxidation", "UNIMOD:1": "Acetyl"}
    
    def get_deeplc_parameters(self, variable: str, default_value):
        value_return = default_value
        if variable in self.get_pipeline_parameters():
            value_return = self.get_pipeline_parameters()[variable]
        elif self.CONFIG_KEY_DeepLCPredict in self.get_default_parameters() and \
                variable in self.get_default_parameters()[self.CONFIG_KEY_DeepLCPredict]:
            value_return = self.get_default_parameters()[self.CONFIG_KEY_DeepLCPredict][variable]
        return value_return

    def _replace_mod(self, x, mod_dict):
        s = ""
        for mod in x.split(','):
            nums = mod.split("-")[0]
            m = mod.split("-")[1]
            if s:
                    s += "|"
            s += nums + "|" + mod_dict.get(m)
        return s
    
    def predict(self, input_to_deeplc, output):
        start_time = datetime.datetime.now()
        print("Start time :", start_time)

        validate_out = pd.read_table(input_to_deeplc, header=0, sep="\t", index_col=0)
        
        validate_out.loc[:, "seq"] = validate_out.apply(lambda x: x["sequence"], axis=1)
        validate_out["modifications"] = validate_out["modifications"].fillna("")
        validate_out.loc[:, "modifications"] = validate_out.apply(lambda x: self._replace_mod(x["modifications"], self.mod_rep) if x["modifications"] != "" else "", axis=1)
        validate_out.loc[:, "tr"] = validate_out.apply(lambda x: x["retention_time"], axis=1)

        filter_out = validate_out[(validate_out["position"] == "non-canonical") | (validate_out["position"] == "canonical") | (validate_out["flanking_ions_support"] == "YES")]
        print("all(no_drop_by_PSMid):" + str(len(filter_out)))
        filter_out = filter_out.drop_duplicates(['PSM_ID'])

        can = filter_out[filter_out["position"]=="canonical"]
        non = filter_out[filter_out["position"]!="canonical"]
        print("can(after_drop_by_PSMid):" + str(len(can)))
        print("non(after_drop_by_PSMid):" + str(len(non)))

        can_after_drop = can
        non_after_drop = non
        if self._redundancy_removal_strategy == 'median':
            can_after_drop = can.groupby(["seq","modifications"]).median().reset_index()
            non_after_drop = non.groupby(["seq","modifications"]).median().reset_index()
        elif self._redundancy_removal_strategy == 'mean':
            can_after_drop = can.groupby(["seq","modifications"]).mean().reset_index()
            non_after_drop = non.groupby(["seq","modifications"]).mean().reset_index()
        elif self._redundancy_removal_strategy == 'first':
            can_after_drop = can.sort_values("tr").drop_duplicates(["seq","modifications"])
            non_after_drop = non.sort_values("tr").drop_duplicates(["seq","modifications"])
        else:
            raise ValueError(
                "You need to specify a strategy for removing redundancy.")    
        print("can(after_drop_by_same_sequence+modification):" + str(len(can_after_drop)))
        print("non(after_drop_by_same_sequence+modification):" + str(len(non_after_drop)))

        if self._transfer_learning == 1:
            dlc = DeepLC(
                deeplc_retrain=True,
                n_epochs=30,
                n_jobs=self._number_of_processes,
            )
        else:
            dlc = DeepLC(
                pygam_calibration=True,
                n_jobs=self._number_of_processes,
            )
        
        dlc.calibrate_preds(seq_df=can_after_drop)
        preds_calib = dlc.make_preds(seq_df=non_after_drop)

        non_after_drop["preds_tr"] = preds_calib
        non_after_drop.index = non_after_drop["seq"]+"+"+non_after_drop["modifications"]
        pred_dict = non_after_drop["preds_tr"].to_dict()

        non.index = non["seq"]+"+"+non["modifications"]
        non["preds_tr"] = [pred_dict[v] for v in non.index]

        non["error"] = non["tr"]-non["preds_tr"]
        non["abserror"] = abs(non["error"])

        plt.hist(non["abserror"], bins=1000)
        y_max = plt.gca().get_ylim()[1]
        plt.vlines(np.percentile(non["abserror"], 95), 0, y_max, color="grey", label="95th percentile")
        plt.vlines(np.percentile(non["abserror"], 99), 0, y_max, color="black", label="99th percentile")
        plt.vlines(np.percentile(non["abserror"], (1-self._filtration_ratio)*100), 0, y_max, color="red", label= str((1-self._filtration_ratio)*100) + "th percentile")
        plt.legend()
        plt.show()
        plt.savefig('percentile.png')

        non_after_filter = non[non["abserror"] < np.percentile(non["abserror"], (1-self._filtration_ratio)*100)]
        print("non(after_filter):" + str(len(non_after_filter)))

        all = pd.concat([can, non_after_filter], axis=0)
        all.to_csv(output,sep="\t",index = None)

        end_time = datetime.datetime.now()
        print("End time :", end_time)
        time_taken = end_time - start_time
        print("Time consumption :", time_taken)






        
            


        
        