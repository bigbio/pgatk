import string

import pandas as pd
from pyspark.sql.types import StructType, StructField, StringType

# SCHEMAS for MZTAB
MZTAB_PSMS_COLUMNS = ["sequence",	"PSM_ID",	"accession",	"unique",	"database",
                      "database_version",	"search_engine	search_engine_score[1]",
                      "modifications",	"retention_time",	"charge",	"exp_mass_to_charge",
                      "calc_mass_to_charge"	"spectra_ref", "pre",	"post",	"start"]

class MzTabColumn:
  def __init__(self, name: str=None, pandas_type =string, spark_type = StringType()) -> None:
    self._name = name
    self._pandas_type = pandas_type
    self._spark_type = spark_type
