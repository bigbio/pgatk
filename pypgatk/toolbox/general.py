"""
This module implements some useful functions for the pipeline runner.

@author ypriverol
"""
import logging
import os
import json
import shutil
import subprocess
from typing import List
from urllib import request
from urllib.error import URLError, ContentTooShortError
import yaml
import gzip

# Logging defaults
from requests import HTTPError

from pypgatk.toolbox.exceptions import ToolBoxException

REMAINING_DOWNLOAD_TRIES = 4

class ParameterConfiguration:
  """
    This class is a helper class for those submodules having to manage configuration files themselves, that are specific
    to them
    """

  _logger_formatters = {
    "DEBUG": "%(asctime)s [%(levelname)7s][%(name)28s][%(module)18s, %(lineno)4s] %(message)s",
    "INFO": "%(asctime)s [%(levelname)7s][%(name)28s] %(message)s"
  }
  _log_level = 'DEBUG'

  _CONFIG_LOGGER = 'logger'
  _CONFIG_LOGGER_FORMATTER = 'formatters'
  _CONFIG_LOGGER_LEVEL = 'loglevel'

  def __init__(self, root_config_name, yaml_configuration, pipeline_parameters):
    """
        This function creates a parameter structure from a yaml config file and the pipeline paramters provdided
        in the commandline
        :param yaml_configuration_file: yaml config file
        :param pipeline_parameters: commandline parameters.
        """

    self._ROOT_CONFIG_NAME = root_config_name

    self._default_params = yaml_configuration
    self._pipeline_parameters = pipeline_parameters

    # Prepare Logging subsystem
    if self._ROOT_CONFIG_NAME in self._default_params:
      if self._CONFIG_LOGGER in self._default_params[self._ROOT_CONFIG_NAME]:
        if self._CONFIG_LOGGER_LEVEL in self._default_params[self._ROOT_CONFIG_NAME][self._CONFIG_LOGGER]:
          self._log_level = self._default_params[self._ROOT_CONFIG_NAME][self._CONFIG_LOGGER][
            self._CONFIG_LOGGER_LEVEL]
        if self._CONFIG_LOGGER_FORMATTER in self._default_params[self._ROOT_CONFIG_NAME][self._CONFIG_LOGGER]:
          self._logger_formatters = self._default_params[self._ROOT_CONFIG_NAME][self._CONFIG_LOGGER][
            self._CONFIG_LOGGER_FORMATTER]

    self._log_handlers = []
    log_handlers_prefix = self._ROOT_CONFIG_NAME + '-'
    log_handlers_extension = '.log'

    self._logger = logging.getLogger(__name__)
    self._logger.setLevel(getattr(logging, self._log_level))
    self._log_files = []
    for llevel, lformat in self._logger_formatters.items():
      logfile = os.path.join(log_handlers_prefix + llevel.lower() + log_handlers_extension)
      lformatter = logging.Formatter(lformat)
      lhandler = logging.FileHandler(logfile, mode='w')
      lhandler.setLevel(getattr(logging, llevel))
      lhandler.setFormatter(lformatter)
      self._log_handlers.append(lhandler)
      # Add the handlers to my own logger
      self._logger.addHandler(lhandler)
      # Keep the path to the log file
      self._log_files.append(logfile)
    self.get_logger().debug("Logging system initialized")

  def get_pipeline_parameters(self):
    return self._pipeline_parameters

  def get_default_parameters(self):
    return self._default_params

  def get_log_handlers(self):
    return self._log_handlers

  def get_logger(self):
    # Get own logger
    return self._logger

  def get_session_log_files(self):
    log_files = []
    # Add the application logs
    log_files.extend(self._log_files)
    return log_files

  def get_logger_for(self, name):
    """
        Create a logger on demand
        :param name: name to be used in the logger
        :return: a new logger on that name
        """
    self.get_logger().debug("Creating logger with name {}".format(name))
    lg = logging.getLogger(name)
    for handler in self.get_log_handlers():
      lg.addHandler(handler)
    lg.setLevel(self._log_level)
    return lg

def read_json(json_file="json_file_not_specified.json"):
  """
    Reads a json file and it returns its object representation, no extra checks
    are performed on the file so, in case anything happens, the exception will
    reach the caller
    :param json_file: path to the file in json format to read
    :return: an object representation of the data in the json file
    """
  with open(json_file) as jf:
    return json.load(jf)

def read_yaml_from_file(yaml_file):
  """
  This function allows to read a yaml file with the configuration
  :param yaml_file: yaml file.
  :return: resturn the yaml data.
  """
  data = None
  with open(yaml_file, 'r') as f:
    data = yaml.safe_load(f.read())
  return data

def read_yaml_from_text(yaml_text):
  """
  Read the content of the yaml text into a data object
  :param yaml_text: yaml text
  :return:
  """
  return yaml.safe_load(yaml_text)

def check_create_folders(folders: List):
  """
    Check if folders exist, create them otherwise
    :param folders: list of folder paths to check
    :return: no return value
    """
  for folder in folders:
    if not os.path.exists(folder):
      try:
        os.makedirs(folder)
      except Exception as e:
        raise ToolBoxException(str(e))
    else:
      if not os.path.isdir(folder):
        raise ToolBoxException("'{}' is not a folder".format(folder))

def clear_cache():
  request.urlcleanup()

def download_file(file_url: str, file_name: str, log: logging, url_file = None) -> str:
  """
     Download file_url and move it to file_name, do nothing if file_name already exists.

    :param log: log to be use.
    :param file_url: file url to be download
    :param file_name: file name where the data will be downloaded
    :param url_file: the url file is used to write the urls to be downloaded, if None, the function will download the file
    :return: name of the file if the file can be download.
    """
  if os.path.isfile(file_name):
    return file_name

  if log is not None:
    log = logging

  if url_file is not None:
    url_file.write("{}\t{}\n".format(file_url, file_name))
    return file_name

  remaining_download_tries = REMAINING_DOWNLOAD_TRIES
  downloaded_file = None
  while remaining_download_tries > 0:
    try:
      downloaded_file, error_code = request.urlretrieve(file_url, file_name)
      log.debug("File downloaded -- " + downloaded_file)
      if downloaded_file.endswith('.gz'):
        extracted_file = downloaded_file.replace('.gz','')
        with open(extracted_file, 'wb') as outfile:
          with gzip.open(downloaded_file, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)
          os.remove(downloaded_file)
          downloaded_file = extracted_file
          log.debug("File extracted-- " + downloaded_file)
      break
    except (HTTPError, URLError, ContentTooShortError, ) as error:
      logging.error("Error downloading -- Incorrect URL or file not found: " + file_url + " on trial no: " + str(
        REMAINING_DOWNLOAD_TRIES - remaining_download_tries))
      log.error("Error code: " + str(error))
      remaining_download_tries = remaining_download_tries - 1
      downloaded_file = None
      continue
    except Exception as error:
      remaining_download_tries = remaining_download_tries - 1
      log.error("Error code: " + str(error))
      downloaded_file = None

  return downloaded_file

  # move the pep file to the desired name
  # if os.path.isfile(downloaded_file):
  #     if os.stat(downloaded_file).st_size > 1000:
  #         shutil.move(downloaded_file, file_name)
  #         logging.debug("File copy to filesystem -- " + downloaded_file)
  #         return file_name
  #     else:
  #         print("Corrupt File (size<1kb): ", file_url)
  #         os.remove(downloaded_file)
  #         return None
  # else:
  #     print("Failed to download the file: ", file_url)
  #     return None

def check_create_folders_overwrite(folders):
  """
    Given a list of folders, this method will create them, overwriting them in case they exist
    :param folders: list of folders to create
    :return: no return value
    :except: if any element in the list of folders is not a folder, an exception will be raised
    """
  invalid_folders = []
  for folder in folders:
    if os.path.exists(folder):
      if not os.path.isdir(folder):
        invalid_folders.append(folder)
  if invalid_folders:
    # If there's any invalid folder, we don't make any change, and we report the situation by raising an exception
    raise ToolBoxException("The following folders ARE NOT FOLDERS - '{}'"
                           .format(invalid_folders))
  for folder in folders:
    try:
      shutil.rmtree(folder)
    except FileNotFoundError as e:
      # It is find if the folder is not there
      pass
  check_create_folders(folders)

def create_latest_symlink(destination_path):
  """
    Create a symlink 'latest' to the given destination_path in its parent folder, i.e. if the given path is
    '/nfs/production/folder', the symlink will be
            /nfs/production/latest -> /nfs/production/folder
    :param destination_path: destination path where the symlink will point to
    :return: no return value
    """
  symlink_path = os.path.join(os.path.dirname(destination_path), 'latest')
  os.symlink(destination_path, symlink_path)

def create_latest_symlink_overwrite(destination_path):
  """
    Create a symlink 'latest' to the given destination_path in its parent folder, i.e. if the given path is
    '/nfs/production/folder', the symlink will be
            /nfs/production/latest -> /nfs/production/folder
    If there already is a 'latest' symlink, it will be overwritten
    :param destination_path: destination path where the symlink will point to
    :return: no return value
    """
  symlink_path = os.path.join(os.path.dirname(destination_path), 'latest')
  if os.path.islink(symlink_path):
    os.unlink(symlink_path)
  os.symlink(destination_path, symlink_path)

def gunzip_files(files):
  """
    Given a list of paths for Gzip compressed files, this method will uncompress them, returning a list with the files
    that could not be gunzipped and the reason why that happened
    :param files: list of paths to files that will be un-compressed
    :return: a list of possible failing to uncompress files
    """
  gunzip_command_template = "gunzip {}"
  files_with_error = []
  for file in files:
    if os.path.isfile(file):
      try:
        gunzip_subprocess = subprocess.Popen(gunzip_command_template.format(file),
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             shell=True)
        # Timeout, in seconds, is either 10 seconds or the size of the file in MB * 10, e.g. 1MB -> 10 seconds
        file_size_mb = os.path.getsize(file) / (1024 * 1024)
        timeout = max(10, int(file_size_mb))
        (stdout, stderr) = gunzip_subprocess.communicate(timeout=timeout)
        if gunzip_subprocess.poll() is not None:
          if gunzip_subprocess.returncode != 0:
            # ERROR - Report this
            err_msg = "ERROR uncompressing file '{}' output from subprocess STDOUT: {}\nSTDERR: {}" \
              .format(file, stdout.decode('utf8'), stderr.decode('utf8'))
            files_with_error.append((file, err_msg))
      except subprocess.TimeoutExpired as e:
        err_msg = "TIMEOUT ERROR uncompressing file '{}', size {}MB, given timeframe of '{}seconds', output from subprocess STDOUT: {}\nSTDERR: {}" \
          .format(file_size_mb,
                  timeout,
                  file,
                  stdout.decode('utf8'),
                  stderr.decode('utf8'))
        files_with_error.append((file, err_msg))
      except Exception as e:
        err_msg = "UNKNOWN ERROR uncompressing file '{}' ---> {}\nOutput from subprocess STDOUT: {}\nSTDERR: {}" \
          .format(file,
                  e,
                  stdout.decode('utf8'),
                  stderr.decode('utf8'))
        files_with_error.append((file, err_msg))
    else:
      files_with_error.append((file, "it IS NOT A FILE"))
  return files_with_error

def parse_peptide_groups(peptide_groups_prefix):
  peptide_groups = {}
  for group in peptide_groups_prefix.split(";"):
    lt = group.split(":")
    class_group = lt[0].replace("{","")
    classes = [x.replace("[","").replace("]","") for x in lt[1].split(",")]
    peptide_groups[class_group] = classes
  return peptide_groups

def parse_peptide_classes(peptide_classes_prefix):
  peptide_groups = {}
  for class_peptide in peptide_classes_prefix.split(","):
    peptide_groups[class_peptide] = [class_peptide]
  return peptide_groups

def is_peptide_group(peptide_group_members, accessions):
  """
  Given a group of classes and a list of accessions of a peptide. Returns True if all accessions match to exactly one class in the group.
  :param peptide_group_members: all protein classes
  :param accessions:  all protein accessions associated with the peptide.
  :return: True if all protein accessions belows to one of these peptide_group_members.
  """

  accession_group = 0
  for accession in accessions:
    for class_peptide in peptide_group_members:
      if class_peptide in accession:
        accession_group += 1
  return len(accessions) == accession_group

def is_peptide_decoy(accessions, prefix):
  return  any(prefix in s for s in accessions)

if __name__ == '__main__':
  print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
