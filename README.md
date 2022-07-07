# ProteoGenomics Analysis Toolkit


![Python application](https://github.com/bigbio/py-pgatk/workflows/Python%20application/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pypgatk/README.html)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f6d030fd7d69413987f7265a01193324)](https://www.codacy.com/gh/bigbio/py-pgatk/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bigbio/py-pgatk&amp;utm_campaign=Badge_Grade)
[![PyPI version](https://badge.fury.io/py/pypgatk.svg)](https://badge.fury.io/py/pypgatk)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pypgatk)

**pypgatk** is a Python library - part of the [ProteoGenomics Analysis Toolkit](https://pgatk.readthedocs.io/en/latest). It provides different bioinformatics tools for proteogenomics data analysis.

# Requirements:

The package requirements vary depending on the way that you want to install it (you need one of the following):

- pip: if installation goes through pip, you will require Python3 and pip3 installed.
- Bioconda: if installation goes through Bioconda, you will require that [conda is installed and configured to use bioconda channels](https://bioconda.github.io/user/index.html).
- Docker container: to use pypgatk from its docker container you will need [Docker](https://docs.docker.com/install/) installed.
- Source code: to use and install from the source code directly, you will need to have git, Python3 and pip.

# Installation

## pip

You can install pypgatk with pip:

```
pip install pypgatk
```

## Bioconda

You can install pypgatk with bioconda (please setup conda and the bioconda channel if you haven't first, as explained [here](https://bioconda.github.io/user/index.html)):

```
conda install pypgatk
```

## Available as a container

You can use the pypgatk tool already setup on a Docker container. You need to choose from the available tags [here](https://quay.io/repository/biocontainers/pypgatk?tab=tags) and replace it in the call below where it says `<tag>`.

```
docker pull quay.io/biocontainers/pypgatk:<tag>
```

**NOTE**: Please note that Biocontainers containers do not have a latest tag, as such a docker pull/run without defining the tag will fail. For instance, a valid call would be (for version 0.0.2):

```
docker run -it quay.io/biocontainers/pypgatk:0.0.2--py_0
```

Inside the container, you can either use the Python interactive shell or the command line version (see below).


## Use latest source code

Alternatively, for the latest version, clone this repo and go into its directory, then execute `pip3 install .` :

```
git clone https://github.com/bigbio/py-pgatk
cd py-pgatk
# you might want to create a virtualenv for pypgatk before installing
pip3 install .
```

# Usage

The pypgatk design combines multiple modules and tools into one framework. All the possible commands are accessible using the commandline tool `pypgatk_cli.py`.

The library provides multiple commands to download, translate and generate protein sequence databases from reference and mutation genome databases.

```
$: pypgatk_cli -h

Usage: pypgatk [OPTIONS] COMMAND [ARGS]...

  This is the main tool that give access to all commands and options
  provided by the pypgatk

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  cbioportal-downloader    Command to download the the cbioportal studies
  cbioportal-to-proteindb  Command to translate cbioportal mutation data into
                           proteindb
  cosmic-downloader        Command to download the cosmic mutation database
  cosmic-to-proteindb      Command to translate Cosmic mutation data into
                           proteindb
  dnaseq-to-proteindb      Generate peptides based on DNA sequences
  ensembl-check            Command to check ensembl database for stop codons,
                           gaps
  ensembl-downloader       Command to download the ensembl information
  generate-decoy           Create decoy protein sequences using multiple
                           methods DecoyPYrat, Reverse/Shuffled Proteins.
  generate-deeplc          Generate input for deepLC tool from idXML,mzTab or
                           consensusXML
  msrescore-configuration  Command to generate the msrescore configuration
                           file from idXML
  peptide-class-fdr        Command to compute the Peptide class FDR
  threeframe-translation   Command to perform 3'frame translation
  vcf-to-proteindb         Generate peptides based on DNA variants VCF files

```

# Full Documentation

[https://pgatk.readthedocs.io/en/latest/pypgatk.html](https://pgatk.readthedocs.io/en/latest/pypgatk.html)

## Cite as
Husen M Umer, Enrique Audain, Yafeng Zhu, Julianus Pfeuffer, Timo Sachsenberg, Janne Lehtiö, Rui M Branca, Yasset Perez-Riverol
Generation of ENSEMBL-based proteogenomics databases boosts the identification of non-canonical peptides
Bioinformatics, Volume 38, Issue 5, 1 March 2022, Pages 1470–1472
https://doi.org/10.1093/bioinformatics/btab838

