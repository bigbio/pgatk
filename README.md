# Python tools for ProteoGenomics Analysis Toolkit


[![Build Status](https://travis-ci.org/bigbio/py-pgatk.svg?branch=master)](https://travis-ci.org/bigbio/py-pgatk)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pypgatk/README.html)


**pypgatk** is a Python library part of the [ProteoGenomics Analysis Toolkit](https://pgatk.readthedocs.io/en/latest). It provides different bioinformatics tools for proteogenomics data analysis.

# Requirements:

This package requirements vary depending on the way that you want to install it (all three are independent, you don't need all these requirements):

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

You can install SCCAF with bioconda (please setup conda and the bioconda channel if you haven't first, as explained [here](https://bioconda.github.io/user/index.html)):

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
docker run -it quay.io/biocontainers/pypgatk:0.0.3--py_0
```

Inside the container, you can either use the Python interactive shell or the command line version (see below).


## Use latest source code

Alternatively, for the latest version, clone this repo and go into its directory, then execute `pip3 install .`:

```
git clone https://github.com/bigbio/py-pgatk
cd py-pgatk
# you might want to create a virtualenv for pypgatk before installing
pip3 install .
```


In order to execute a task in `pypgatk_cli` the user should use a `COMMAND` that perform the specific task and the specific task arguments/options:

```bash
$: python3.7 pypgatk_cli.py -h
Usage: pypgatk_cli.py [OPTIONS] COMMAND [ARGS]...

  This is the main tool that give access to all commands and options provided by the pypgatk

Options:
  -h, --help  Show this message and exit.
```

### Please read the docs here: <https://pgatk.readthedocs.io/en/latest/pypgatk.html>
