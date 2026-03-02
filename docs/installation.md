# Installation

## Requirements

The package requirements vary depending on the way that you want to install it (you need one of the following):

- **pip**: if installation goes through pip, you will require Python3 and pip3 installed.
- **Bioconda**: if installation goes through Bioconda, you will require that [conda is installed and configured to use bioconda channels](https://bioconda.github.io/user/index.html).
- **Docker container**: to use pgatk from its docker container you will need [Docker](https://docs.docker.com/install/) installed.
- **Source code**: to use and install from the source code directly, you will need to have git, Python3 and pip.

## pip

You can install pgatk with pip:

```bash
pip install pgatk
```

## Bioconda

You can install pgatk with bioconda (please setup conda and the bioconda channel if you haven't first, as explained [here](https://bioconda.github.io/user/index.html)):

```bash
conda install pgatk
```

## Available as a container

You can use the pgatk tool already setup on a Docker container. You need to choose from the available tags [here](https://quay.io/repository/biocontainers/pgatk?tab=tags) and replace it in the call below where it says `<tag>`.

```bash
docker pull quay.io/biocontainers/pgatk:<tag>
```

!!! note
    Biocontainers containers do not have a `latest` tag, so a docker pull/run without defining the tag will fail. For instance, a valid call would be (for version 0.0.2):

    ```bash
    docker run -it quay.io/biocontainers/pgatk:0.0.2--py_0
    ```

Inside the container, you can either use the Python interactive shell or the command line version.

## Use latest source code

Alternatively, for the latest version, clone this repo and go into its directory, then execute `pip3 install .`:

```bash
git clone https://github.com/bigbio/pgatk
cd pgatk
# you might want to create a virtualenv for pgatk before installing
pip3 install .
```
