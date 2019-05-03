.. _pypgatk


Pypgatk: Python Tools for ProteoGenomics
===========================

The Pypgatk framework and library provides a set of tools and functionalities to perform proteogenomics analysis. In order to execute a task in `pypgatk` the user should use a `COMMAND` that perform the specific task and the
specific task arguments/options:

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk -h
      Usage: pypgatk.py [OPTIONS] COMMAND [ARGS]...

      This is the main tool that give access to all commands and options provided by the pypgatk

      Options:
         -h, --help  Show this message and exit.

      Commands:
        cbioportal-downloader    Command to download the the cbioportal studies
        cbioportal-to-proteindb  Command to translate cbioportal mutation data into proteindb
        cosmic-downloader        Command to download the cosmic mutation database
        cosmic-to-proteindb      Command to translate Cosmic mutation data into proteindb
        ensembl-downloader       Command to download the ensembl information


Data downloader Tool
------------------

The `Data downloader` is a set of `COMMANDs` to download data from different Genomics data providers such as ENSEMBL, COSMIC or cBioPortal.

Downloading ENSEMBL data.
~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading data from `ENSEMBL <https://www.ensembl.org/info/data/ftp/index.html>`_ can be done using the command `ensembl-downloader`. The current tool enables to download the following files for each taxonomy:

- GTF
- Protein Sequence (FASTA),
- CDS (FASTA)
- Variation (VCF))

.. hint:: By default the command `ensembl-downloader` download all file types for all the ENSEMBL species.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py ensembl-downloader -h
      Usage: pypgatk.py ensembl-downloader [OPTIONS]

      This tool enables to download from enseml ftp the FASTA and GTF files

      Options:
        -c, --config-file TEXT          Configuration file for the ensembl data downloader pipeline
        -o, --output-directory TEXT     Output directory for the peptide databases
        -fp, --folder-prefix-release TEXT Output folder prefix to download the data
        -t, --taxonomy TEXT             Taxonomy List (comma separated) that will be use to download the data from Ensembl
        -sg, --skip-gtf                 Skip the gtf file during the download
        -sp, --skip-protein             Skip the protein fasta file during download
        -sc, --skip-cds                 Skip the CDS file download
        -snr, --skip-ncrna              Skip the ncRNA file download
        -h, --help                      Show this message and exit.


Each of the file types can be skip using the corresponding option. For example, if the user do not want to download the the protein sequence fasta file, it can be skip by using the argument `pypgatk.py ensembl-downloader --skip-protein`

Downloading COSMIC data.
~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading mutation data from `COSMIC <https://cancer.sanger.ac.uk/cosmic>`_ is performed using the COMMAND `cosmic-downloader`. The current COMMAND allows users to download the following files:

- Cosmic mutation file (CosmicMutantExport)
- Cosmic all genes (All_COSMIC_Genes)

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cosmic-downloader -h
      Usage: pypgatk.py cosmic-downloader [OPTIONS]

      Options:
        -c, --config-file TEXT       Configuration file for the ensembl data downloader pipeline
        -o, --output-directory TEXT  Output directory for the peptide databases
        -u, --username TEXT          Username for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
        -p, --password TEXT          Password for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
        -h, --help                   Show this message and exit.

.. note:: In order to be able to download COSMIC data, the user should provide a user and password. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register).__

Downloading cBioPortal data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading mutation data from `cBioPortal <https://www.cbioportal.org/>`_ is performed using the command `cbioportal-downloader`. cBioPortal store multiple studies (https://www.cbioportal.org/datasets) containing mutation data.
Currently is not possible to search the studies by PubMedID, only can be search by study_id.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cbioportal-downloader -h
      Usage: pypgatk.py cbioportal-downloader [OPTIONS]

      Options:
        -c, --config-file TEXT Configuration file for the ensembl data downloader pipeline
        -o, --output-directory TEXT  Output directory for the peptide databases
        -l, --list-studies           Print the list of all the studies in cBioPortal (https://www.cbioportal.org)
        -d, --download-study TEXT    Download an specific Study from cBioPortal -- (all to download all studies)
        -h, --help                   Show this message and exit.


The argument `-l` (`--list-studies`) allow the users to list all the studies store in cBioPortal. If the user is interested in only one study, it can use the argument `-d` (`--download-study`).

From Genome information to protein sequence databases
----------------------------

The **Pypgatk** framework provides a set of tools (COMMAND) to convert genome mutation and variant databases to protein sequence databases (FASTA). In order to perform this task, we have implemented multiple
commands depending on the mutation provider (cBioPortal or COSMIC).

Cosmic Mutations to Proitein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`COSMIC <https://cancer.sanger.ac.uk/cosmic/>`_ the Catalogue of **Human** Somatic Mutations in Cancer – is the world's largest source of expert manually curated somatic mutation information relating to human cancers. The current tool use the command `cosmic-to-proteindb` to convert the cosmic somatic mutations file into a protein sequence database file.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cosmic-to-proteindb -h
      Usage: pypgatk.py cosmic-to-proteindb [OPTIONS]

      Options:
        -c, --config-file TEXT      Configuration file for the cosmic data pipelines
        -in, --input-mutation TEXT  Cosmic Mutation data file
        -fa, --input-genes TEXT     All Cosmic genes
        -out, --output-db TEXT      Protein database including all the mutations
        -h, --help                  Show this message and exit.

The file input of the tool `-in` (`--input-mutation`) is the cosmic mutation data file. The genes file `-fa` (`--input-genes`) contains the original CDS sequence for all genes used by the COSMIC team to annotate the mutations.
The output of the tool is a protein fasta file and will be written in the following path `-out` (--output-db)

cBioPortal Mutations to Protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. All datasets can be viewed in this web page (https://www.cbioportal.org/datasets). The current tool
use the command `cbioportal-to-proteindb` to convert the bcioportal mutations file into a protein sequence database file.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cbioportal-to-proteindb -h
      Usage: pypgatk.py cbioportal-to-proteindb [OPTIONS]

      Options:
        -c, --config-file TEXT      Configuration for
        -in, --input-mutation TEXT  Cbioportal mutation file
        -fa, --input-cds TEXT       CDS genes from ENSEMBL database
        -out, --output-db TEXT      Protein database including all the mutations
        -h, --help                  Show this message and exit.

The file input of the tool `-in` (`--input-mutation`) is the cbioportal mutation data file. The CDS sequence for all genes input file `-fa` (`--input-genes`) can be provided using the ENSEMBL CDS files. In order to download the CDS files, the user can use the `ensembl-downloader` command.
The output of the tool is a protein fasta file and will be written in the following path `-out` (--output-db)

Contributions
-----------------------

- Yafeng Zhu ([yafeng](http://github.com/yafeng))
- Husen M. Umer ([husensofteng](https://github.com/husensofteng))
- Enrique Audain ([enriquea](https://github.com/enriquea))
- Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
