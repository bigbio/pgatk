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
        vcf-to-proteindb         Command to translate genomic variatns to protein sequences
        dnaseq-to-proteindb      Command to translate sequences generated from RNA-seq and DNA sequences


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

      This tool enables to download from ENSEMBL ftp the FASTA, GTF and VCF files

      Options:
        -c, --config_file TEXT          Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT     Output directory for the peptide databases
        -fp, --folder_prefix_release TEXT Output folder prefix to download the data
        -t, --taxonomy TEXT             Taxonomy List (comma separated) that will be use to download the data from Ensembl
        -sv, --skip_vcf                 Skip the vcf file during the download
        -sg, --skip_gtf                 Skip the gtf file during the download
        -sp, --skip_protein             Skip the protein fasta file during download
        -sc, --skip_cds                 Skip the CDS file download
        -snr, --skip_ncrna              Skip the ncRNA file download
        -h, --help                      Show this message and exit.


Each of the file types can be skip using the corresponding option. For example, to avoid downloading the protein sequence fasta file, use the argument `pypgatk.py ensembl-downloader --skip_protein`

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
        -c, --config_file TEXT       Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT  Output directory for the peptide databases
        -u, --username TEXT          Username for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
        -p, --password TEXT          Password for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
        -h, --help                   Show this message and exit.

.. note:: In order to be able to download COSMIC data, the user should provide a user and password. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register).__

Downloading cBioPortal data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading mutation data from `cBioPortal <https://www.cbioportal.org/>`_ is performed using the command `cbioportal-downloader`. cBioPortal stores mutation data from multiple studies (https://www.cbioportal.org/datasets).
Currently, it is not possible to search the studies by PubMedID, they can only be search by study_id.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cbioportal-downloader -h
      Usage: pypgatk.py cbioportal-downloader [OPTIONS]

      Options:
        -c, --config_file TEXT Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT  Output directory for the peptide databases
        -l, --list_studies           Print the list of all the studies in cBioPortal (https://www.cbioportal.org)
        -d, --download_study TEXT    Download an specific Study from cBioPortal -- (all to download all studies)
        -h, --help                   Show this message and exit.


The argument `-l` (`--list_studies`) allow the users to list all the studies stored in cBioPortal. The `-d` (`--download_study`) argument can be used to obtain mutation data from a particular study.

From Genome information to protein sequence databases
----------------------------

The **Pypgatk** framework provides a set of tools (COMMAND) to convert genome mutation and variant databases to protein sequence databases (FASTA). In order to perform this task, we have implemented multiple
commands depending on the mutation provider (cBioPortal or COSMIC, ENSEMBL).

Cosmic Mutations to Protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`COSMIC <https://cancer.sanger.ac.uk/cosmic/>`_ the Catalogue of **Human** Somatic Mutations in Cancer – is the world's largest source of expert manually curated somatic mutation information relating to human cancers. The current tool uses the command `cosmic-to-proteindb` to convert the cosmic somatic mutations file into a protein sequence database file.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cosmic-to-proteindb -h
      Usage: pypgatk.py cosmic-to-proteindb [OPTIONS]

      Options:
        -c, --config_file TEXT      Configuration file for the cosmic data pipelines
        -in, --input_mutation TEXT  Cosmic Mutation data file
        -fa, --input_genes TEXT     All Cosmic genes
        -out, --output_db TEXT      Protein database including all the mutations
        -h, --help                  Show this message and exit.

The file input of the tool `-in` (`--input_mutation`) is the cosmic mutation data file. The genes file `-fa` (`--input_genes`) contains the original CDS sequence for all genes used by the COSMIC team to annotate the mutations.
The output of the tool is a protein fasta file and is written in the following path `-out` (--output-db)

cBioPortal Mutations to Protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. The available datasets can be viewed in this web page (https://www.cbioportal.org/datasets). The current tool
uses the command `cbioportal-to-proteindb` to convert the bcioportal mutations file into a protein sequence database file.

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py cbioportal-to-proteindb -h
      Usage: pypgatk.py cbioportal-to-proteindb [OPTIONS]

      Options:
        -c, --config_file TEXT      Configuration for cBioportal
        -in, --input_mutation TEXT  Cbioportal mutation file
        -fa, --input_cds TEXT       CDS genes from ENSEMBL database
        -out, --output_db TEXT      Protein database including the mutations
        -h, --help                  Show this message and exit.

The file input of the tool `-in` (`--input_mutation`) is the cbioportal mutation data file. The CDS sequence for all genes input file `-fa` (`--input_genes`) can be provided using the ENSEMBL CDS files. In order to download the CDS files, the user can use the `ensembl-downloader` command. Please note that the cBioportal mutations are aligned to the hg19 assembly, make sure that the correct genome assembly is selected for the download.
The output of the tool is a protein fasta file and it is written in the following path `-out` (--output_db)

Annotated variants (VCF) to protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Variant Calling Format (VCFv4.1) is a text file to represent genomic variants. Variant calling methods generate a VCF file that could be used as input with VEP for variant annotation. VEP reports the trasncripts that are affected by the each variant along with the consequences of the effect. The vcf_to_proteindb COMMAND takes the VEP-annotated VCF and translates all genomic variatns in the VCF that affect protein-coding transcripts. It also allows for other variants to be translated by selected the desired biotypes and consequences. 

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk.py vcf-to-proteindb -h
      Usage: pypgatk.py vcf-to-proteindb [OPTIONS]

      Options:
        -c, --config_file TEXT      Configuration for VCF conversion parameters
        --vep_annotated_vcf         VCF file containing the annotated genomic variants
        --gene_annotations_gtf        Gene models in the GTF format that is used with VEP
        --transcripts_fasta         Fasta sequences for the transripts in the GTF file used to annotated the VCF
        --output_proteindb          Output file to write the resulting variant protein sequences
        -h, --help                  Show this message and exit.

The file input of the tool `--vep_annotated_vcf` is the VCF file that can be obtained with the data_downloader COMMAND, for instance. The `gene_annotations_gtf` file can be obtained with the data_downloader COMMAND, for instance. The GTF file should match the one used for the variant annotation in VEP. The `--transcripts_fasta` file contains the CDS and DNA sequences for all genes present in the GTF file. This file can be generated from the GTF file using the gffread tool.

.. code-block:: bash
   :linenos:
   
   $: gffread -F -w transcripts_fasta.fa -g genome.fa gene_annotations_gtf
   
The output of the tool is a protein fasta file and is written in the following path `--output_proteindb`.


Transcripts (DNA) to Protein sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Contributions
-----------------------

- Yafeng Zhu ([yafeng](http://github.com/yafeng))
- Husen M. Umer ([husensofteng](https://github.com/husensofteng))
- Enrique Audain ([enriquea](https://github.com/enriquea))
- Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
