.. _pypgatk


Pypgatk: Python Tools for ProteoGenomics
========================================

The Pypgatk framework and library provides a set of tools to perform proteogenomics analysis. 
In order to execute a task in ``pypgatk`` the user should use a ``COMMAND`` to perform the specific task and specify the
task arguments:


.. code-block:: bash
   :linenos:

   $: python pypgatk -h
      Usage: pypgatk_cli.py [OPTIONS] COMMAND [ARGS]...

      This is the main tool that gives access to all commands and options provided by the pypgatk_cli

      Options:
         -h, --help  Show this message and exit.

      Commands:
        ensembl-downloader       Command to download the ensembl information
        cbioportal-downloader    Command to download the the cbioportal studies
        cosmic-downloader        Command to download the cosmic mutation database
        dnaseq-to-proteindb      Command to translate sequences generated from RNA-seq and DNA sequences
        vcf-to-proteindb         Command to translate genomic variatns to protein sequences
        cbioportal-to-proteindb  Command to translate cbioportal mutation data into proteindb
        cosmic-to-proteindb      Command to translate Cosmic mutation data into proteindb
        generate-decoy      	 Command to generate decoy database from a proteindb


.. _installation:

Installation
------------

Clone the source code for pypgatk from source:

.. code-block:: bash
   
   git clone https://github.com/bigbio/py-pgatk.git
 

``pypgatk`` depends on several ``Python3`` packages that are listed in ``requirements.txt``, once in the downloaded directory install the dependencies using ``pip``:

.. code-block:: bash
   
   pip install -r requirements.txt


Install the ``pypgatk`` package from source:

.. code-block:: bash
   
   python setup.py install


.. _data-downloader:

Data Downloader Tools
------------------

The Data downloader is a set of COMMANDs to download data from different Genomics data providers including ENSEMBL, COSMIC and cBioPortal.

.. _ensembl-downloader:

Downloading ENSEMBL Data.
~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading data from `ENSEMBL <https://www.ensembl.org/info/data/ftp/index.html>`_ can be done using the command :ref:`ensembl_downloader <ensembl-downloader>`. 
The current tool enables downloading the following files for any taxonomy that is available ENSEMBL:

- GTF
- Protein Sequence (FASTA)
- CDS (FASTA)
- CDNA sequences (FASTA)
- Non-coding RNA sequences (FASTA)
- Nucleotide Variation (VCF)

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk_cli.py ensembl-downloader -h
      Usage: pypgatk_cli.py ensembl-downloader [OPTIONS]

      This tool enables to download from ENSEMBL ftp the FASTA, GTF and VCF files

      Options:
        -c, --config_file TEXT          Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT     Output directory for the peptide databases
        -fp, --folder_prefix_release TEXT Output folder prefix to download the data
        -t, --taxonomy TEXT             Taxonomy identifiers (comma separated) that will be use to download the data from Ensembl
        -l, --list_taxonomies TEXT             List the available species from Ensembl, users can find the desired taxonomy identifier from this list.
		-sv, --skip_vcf                 Skip the vcf file during the download
        -sg, --skip_gtf                 Skip the gtf file during the download
        -sp, --skip_protein             Skip the protein fasta file during download
        -sc, --skip_cds                 Skip the CDS file download
		-sd, --skip_cdna              	Skip the cDNA file download
        -sn, --skip_ncrna              Skip the ncRNA file download
        -h, --help                      Show this message and exit.


.. _ensembl-downloader_example:

**Examples**

- List all species without downloading any data::

	python pypgatk_cli.py ensembl-downloader -l -sv -sg -sp -sc -sd -sn

- Download all files except cDNA for Tureky (species id=9103, note that th species id cab be obtained from the list above):: 

	python pypgatk_cli.py ensembl-downloader -t 9103 -sd -o ensembl_files

- *[To be implemented]* Download CDS file for Humans (species id=9606) from release 94 and genome assembly GRCh37 :: 

	python pypgatk_cli.py ensembl-downloader -t 9606 -sv -sg -sp -sd -sn -o ensembl_files --release 94 --assembly GRCh37

.. note:: By default the command ``ensembl-downloader`` downloads all datasets for all species from the latest ENSEMBL release. To limit the download to a particular species specify the species identifier using the ``-t`` option. To list all available species run the command with ``-l (--list_taxonomies)`` option.

.. note:: Any of the file types can be skipped using the corresponding option. For example, to avoid downloading the protein sequence fasta file, use the argument ``--skip_protein``. Also, note that not all file types exists for all species so obviously the downloaded files depends on availabiliy of the dataset in ENSEMBL.

.. hint:: a VCF file per chromosome is downloaded for homo sapiens due to the large file size they have been distributed this way by ENSEMBL. For other species, a single VCF including all chromosomes is downloaded.  

.. _cosmic-downloader:


Downloading COSMIC Data.
~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading mutation data from `COSMIC <https://cancer.sanger.ac.uk/cosmic>`_ is performed using the COMMAND ``cosmic-downloader``. 
The current COMMAND allows users to download the following files:

- Cosmic mutation file (CosmicMutantExport)
- Cosmic all genes (All_COSMIC_Genes)

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk_cli.py cosmic-downloader -h
      Usage: pypgatk_cli.py cosmic-downloader [OPTIONS]

      Required parameters:
        -u, --username TEXT          Username for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
        -p, --password TEXT          Password for cosmic database -- please if you dont have one register here (https://cancer.sanger.ac.uk/cosmic/register)
	  
	  Optional parameters:
        -c, --config_file TEXT       Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT  Output directory for the peptide databases
        -h, --help                   Show this message and exit.
        
.. note:: In order to be able to download COSMIC data, the user should provide a user and password. Please first register in COSMIC database (https://cancer.sanger.ac.uk/cosmic/register).

.. _cosmic-downloader_example:

**Examples**

- Downlaod ``CosmicMutantExport.tsv.gz`` and ``All_COSMIC_Genes.fasta.gz``::
	
	python pypgatk_cli.py cosmic-downloader -u userName -p passWord -c config/cosmic_config.yaml -o cosmic_files

.. _cbioportal-downloader:


Downloading cBioPortal Data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Downloading mutation data from `cBioPortal <https://www.cbioportal.org/>`_ is performed using the command ``cbioportal-downloader``. 
cBioPortal stores mutation data from multiple studies (https://www.cbioportal.org/datasets). Each dataset in cBioPortal has an associated study_id.

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python3.7 pypgatk_cli.py cbioportal-downloader -h
      Usage: pypgatk_cli.py cbioportal-downloader [OPTIONS]

      Options:
        -c, --config_file TEXT Configuration file for the ensembl data downloader pipeline
        -o, --output_directory TEXT  Output directory for the peptide databases
        -l, --list_studies           Print the list of all the studies in cBioPortal (https://www.cbioportal.org)
        -d, --download_study TEXT    Download a specific Study from cBioPortal -- (all to download all studies)
        -h, --help                   Show this message and exit.


.. note:: 
	The argument ``-l`` (``--list_studies``) allows the user to list all the studies stored in cBioPortal. 
	The ``-d`` (``--download_study``) argument can be used to obtain mutation data from a particular study.

.. _cbioportal-downloader_example:

**Examples**

- Download data for study ID `blca_mskcc_solit_2014 <https://www.cbioportal.org/study/summary?id=blca_mskcc_solit_2014>`_::
	
	python pypgatk_cli.py cbioportal-downloader -d blca_mskcc_solit_2014 -o cbiportal_files
   
- Download data for all studies in cBioPortal::

	python pypgatk_cli.py cbioportal-downloader -d all -o cbioportal_files

If you face issues downloading all studies from cBioPortal using the ``cbioportal-downloader``, please download the studies from the `data hub <https://github.com/cBioPortal/datahub/tree/master/public>`_ through ``git-lfs`` 
which is used to download large files from gitHub repositories, see `installation instructions: <https://github.com/git-lfs/git-lfs/wiki/Installation>`_.

Following `instructions given on the datahub repositority <https://github.com/cBioPortal/datahub>`_, download the entire list of datasets using:: 
	
	git clone https://github.com/cBioPortal/datahub.git
	cd datahub
	git lfs install --local --skip-smudge
	git lfs pull -I public --include "data_clinical_sample.txt"
	git lfs pull -I public --include "data_mutations_mskcc.txt"
	
	
.. _generate-proteindb:


Generate Protein Databases
--------------------------

The **Pypgatk** framework provides a set of tools (COMMAND) to generate protein databaseas in ``FASTA`` format from DNA sequences, variants, and mutations. In order to perform this task, we have implemented multiple
commands depending on data type provided by the user and the public data providers (cBioPortal, COSMIC and ENSEMBL).

.. _cosmic-to-proteindb:

Cosmic Mutations to Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`COSMIC <https://cancer.sanger.ac.uk/cosmic/>`_ the Catalogue of **Human** Somatic Mutations in Cancer – is the world's largest source of expert manually curated somatic mutation information relating to human cancers. 
The command ``cosmic-to-proteindb`` converts the cosmic somatic mutations file into a protein sequence database file.

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk_cli.py cosmic-to-proteindb -h
      Usage: pypgatk_cli.py cosmic-to-proteindb [OPTIONS]

      Required parameters:
        -in, --input_mutation TEXT  Cosmic Mutation data file
        -fa, --input_genes TEXT     All Cosmic genes
        -out, --output_db TEXT      Protein database including all the mutations
      
      Optional parameters:
        -c, --config_file TEXT      Configuration file for the cosmic data pipelines
        -t, --tissue_type           Only consider mutations from these tissue tyoes, by default mutations from all tissue types are considered (default ``all``)
        -s,	--split_by_tissue_type  Generate a proteinDB output file for each tissue type in the mutations file (affected by ``--tissue_type``) (default ``False``)
        -h, --help                  Show this message and exit.

The file input of the tool ``-in`` (``--input_mutation``) is the cosmic mutation data file. 
The genes file ``-fa`` (``--input_genes``) contains the original CDS sequence for all genes used by the COSMIC team to annotate the mutations. 
:ref:`Use cosmic-downloader <cosmic-downloader_example>` to obtain the input files from COSMIC.

The output of the tool is a protein fasta file and is written in the following path `-out` (``--output_db``)

.. _cosmic-to-proteindb_example:

**Examples** 

- Generate cancer-type specific protein databases. For each cancer type in COSMIC generate a protein database based on the Primary site given in the mutations file::
  
   python pypgatk_cli.py cosmic-to-proteindb -in CosmicMutantExport.tsv -fa All_COSMIC_Genes.fasta -out cosmic_proteinDB.fa --split_by_tissue_type


.. _cbioportal-to-proteindb:

cBioPortal Mutations to Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. 
The available datasets can be viewed in this web page (https://www.cbioportal.org/datasets). 
The command ``cbioportal-to-proteindb`` converts the bcioportal mutations file into a protein sequence database file.

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk_cli.py cbioportal-to-proteindb -h
      Usage: pypgatk_cli.py cbioportal-to-proteindb [OPTIONS]

      Options:
        -c, --config_file TEXT           Configuration for cBioportal
        -in, --input_mutation TEXT       Cbioportal mutation file
        -fa, --input_cds TEXT            CDS genes from ENSEMBL database
        -out, --output_db TEXT           Protein database including the mutations
        -t, --tissue_type TEXT           Only consider mutations from these tissue tyoes, by default mutations from all tissue types are considered (default ``all``)
        -s,	--split_by_tissue_type BOOL  Generate a proteinDB output file for each tissue type in the mutations file (affected by ``--tissue_type``) (default ``False``)
        -c, --clinical_sample_file TEXT  Clinical sample file that contains the cancery type per sample identifier 
        -h, --help                       Show this message and exit.

.. note:: The clinical sample file for each mutation file can be found under the same directory as the mutation file downloaded from cBioportal (It should have at least two columns named: Cancer Type and Sample Identifier). The file is only needed if generating tissue type databases is desired (that is when -s or -t is given).

The file input of the tool ``-in`` (``--input_mutation``) is the cbioportal mutation data file. 
An example is given in :ref:`cbioportal-downloader <cbioportal-downloader_example>` showing how to obtain the mutations file for a particular study.
The CDS sequence for all genes input file ``-fa`` (``--input_genes``) can be obtained using the ENSEMBL CDS files, see :ref:`this example <ensembl-downloader_example>`. 
The output of the tool is a protein fasta file and it is written in the following path ``-out`` (``--output_db``)

.. note:: The cBioportal mutations are aligned to the hg19 assembly, make sure that the correct genome assembly is selected for the download.

.. _cbioportal-to-proteindb_example:

**Examples**

- translate mutations from ``Bladder`` samples in studyID: ``blca_mskcc_solit_2014`` (:ref:`use cbioportal-downloader <cbioportal-downloader_example>` to download the study, then extract the content of the downloaded file)::
	
	python pypgatk_cli.py cbioportal-to-proteindb --config_file config/cbioportal_config.yaml --input_cds human_hg19_cds.fa  --input_mutation data_mutations_mskcc.txt --clinical_sample_file data_clinical_sample.txt --output_db bladder_proteindb.fa

.. _vcf-to-proteindb:


Annotated Variants (VCF) to Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Variant Calling Format (VCFv4.1) is a text file representing genomic variants. 
Variant calling methods generate a VCF file that can be used as input to VEP for variant annotation. 
VEP reports the trasncripts that are affected by each variant along with the consequences of the effect. 

The ``vcf_to_proteindb`` COMMAND takes a VEP-annotated VCF and translates the genomic variants in the VCF that affect protein-coding transcripts. It also allows for other variants to be translated by selecting the desired biotypes and consequences.

Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk_cli.py vcf-to-proteindb -h
      Usage: pypgatk_cli.py vcf-to-proteindb [OPTIONS]

      Required parameters:
        -c, --config_file TEXT      Configuration for VCF conversion parameters
        --vep_annotated_vcf         VCF file containing the annotated genomic variants
        --gene_annotations_gtf        Gene models in the GTF format that is used with VEP
        --input_fasta         Fasta sequences for the transripts in the GTF file used to annotated the VCF
        --output_proteindb          Output file to write the resulting variant protein sequences
      
      Options:
        --translation_table INTEGER     Translation table (Default 1). Please see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi> for identifiers of translation tables.
        --mito_translation_table INTEGER	Mito_trans_table (default 2), also from <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi> 
        --var_prefix TEXT 	String to add before the variant peptides
        --report_ref_seq	In addition to variant peptides, also report the reference peptide from the transcript overlapping the variant 
        --output_proteindb TEXT	Output file name, exits if already exists
        --annotation_field_name TEXT	Annotation Field name found in the INFO column, e.g CSQ or vep
      	--af_field TEXT	Field name in the VCF INFO column that shows the variant allele frequency (VAF, default is AF).
      	--af_threshold FLOAT      Minium allele frequency threshold for considering the variants
  		--transcript_index INTEGER	Index of transcript ID in the annotated columns in the VCF INFO field (separated by |) (default is 3)
 		--consequence_index INTEGER	Index of consequence in the annotated columns in the VCF INFO field (separated by |) (default is 1)
 		--exclude_biotypes TEXT         Variants affecting gene/transcripts in these biotypes will not be considered for translation (affected by include_biotypes). 
  		--exclude_consequences TEXT     Variants with these consequences will not be considered for translation (default: downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant)
        --skip_including_all_cds	By default any affected transcript that has a defined CDS will be translated, this option disables this features instead it only depends on the specified biotypes
  		--include_biotypes TEXT	Translate affected transcripts that have one of these biotypes
  		--include_consequences TEXT	Consider variants that have one of these consequences (default is all) (for the list of consequences see: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html.
  		--biotype_str TEXT	String used to identify gene/transcript biotype in the gtf file (default transcript_biotype).
  		--ignore_filters	Enabling this option causes all variants to be parsed. By default only variants that have not failed any filters will be processed (FILTER field is PASS, None, .) or if the filters are subset of the accepted_filters (default is False)
  		--accepted_filters TEXT	Accepted filters for variant parsing
        -h, --helP		Show this message and exit.

The file input of the tool ``--vcf_annotated_vcf`` is a VCF file that can be provided by the user or obtained from ENSEMBL using :ref:`ensembl_downloader <ensembl-downloader>`, see :ref:`an example here <ensembl-downloader_example>`. 
The ``gene_annotations_gtf`` file can also be obtained with the :ref:`ensembl_downloader <ensembl-downloader>`. 
The GTF file should match the one used for the variant annotation in VEP. 
The ``--input_fasta`` file contains the ``CDS`` and DNA sequences for all genes present in the GTF file. 
This file can be generated from the GTF file using the ``gffread`` tool as follows::
	
	$: gffread -F -w input_fasta.fa -g genome.fa gene_annotations_gtf

The output of the tool is a protein fasta file and is written in the following path ``--output_proteindb``.


.. _vcf-to-proteindb_examples:

**Examples**

- Translate human *missense* variants from ENSEMBL VCFs that have a minimum *AF 5%* and affect any *protein_coding* gene or *lincRNAs*::
	
	python pypgatk.py vcf-to-proteindb 
 		--vep_annotated_vcf homo_sapiens_incl_consequences.vcf 
 		--include_biotypes lncRNA 
 		--include_consequences missense 
 		--af_threshold 0.05

.. note:: 
	- By default  vcf-to-proteindb considers transcript that have a coding sequence that includes all protein_coding genes. In order to also include lincRNAs we use the ``--include_biotypes`` option that accepts multiple entries separated by comma. The biotypes can be any of the ENSEMBL gene/transcript biotypes: https://www.ensembl.org/info/genome/genebuild/biotypes.html. 
	- The choice of using gene or transcript biotype can be specified using the ``--biotype_str option``. 
	- Also, by default all consequences are accepted except those given with ``--exclude_biotypes``. See the list consequences of consequences generated by VEP: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

- Translate human *missense* variants or *inframe_insertion* from gnoMAD VCFs that have a minmum 1% allele frquency in control samples and affect any protein_coding gene::
	
	$: python pypgatk.py vcf-to-proteindb 
 		--vep_annotated_vcf gnmad_genome.vcf 
 		--include_consequences missense, frameshift_insert 
 		--annotation_field_name vep --af_threshold 0.01 
 		--af_field control_af 
 		--biotype_str transcript_type 
 		--transcript_index 6

.. hint:: 
	- By default  ``vcf-to-proteindb`` considers transcript that have a coding sequence which includes all *protein_coding* transcripts and since the required biotype is protein coding transcripts thereore there is no need to specify any biotypes.  
	- The provided VCF file has some specific properties: the annotation field is specified with the string *vep* hence the ``--annotation_field_name parameter``,  the transcriptat the sixth position in the annotation field, and since gnomAD collects variants from many sources it provides allele frequencies across many many sub-populations and sub-groups, in this case the goal is to use only variants that are common within control samples therefroe the ``--af_field`` is set to ``control_af``. 
	- Since gnomAD uses GENCODE gene annotations for annotation the variants we need to change the default ``biotype_str`` from *transcript_biotype* to *transcript_type* (as written in the GTF file).

.. note:: 
	As shown in the two examples above, when ENSEMBL data is used, the default options should work. 
	However, for using other data sources such as variants from gnomAD, GTF from GENOCODE and others one or more of the following parameters need to be changed:
	
		--af_field (from the VCF INFO field)
			
		--annotation_field_name (from the VCF INFO field)
			
		--transcript_index (from the annotation field in the VCF INFO field)
			
		--consequence_index (from the annotation field in the VCF INFO field)
			
		--biotype_str (from the GTF INFO field)


.. _dnaseq-to-proteindb:

Transcripts (DNA) to Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DNA sequences given in a fasta format can be translated using the ``dnaseq-to-proteindb`` tool. This tool allows for translation 
of all kinds of transcripts (coding and noncoding) by specifying the desired biotypes.
The most suited ``--input_fasta`` file can be generated from a given GTF file using the ``gffread`` commad as follows::
	
	$: gffread -F -w transcript_sequences.fa -g genome.fa gene_annotations_gtf

The fasta file that is generated from the GTF file would contain DNA sequences for all transcripts regardless of their biotypes. Also, it specifies the CDS positions for the protein coding transcripts.
The ``dnaseq-to-proteindb`` command recognizes the features such as biotype and expression values in the fasta header that are taken from the GTF INFO filed (if available).
However, it is not required to have those information in the fasta header but their presence enables the user to filter by biotype and expression values during the translation step. 


Command Options
^^^^^^^^^^^^^^

.. code-block:: bash
   :linenos:

   $: python pypgatk.py dnaseq-to-proteindb -h
      Usage: pypgatk.py dnaseq-to-proteindb [OPTIONS]

      Required parameters:
        -c, --config_file TEXT      Configuration for VCF conversion parameters
        --input_fasta         Fasta sequences for the transripts in the GTF file used to annotated the VCF
        --output_proteindb          Output file to write the resulting variant protein sequences
        
      Optional parameters:	
  		--translation_table INTEGER    Translation Table (default 1)
  		--num_orfs INTEGER             Number of ORFs (default 0)
  		--num_orfs_complement INTEGER  Number of ORFs from the reverse side (default 0)
  		--skip_including_all_cds       By default any transcript that has a defined CDS will be translated, this option disables this features instead it only depends on the biotypes
  		--include_biotypes TEXT        Translate sequences with the spcified biotypes. Multiple biotypes can be given separated by comma. To translate all sequences in the input_fasta file set this option to ``all`` (default None).
  		--exclude_biotypes TEXT        Skip sequences with unwanted biotypes (affected by --include_biotypes) (default None). 
  		--biotype_str TEXT             String used to identify gene/transcript biotype in the fasta file (default transcript_biotype).
  		--expression_str TEXT          String to be used for extracting expression value (TPM, FPKM, etc) (default None).
  		--expression_thresh FLOAT      Threshold used to filter transcripts based on their expression values (default 5, affected by --expression_str) 
  		-h, --help                     Show this message and exit

.. _dnaseq-to-proteindb_examples:

**Examples**

- Generate the canonical protein database, i.e. translate all *protein_coding* transcripts::
	
	python pypgatk.py dnaseq-to-proteindb 
		--config_file config/ensembl_config.yaml 
		--input_fasta testdata/transcript_sequences.fa 
		--output_proteindb testdata/proteindb_from_CDSs_DNAseq.fa

- Generate a protein database from lincRNA and canonical proteins::

	python pypgatk.py dnaseq-to-proteindb 
		--config_file config/ensembl_config.yaml 
		--input_fasta testdata/transcript_sequences.fa 
		--output_proteindb testdata/proteindb_from_processed_pseudogene.fa
		--include_biotypes lincRNA
	
- Generate a protein database from processed pseudogene::

	python pypgatk.py dnaseq-to-proteindb 
		--config_file config/ensembl_config.yaml 
		--input_fasta testdata/transcript_sequences.fa 
		--output_proteindb testdata/proteindb_from_processed_pseudogene.fa
		--include_biotypes processed_pseudogene,transcribed_processed_pseudogene,translated_processed_pseudogene
		--skip_including_all_cds 
	
- Generate alternative ORFs from canonical sequences::	
	
	python pypgatk.py dnaseq-to-proteindb 
		--config_file config/ensembl_config.yaml 
		--input_fasta testdata/transcript_sequences.fa 
		--output_proteindb testdata/proteindb_from_processed_pseudogene.fa
		--include_biotypes altORFs
		--skip_including_all_cds


Contributions
-------------

- Husen M. Umer ([husensofteng](https://github.com/husensofteng))
- Yafeng Zhu ([yafeng](http://github.com/yafeng))
- Enrique Audain ([enriquea](https://github.com/enriquea))
- Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
