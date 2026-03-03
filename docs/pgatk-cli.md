# pgatk: Python Tools for ProteoGenomics

The pgatk framework and library provides a set of tools to perform proteogenomics analysis. In order to execute a task in `pgatk` the user should use a `COMMAND` to perform the specific task and specify the task arguments:

```bash
$ pgatk -h
Usage: pgatk [OPTIONS] COMMAND [ARGS]...

  This is the main tool that gives access to all commands and options provided by the pgatk

Options:
   --version   Show the version and exit.
   -h, --help  Show this message and exit.

Commands:
  cbioportal-downloader    Command to download the the cbioportal studies
  cbioportal-to-proteindb  Command to translate cbioportal mutation data into proteindb
  clinvar-to-proteindb     Generate protein database from ClinVar VCF + RefSeq GTF
  cosmic-downloader        Command to download the cosmic mutation database
  cosmic-to-proteindb      Command to translate Cosmic mutation data into proteindb
  digest-mutant-protein    Digest mutant proteins and filter against canonical proteome
  dnaseq-to-proteindb      Generate peptides based on DNA sequences
  ensembl-check            Command to check ensembl database for stop codons, gaps
  ensembl-downloader       Command to download the ensembl information
  generate-decoy           Create decoy protein sequences using multiple methods
  map-peptide2genome       Map peptides to genomic coordinates (GFF3 output)
  ncbi-downloader          Download NCBI RefSeq and ClinVar reference files
  threeframe-translation   Command to perform 3'frame translation
  vcf-to-proteindb         Generate peptides based on DNA variants VCF files
```

## Data Downloader Tools

The Data downloader is a set of commands to download data from different Genomics data providers including ENSEMBL, COSMIC and cBioPortal.

### Downloading ENSEMBL Data

Downloading data from [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html) can be done using the command `ensembl-downloader`. The current tool enables downloading the following files for any taxonomy that is available in ENSEMBL:

- GTF
- Protein Sequence (FASTA)
- CDS (FASTA)
- CDNA sequences (FASTA)
- Non-coding RNA sequences (FASTA)
- Nucleotide Variation (VCF)
- Genome assembly DNA sequences (FASTA)

#### Command Options

```bash
$ pgatk ensembl-downloader -h
Usage: pgatk ensembl-downloader [OPTIONS]

  This tool enables to download from ENSEMBL ftp the FASTA, GTF and VCF files

  Required parameters:
    -c, --config_file TEXT          Configuration file for the ensembl data downloader pipeline
    -o, --output_directory TEXT     Output directory for the peptide databases

  Optional parameters:
    -l, --list_taxonomies TEXT      List the available species from Ensembl
    -fp, --folder_prefix_release TEXT  Output folder prefix to download the data
    -t, --taxonomy TEXT             Taxonomy identifiers (comma separated)
    -sv, --skip_vcf                 Skip the vcf file during the download
    -sg, --skip_gtf                 Skip the gtf file during the download
    -sp, --skip_protein             Skip the protein fasta file during download
    -sc, --skip_cds                 Skip the CDS file download
    -sn, --skip_ncrna               Skip the ncRNA file download
    -sdn, --skip_cdna               Skip the cDNA file download
    -sd, --skip_dna                 Skip the DNA file download
    -h, --help                      Show this message and exit.
```

#### Examples

- List all species without downloading any data:

    ```bash
    pgatk ensembl-downloader -l -sv -sg -sp -sc -sd -sn
    ```

- Download all files except cDNA for Turkey (species id=9103):

    ```bash
    pgatk ensembl-downloader -t 9103 -sd -o ensembl_files
    ```

!!! note
    By default the command `ensembl-downloader` downloads all datasets for all species from the latest ENSEMBL release. To limit the download to a particular species specify the species identifier using the `-t` option. To list all available species run the command with `-l` (`--list_taxonomies`) option.

!!! note
    Any of the file types can be skipped using the corresponding option. For example, to avoid downloading the protein sequence fasta file, use the argument `--skip_protein`. Also, note that not all file types exist for all species so the downloaded files depend on availability of the dataset in ENSEMBL.

!!! tip "Hint"
    A VCF file per chromosome is downloaded for homo sapiens due to the large file size they have been distributed this way by ENSEMBL. For other species, a single VCF including all chromosomes is downloaded.

### Downloading COSMIC Data

Downloading mutation data from [COSMIC](https://cancer.sanger.ac.uk/cosmic) is performed using the command `cosmic-downloader`. The current command allows users to download the following files:

- Cosmic mutation file (CosmicMutantExport)
- Cosmic all genes (All_COSMIC_Genes)

#### Command Options

```bash
$ pgatk cosmic-downloader -h
Usage: pgatk cosmic-downloader [OPTIONS]

  Required parameters:
    -u, --username TEXT          Username for cosmic database
    -p, --password TEXT          Password for cosmic database

  Optional parameters:
    -c, --config_file TEXT       Configuration file for the ensembl data downloader pipeline
    -o, --output_directory TEXT  Output directory for the peptide databases
    -h, --help                   Show this message and exit.
```

!!! note
    In order to be able to download COSMIC data, the user should provide a user and password. Please first register in the [COSMIC database](https://cancer.sanger.ac.uk/cosmic/register).

#### Examples

- Download `CosmicMutantExport.tsv.gz` and `All_COSMIC_Genes.fasta.gz`:

    ```bash
    pgatk cosmic-downloader -u userName -p passWord -c config/cosmic_config.yaml -o cosmic_files
    ```

### Downloading cBioPortal Data

Downloading mutation data from [cBioPortal](https://www.cbioportal.org/) is performed using the command `cbioportal-downloader`. cBioPortal stores mutation data from multiple studies (https://www.cbioportal.org/datasets). Each dataset in cBioPortal has an associated study_id.

#### Command Options

```bash
$ pgatk cbioportal-downloader -h
Usage: pgatk cbioportal-downloader [OPTIONS]

  Parameters:
    -c, --config_file TEXT          Configuration file for the ensembl data downloader pipeline
    -o, --output_directory TEXT     Output directory for the peptide databases
    -l, --list_studies              Print the list of all the studies in cBioPortal
    -d, --download_study TEXT       Download a specific Study from cBioPortal (use "all" to download all)
    -h, --help                      Show this message and exit.
```

!!! note
    The argument `-l` (`--list_studies`) allows the user to list all the studies stored in cBioPortal. The `-d` (`--download_study`) argument can be used to obtain mutation data from a particular study.

#### Examples

- Download data for study ID [blca_mskcc_solit_2014](https://www.cbioportal.org/study/summary?id=blca_mskcc_solit_2014):

    ```bash
    pgatk cbioportal-downloader -d blca_mskcc_solit_2014 -o cbioportal_files
    ```

- Download data for all studies in cBioPortal:

    ```bash
    pgatk cbioportal-downloader -d all -o cbioportal_files
    ```

If you face issues downloading all studies from cBioPortal using the `cbioportal-downloader`, please download the studies from the [data hub](https://github.com/cBioPortal/datahub/tree/master/public) through `git-lfs` which is used to download large files from GitHub repositories, see [installation instructions](https://github.com/git-lfs/git-lfs/wiki/Installation).

Following [instructions given on the datahub repository](https://github.com/cBioPortal/datahub), download the entire list of datasets using:

```bash
git clone https://github.com/cBioPortal/datahub.git
cd datahub
git lfs install --local --skip-smudge
git lfs pull -I public --include "data_clinical_sample.txt"
git lfs pull -I public --include "data_mutations_mskcc.txt"
```

### Downloading NCBI / ClinVar Data

Downloading NCBI RefSeq annotations and ClinVar variants for human (GRCh38) is performed using the command `ncbi-downloader`. The tool downloads four files:

- RefSeq gene annotations (GTF)
- RefSeq transcript nucleotide sequences (FASTA)
- Assembly report (chromosome name mapping)
- ClinVar variant calls (VCF)

#### Command Options

```bash
$ pgatk ncbi-downloader -h
Usage: pgatk ncbi-downloader [OPTIONS]

  Required parameters:
    -o, --output-dir TEXT     Output directory for downloaded files

  Optional parameters:
    -c, --config_file TEXT    Configuration YAML file
    --force                   Re-download files even if they exist
    -h, --help                Show this message and exit.
```

#### Examples

- Download all NCBI/ClinVar reference files:

    ```bash
    pgatk ncbi-downloader -o ncbi_data
    ```

- Force re-download of all files:

    ```bash
    pgatk ncbi-downloader -o ncbi_data --force
    ```

## Generate Protein Databases

The **pgatk** framework provides a set of tools to generate protein databases in `FASTA` format from DNA sequences, variants, and mutations. Multiple commands are available depending on the data type provided by the user and the public data providers (cBioPortal, COSMIC and ENSEMBL).

### Cosmic Mutations to Protein Sequences

[COSMIC](https://cancer.sanger.ac.uk/cosmic/) the Catalogue of **Human** Somatic Mutations in Cancer is the world's largest source of expert manually curated somatic mutation information relating to human cancers. The command `cosmic-to-proteindb` converts the cosmic somatic mutations file into a protein sequence database file.

#### Command Options

```bash
$ pgatk cosmic-to-proteindb -h
Usage: pgatk cosmic-to-proteindb [OPTIONS]

  Required parameters:
    -in, --input_mutation TEXT   Cosmic Mutation data file
    -fa, --input_genes TEXT      All Cosmic genes
    -out, --output_db TEXT       Protein database including all the mutations

  Optional parameters:
    -c, --config_file TEXT       Configuration file for the cosmic data pipelines
    -f, --filter_column          Column name to use for filtering or splitting mutations (default: "Primary site")
    -a, --accepted_values        Only consider mutations from these groups (default: all)
    -s, --split_by_filter_column Generate a proteinDB output file for each group (default: False)
    -h, --help                   Show this message and exit.
```

The file input of the tool `-in` (`--input_mutation`) is the cosmic mutation data file. The genes file `-fa` (`--input_genes`) contains the original CDS sequence for all genes used by the COSMIC team to annotate the mutations. [Use cosmic-downloader](#downloading-cosmic-data) to obtain the input files from COSMIC.

The output of the tool is a protein fasta file and is written in the path specified by `-out` (`--output_db`).

#### Examples

- Generate cancer-type specific protein databases. For each cancer type in COSMIC generate a protein database based on the Primary site:

    ```bash
    pgatk cosmic-to-proteindb -in CosmicMutantExport.tsv -fa All_COSMIC_Genes.fasta -out cosmic_proteinDB.fa --split_by_filter_column
    ```

- Generate cell-line specific protein databases:

    ```bash
    pgatk cosmic-to-proteindb -in CosmicCLP_MutantExport.tsv -fa All_CellLines_Genes.fasta -out cosmicCLP_proteinDB.fa --split_by_filter_column --filter_column 'Sample name'
    ```

### cBioPortal Mutations to Protein Sequences

The cBioPortal for Cancer Genomics provides visualization, analysis and download of large-scale cancer genomics data sets. The available datasets can be viewed at [https://www.cbioportal.org/datasets](https://www.cbioportal.org/datasets). The command `cbioportal-to-proteindb` converts the cBioPortal mutations file into a protein sequence database file.

#### Command Options

```bash
$ pgatk cbioportal-to-proteindb -h
Usage: pgatk cbioportal-to-proteindb [OPTIONS]

  Required parameters:
    -c, --config_file TEXT           Configuration for cBioportal
    -in, --input_mutation TEXT       Cbioportal mutation file
    -fa, --input_cds TEXT            CDS genes from ENSEMBL database
    -out, --output_db TEXT           Protein database including the mutations

  Optional parameters:
    -f, --filter_column TEXT         Column in the VCF file to be used for filtering or splitting mutations
    -a, --accepted_values TEXT       Limit mutations to specific groups (tissue type, sample name, etc)
    -s, --split_by_filter_column     Generate a proteinDB per group as specified in the filter_column
    -cl, --clinical_sample_file TEXT Clinical sample file with cancer type per sample identifier
    -h, --help                       Show this message and exit.
```

!!! note
    The clinical sample file for each mutation file can be found under the same directory as the mutation file downloaded from cBioportal (It should have at least two columns named: Cancer Type and Sample Identifier). The file is only needed when generating tissue type databases (when `-s` or `-a` is given).

The file input of the tool `-in` (`--input_mutation`) is the cBioPortal mutation data file. An example is given in [cBioPortal downloader](#downloading-cbioportal-data) showing how to obtain the mutations file for a particular study. The CDS sequence for all genes input file `-fa` (`--input_genes`) can be obtained using the ENSEMBL CDS files, see [ENSEMBL downloader](#downloading-ensembl-data).

!!! note
    The cBioPortal mutations are aligned to the hg19 assembly, make sure that the correct genome assembly is selected for the download.

#### Examples

- Translate mutations from `Bladder` samples in studyID `blca_mskcc_solit_2014`:

    ```bash
    pgatk cbioportal-to-proteindb --config_file config/cbioportal_config.yaml --input_cds human_hg19_cds.fa --input_mutation data_mutations_mskcc.txt --clinical_sample_file data_clinical_sample.txt --output_db bladder_proteindb.fa
    ```

### Variants (VCF) to Protein Sequences

Variant Calling Format (VCFv4.1) is a text file representing genomic variants. The `vcf-to-proteindb` command takes a VCF file and a GTF (Gene annotations) file to translate the genomic variants in the VCF that affect protein-coding transcripts.

#### Command Options

```bash
$ pgatk vcf-to-proteindb -h
Usage: pgatk vcf-to-proteindb [OPTIONS]

  Required parameters:
    -c, --config_file TEXT             Configuration for VCF conversion parameters
    -v, --vcf                          VCF file containing the genomic variants
    -g, --gene_annotations_gtf         Gene models in the GTF format
    -f, --input_fasta                  Fasta sequences for the transcripts in the GTF file
    -o, --output_proteindb             Output file to write the resulting variant protein sequences

  Options:
    --translation_table INTEGER        Translation table (Default 1)
    --mito_translation_table INTEGER   Mito_trans_table (default 2)
    --var_prefix TEXT                   String to add as prefix for the variant peptides
    --report_ref_seq                   Also report the reference peptide from overlapping transcripts
    --annotation_field_name TEXT       Annotation field name in INFO column (default: CSQ)
    --af_field TEXT                    Field name for variant allele frequency (default: none)
    --af_threshold FLOAT               Minimum allele frequency threshold
    --transcript_index INTEGER         Index of transcript ID in annotation columns (default: 3)
    --consequence_index INTEGER        Index of consequence in annotation columns (default: 1)
    --include_consequences TEXT         Consider variants with these consequences (default: all)
    --exclude_consequences TEXT         Exclude these consequences (default: downstream_gene_variant,
                                        upstream_gene_variant, intergenic_variant, intron_variant,
                                        synonymous_variant)
    --skip_including_all_cds           Disable automatic translation of transcripts with defined CDS
    --ignore_filters                   Parse all variants regardless of FILTER field
    --accepted_filters TEXT            Accepted filters for variant parsing
    -h, --help                         Show this message and exit.
```

The file input `--vcf` is a VCF file that can be provided by the user or obtained from ENSEMBL using the [ensembl-downloader](#downloading-ensembl-data). The `--gene_annotations_gtf` file can also be obtained with the ensembl-downloader.

The `--input_fasta` file contains the `CDS` and DNA sequences for all genes present in the GTF file. This file can be generated from the GTF file using the [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) tool as follows:

```bash
gffread -F -w input_fasta.fa -g genome.fa gene_annotations_gtf
```

The output of the tool is a protein fasta file written to the path specified by `--output_proteindb`.

#### Examples

- Translate human *missense* variants from ENSEMBL VCFs that have a minimum *AF 5%*:

    ```bash
    pgatk vcf-to-proteindb \
        --vcf homo_sapiens_incl_consequences.vcf \
        --input_fasta transcripts.fa \
        --gene_annotations_gtf genes.gtf \
        --include_consequences missense_variant \
        --af_field MAF \
        --af_threshold 0.05 \
        --output_proteindb var_peptides.fa
    ```

!!! note
    - By default `vcf-to-proteindb` considers transcripts that have a coding sequence that includes all protein_coding genes.
    - By default all consequences are accepted except those given with `--exclude_consequences`. See the list of consequences generated by VEP: [https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)

- Translate human *missense* variants or *inframe_insertion* from gnomAD VCFs that have a minimum 1% allele frequency in control samples:

    ```bash
    pgatk vcf-to-proteindb \
        --vcf gnomad_genome.vcf \
        --input_fasta gencode.fa \
        --gene_annotations_gtf gencode.gtf \
        --include_consequences missense_variant,frameshift_insert \
        --annotation_field_name vep \
        --af_threshold 0.01 \
        --af_field control_af \
        --transcript_index 6
    ```

!!! tip "Hint"
    - `vcf-to-proteindb` considers transcripts that have a coding sequence which includes all *protein_coding* transcripts.
    - The provided VCF file has some specific properties: the annotation field is specified with the string *vep* hence the `--annotation_field_name` parameter, the transcript is at the sixth position in the annotation field, and since gnomAD collects variants from many sources it provides allele frequencies across many sub-populations. In this case we use only variants common within control samples therefore the `--af_field` is set to `control_af`.
    - Since gnomAD uses GENCODE gene annotations, you need to change the default `biotype_str` from *transcript_biotype* to *transcript_type* (as written in the GTF file).

!!! note
    When ENSEMBL data is used, the default options should work. However, for other data sources such as variants from gnomAD, GTF from GENCODE and others one or more of the following parameters need to be changed: `--af_field`, `--annotation_field_name`, `--transcript_index`, `--consequence_index`.

- Translate human variants from a custom VCF obtained from sequencing of a sample:

    ```bash
    pgatk vcf-to-proteindb \
        --vcf sample.vcf \
        --input_fasta transcripts.fa \
        --gene_annotations_gtf genes.gtf \
        --annotation_field_name '' \
        --output_proteindb var_peptides.fa
    ```

### ClinVar Variants to Protein Sequences

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is NCBI's public archive of reports on the relationships among human variations and phenotypes. The command `clinvar-to-proteindb` generates a variant protein database from ClinVar VCF files using NCBI RefSeq gene models. Unlike `vcf-to-proteindb`, this command does **not** require VEP annotations.

#### Command Options

```bash
$ pgatk clinvar-to-proteindb -h
Usage: pgatk clinvar-to-proteindb [OPTIONS]

  Required parameters:
    -v, --vcf TEXT              ClinVar VCF file path
    -g, --gtf TEXT              NCBI RefSeq GTF file path
    -f, --fasta TEXT            RefSeq transcript nucleotide FASTA file path
    -a, --assembly-report TEXT  NCBI assembly report file path
    -o, --output TEXT           Output protein FASTA file path

  Optional parameters:
    -c, --config_file TEXT      Configuration YAML file
    -h, --help                  Show this message and exit.
```

The input files can be downloaded using the [ncbi-downloader](#downloading-ncbi--clinvar-data) command.

#### Examples

- Generate a protein database from ClinVar variants:

    ```bash
    pgatk clinvar-to-proteindb \
        --vcf ncbi_data/clinvar.vcf.gz \
        --gtf ncbi_data/GRCh38_latest_genomic.gtf.gz \
        --fasta ncbi_data/GRCh38_latest_rna.fna.gz \
        --assembly-report ncbi_data/GRCh38_latest_assembly_report.txt \
        --output clinvar_proteins.fa
    ```

### Transcripts (DNA) to Protein Sequences

DNA sequences given in a FASTA format can be translated using the `dnaseq-to-proteindb` tool. This tool allows for translation of all kinds of transcripts (coding and noncoding) by specifying the desired biotypes.

The most suited `--input_fasta` file can be generated from a given GTF file using the `gffread` command as follows:

```bash
gffread -F -w transcript_sequences.fa -g genome.fa gene_annotations_gtf
```

The FASTA file generated from the GTF file would contain DNA sequences for all transcripts regardless of their biotypes. It also specifies the CDS positions for the protein coding transcripts. The `dnaseq-to-proteindb` command recognizes features such as biotype and expression values in the FASTA header that are taken from the GTF INFO field (if available). However, it is not required to have those in the FASTA header but their presence enables the user to filter by biotype and expression values during the translation step.

#### Command Options

```bash
$ pgatk dnaseq-to-proteindb -h
Usage: pgatk dnaseq-to-proteindb [OPTIONS]

  Required parameters:
    -c, --config_file TEXT             Configuration for VCF conversion parameters
    --input_fasta                      Fasta sequences for the transcripts
    --output_proteindb                 Output file to write the resulting protein sequences

  Optional parameters:
    --translation_table INTEGER        Translation Table (default 1)
    --num_orfs INTEGER                 Number of ORFs (default 0)
    --num_orfs_complement INTEGER      Number of ORFs from the reverse side (default 0)
    --skip_including_all_cds           Disable automatic translation of transcripts with defined CDS
    --include_biotypes TEXT            Translate sequences with specified biotypes (default: protein coding)
    --exclude_biotypes TEXT            Skip sequences with unwanted biotypes (default: None)
    --biotype_str TEXT                 String used to identify gene/transcript biotype (default: transcript_biotype)
    --expression_str TEXT              String for extracting expression value (default: None)
    --expression_thresh FLOAT          Threshold for expression value filtering (default: 5)
    --var_prefix TEXT                  Prefix to be added to fasta headers (default: none)
    -h, --help                         Show this message and exit.
```

#### Examples

- Generate the canonical protein database (translate all *protein_coding* transcripts):

    ```bash
    pgatk dnaseq-to-proteindb \
        --config_file config/ensembl_config.yaml \
        --input_fasta transcript_sequences.fa \
        --output_proteindb proteindb_from_CDSs_DNAseq.fa
    ```

- Generate a protein database from lincRNA and canonical proteins:

    ```bash
    pgatk dnaseq-to-proteindb \
        --config_file config/ensembl_config.yaml \
        --input_fasta transcript_sequences.fa \
        --output_proteindb proteindb_from_lincRNA_canonical_sequences.fa \
        --var_prefix lincRNA_ \
        --include_biotypes lincRNA
    ```

- Generate a protein database from processed pseudogene:

    ```bash
    pgatk dnaseq-to-proteindb \
        --config_file config/ensembl_config.yaml \
        --input_fasta transcript_sequences.fa \
        --output_proteindb proteindb_from_processed_pseudogene.fa \
        --var_prefix pseudogene_ \
        --include_biotypes processed_pseudogene,transcribed_processed_pseudogene,translated_processed_pseudogene \
        --skip_including_all_cds
    ```

- Generate alternative ORFs from canonical sequences:

    ```bash
    pgatk dnaseq-to-proteindb \
        --config_file config/ensembl_config.yaml \
        --input_fasta transcript_sequences.fa \
        --output_proteindb proteindb_from_altORFs.fa \
        --var_prefix altorf_ \
        --include_biotypes altORFs \
        --skip_including_all_cds
    ```

- Generate protein sequences (six-frame translation) from a Genome assembly:

    ```bash
    pgatk dnaseq-to-proteindb \
        --config_file config/ensembl_config.yaml \
        --input_fasta genome.fa \
        --output_proteindb proteindb_genome.fa \
        --biotype_str '' \
        --num_orfs 3 \
        --num_orfs_complement 3
    ```

### Generate Decoy Database

`generate-decoy` command enables generation of decoy databases for any given protein sequence database. Decoy databases are needed to evaluate significance of spectra-sequence matching scores in proteomics mass spectrometry experiments.

*DecoyPYrat* is integrated into `pgatk` as the standard method for generating decoy sequences. In addition to reversing the target sequences, the tool replaces the cleavage with preceding amino acids. Also, it checks for the presence of the reversed sequence in the target sequences and if found, *DecoyPYrat* shuffles the sequences to avoid target-decoy sequence matches. For more information please read the *DecoyPYrat* manual available at: [https://www.sanger.ac.uk/science/tools/decoypyrat](https://www.sanger.ac.uk/science/tools/decoypyrat).

#### Command Options

```bash
$ pgatk generate-decoy -h
Usage: pgatk generate-decoy [OPTIONS]

  Required parameters:
    -c, --config_file TEXT          Configuration file for decoy generation
    -o, --output TEXT               Output file for decoy database
    -i, --input TEXT                FASTA file of target protein sequences (*.fasta|*.fa)

  Optional parameters:
    -s, --cleavage_sites TEXT       Amino acids at which to cleave (Default: KR)
    -a, --anti_cleavage_sites TEXT  Amino acids at which not to cleave if following cleavage site
    -p, --cleavage_position TEXT    Cleavage position [c, n] (Default: c)
    -l, --min_peptide_length INTEGER  Minimum peptide length to compare (Default: 5)
    -n, --max_iterations INTEGER    Max shuffle iterations (Default: 100)
    -x, --do_not_shuffle TEXT       Turn OFF shuffling of decoy peptides (Default: false)
    -w, --do_not_switch TEXT        Turn OFF switching of cleavage site (Default: false)
    -d, --decoy_prefix TEXT         Accession prefix for decoy proteins (Default: DECOY_)
    -t, --temp_file TEXT            Temporary file for decoys prior to shuffling
    -b, --no_isobaric TEXT          Do not make decoy peptides isobaric (Default: false)
    -m, --memory_save TEXT          Slower but uses less memory (Default: false)
    -h, --help                      Show this message and exit.
```

#### Examples

- Generate decoy sequences for a protein database:

    ```bash
    pgatk generate-decoy -c config/protein_decoy.yaml --input proteindb.fa --output decoy_proteindb.fa
    ```

## Post-Processing Utilities

### Digest and Filter Variant Peptides

The `digest-mutant-protein` command performs *in silico* tryptic digestion of mutant proteins and filters out any peptides that also appear in a canonical reference proteome. The result is a compact FASTA of variant-specific peptides only.

#### Command Options

```bash
$ pgatk digest-mutant-protein -h
Usage: pgatk digest-mutant-protein [OPTIONS]

  Required parameters:
    -i, --input TEXT        Input mutant protein FASTA file(s), comma-separated
    -f, --fasta TEXT        Reference canonical protein FASTA
    -o, --output TEXT       Output file for unique variant peptides

  Optional parameters:
    --prefix TEXT           Header prefix for output entries (default: Mutation)
    --min-len INTEGER       Minimum peptide length (default: 7)
    --max-len INTEGER       Maximum peptide length (default: 40)
    --missed-cleavages INTEGER  Number of missed cleavages (default: 0)
    -h, --help              Show this message and exit.
```

#### Examples

- Digest variant proteins and keep only variant-specific peptides:

    ```bash
    pgatk digest-mutant-protein \
        --input variant_proteins.fa \
        --fasta canonical_proteins.fa \
        --output unique_variant_peptides.fa \
        --min-len 7 \
        --max-len 40 \
        --missed-cleavages 2
    ```

- Combine and filter multiple variant sources:

    ```bash
    pgatk digest-mutant-protein \
        --input cosmic_proteins.fa,clinvar_proteins.fa,ensembl_variants.fa \
        --fasta canonical_proteins.fa \
        --output combined_variant_peptides.fa
    ```

### Map Peptides to Genomic Coordinates

The `map-peptide2genome` command maps proteomics-identified peptides to genomic coordinates and produces a GFF3 file suitable for visualization in genome browsers (UCSC, IGV, Ensembl).

#### Command Options

```bash
$ pgatk map-peptide2genome -h
Usage: pgatk map-peptide2genome [OPTIONS]

  Required parameters:
    -i, --input TEXT        Input peptide identification TSV
    -g, --gtf TEXT          GTF gene annotation file
    -f, --fasta TEXT        Protein FASTA file
    -m, --idmap TEXT        Protein-to-transcript ID mapping file
    -o, --output TEXT       Output GFF3 file

  Optional parameters:
    --pep-col INTEGER       Peptide column index, 0-based (default: 0)
    --prot-col INTEGER      Protein column index, 0-based (default: 1)
    -h, --help              Show this message and exit.
```

#### Examples

- Map identified peptides to genomic coordinates:

    ```bash
    pgatk map-peptide2genome \
        --input peptide_identifications.tsv \
        --gtf genes.gtf \
        --fasta proteins.fa \
        --idmap protein_to_transcript.tsv \
        --output peptides.gff3
    ```
