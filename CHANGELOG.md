Changelog
=========


(unreleased)
------------
- Change command config path. [ypriverol]
- Ensembl downloader updated. [ypriverol]
- Increase package version. [ypriverol]
- Code clean. [ypriverol]
- New downloader for previous versions of ensembl. [ypriverol]
- Change command config path. [ypriverol]
- Better error handling for cosmic downloads. [ypriverol]
- Commandline ordered. [ypriverol]
- Handle the stop codons into split proteins instead of removing the
  stop codons character * [ypriverol]
- Change version. [Yasset Perez-Riverol]

  change version
- Update setup.py. [Yasset Perez-Riverol]
- Fixed typo. [husensofteng]
- Updated cosmic and cbioportal to proteindb tests using filter_column
  instread of -t. [husensofteng]
- Made changes to mutations to protein db allowing filtering mutations
  file by any column instead of just tissue type, this also enables
  translating cosmic cell lines data. [husensofteng]
- Merge remote-tracking branch 'origin/master' [yperez]
- Merge pull request #27 from bigbio/dependabot/pip/cryptography-3.2.
  [Yasset Perez-Riverol]

  Bump cryptography from 2.6.1 to 3.2
- Bump cryptography from 2.6.1 to 3.2. [dependabot[bot]]

  Bumps [cryptography](https://github.com/pyca/cryptography) from 2.6.1 to 3.2.
  - [Release notes](https://github.com/pyca/cryptography/releases)
  - [Changelog](https://github.com/pyca/cryptography/blob/master/CHANGELOG.rst)
  - [Commits](https://github.com/pyca/cryptography/compare/2.6.1...3.2)
- Cosmic cell-lines downloaded implemented. [yperez]
- Fixed error in parameter reading for biotype, num_orf in ensembl.py.
  [husensofteng]
- Merge remote-tracking branch 'origin/master' [yperez]
- Create pythonpublish.yml. [Yasset Perez-Riverol]
- Travis removed. [yperez]
- Merge remote-tracking branch 'origin/master' [yperez]
- Create pythonpackage.yml. [Yasset Perez-Riverol]
- Major changes in package structure. [yperez]
- Major changes in package structure. [yperez]
- Major changes in package structure. [yperez]
- Major changes in package structure. [yperez]
- Create pythonapp.yml. [Yasset Perez-Riverol]
- Major changes in package structure. [yperez]
- Added the ESEMBL name as parameter. [Yasset Perez-Riverol]
- Minor changes. [Yasset Perez-Riverol]
- Added package to dependencies. [Yasset Perez-Riverol]
- Skip records that have either alt or ref equal to None. [husensofteng]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk.
  [husensofteng]
- Added package to dependencies. [Yasset Perez-Riverol]
- Added package to dependencies. [Yasset Perez-Riverol]
- Added package to dependencies. [Yasset Perez-Riverol]
- Implement better call retrieve system for ensembl. [Yasset Perez-
  Riverol]
- Changes in the setup file to be installable. [Yasset Perez-Riverol]
- Changes in the setup file to be installable. [Yasset Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Print when error in alt-debugging. [husensofteng]
- Updated readme. [husensofteng]
- Updated readme. [husensofteng]
- Updated readme. [husensofteng]
- Updated readme. [husensofteng]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk.
  [husensofteng]
- Rename README.rst to README.md. [Husen M. Umer]
- Rename README.md to README.rst. [Husen M. Umer]
- Updated readme. [husensofteng]
- Updated readme. [husensofteng]
- Updated readme refer to readthedocs. [husensofteng]
- Updated readme to match the docs. [husensofteng]
- Skip vcf recrods that have transcripts not found in gtf file.
  [husensofteng]
- Test result files. [husensofteng]
- Minor. [husensofteng]
- Re-ran test cases. [husensofteng]
- Re-ran test cases. [husensofteng]
- Added var_prefix param to the test cases. [husensofteng]
- Added var_prefix param to dnaseq-to-proteindb command to allow for
  prefix in fasta heades. [husensofteng]
- Keep only alpha chars in tissue types for file names in cosmic
  proteindb. [husensofteng]
- Fixed bug in cbiportal tissue specific databases. [husensofteng]
- Fixed bug. [husensofteng]
- Minor fix. [husensofteng]
- Check for header continously in cbioportal files. [husensofteng]
- Minor fix. [husensofteng]
- Minor fix in cbiportal-to-proteindb. [husensofteng]
- Removed print message. [husensofteng]
- Updated tests accordingly in pypgatk_tests.py. [husensofteng]
- Allowed handling gencode fasta headers in vcf_to_proteindb by applying
  key_function to fasta input as well as in dnaseq_to_proteindb by
  adding a '|' separator when space is not found in the fasta headers.
  [husensofteng]
- Changed default include_biotypes to all protein coding genes, changed
  default af_field to None. [husensofteng]
- Added alias for parameters in vcf_to_proteindb, changed default
  include_biotypes to all protein coding genes. [husensofteng]
- Error printing update in cgenomes_proteindb.py. [husensofteng]
- More error reporting in cgenomes_proteindb.py and skip comment lines.
  [husensofteng]
- Minor fix. [husensofteng]
- Ran tests case using Turkey as the species for ensemble data download.
  [husensofteng]
- Implemented get vcf function in ensembl-downloader. [husensofteng]
- Fixed error, refactor and tests cases for data download.
  [husensofteng]
- Minor code refactory. [husensofteng]
- Fixed error in data_downloader. [husensofteng]
- Fixed cbioportal-to-proteindb error and ran the related test case.
  [husensofteng]
- Changed error.code to error_code. [husensofteng]
- Added package name to __init. [husensofteng]
- Changes in cbioportal download. [Yasset Perez-Riverol]
- Changes in cbioportal download. [Yasset Perez-Riverol]
- Minor changes in the pypi release. [Yasset Perez-Riverol]
- Minor changes in the cbioportal downlaoder. [Yasset Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Small change in the code to remove _ from prefix concat. [Yasset
  Perez-Riverol]
- Changes in the code to correct bug. [Yasset Perez-Riverol]
- Keep original protein id in decoy sequences. [Yafeng Zhu]
- Decoy generation working ok - with tests included. [Yasset Perez-
  Riverol]
- Decoy commandline fixed. [Yasset Perez-Riverol]
- Minor changes in the pipeline to add decoy. [Yasset Perez-Riverol]
- First iteration decoy db. [Yasset Perez-Riverol]
- Minor changes in the commandline structure - code clean. [Yasset
  Perez-Riverol]
- Minor changes in the code to download cdna. [Yasset Perez-Riverol]
- Minro changes. [Yasset Perez-Riverol]
- Minor changes. [Yasset Perez-Riverol]
- Minor changes. [Yasset Perez-Riverol]
- Minor changes in the string format. [Yasset Perez-Riverol]
- Minor changes in the string format. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Travis done. [Yasset Perez-Riverol]
- Fixed list_studies issue in cbioportal_downloader. [husensofteng]
- Added supported for tissue type databases. [husensofteng]
- Added support for generating tissue type specific db from COSMIC.
  [husensofteng]
- Minor update. [husensofteng]
- Corrected paths in pypgatk_tests.py. [husensofteng]
- Minor changes in code. [Yasset Perez-Riverol]
- Minor changes in the library, travis added. [Yasset Perez-Riverol]
- Major refactoring of the structure of the project. [Yasset Perez-
  Riverol]
- Major refactoring of the structure of the project. [Yasset Perez-
  Riverol]
- Major refactoring of the structure of the project. [Yasset Perez-
  Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Updated ensembl.py: wrote _i to seq ID number ORFs. [husensofteng]
- Fixed errors and ran test cases. [husensofteng]
- Corrected the merge issues. [husensofteng]
- Added support vcf variant filtering, testdata and testcases.
  [husensofteng]
- Minor code changes. [Yasset Perez-Riverol]
- Merge changes. [Yasset Perez-Riverol]
- Added test cases for generating lncRNAs, pseudogenes, altORfs,
  gnomad_vcf. [husensofteng]
- Added support for generating altORF in dnaseq_to_proteindb.
  [husensofteng]
- Minor changes in the ensembl.py. [Yasset Perez-Riverol]
- Changed gtf variable. [husensofteng]
- Corrected vep to vcf. [husensofteng]
- Upadted vcf_to_proteindb and dnaseq_to_proteindb. [husensofteng]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Deleted old TranscriptToProteinDB.py. [Yasset Perez-Riverol]
- Minor changes, new VCF transcript to proteindb added. [Yasset Perez-
  Riverol]
- Merge changes in the tanslation table. [Yasset Perez-Riverol]
- Added back TRANSLATION_TABLE option to vep_to_proteindb.
  [husensofteng]
- Minor changes. [Yasset Perez-Riverol]
- Minor changes in the tests that enables running default vep to
  proteindb. [Yasset Perez-Riverol]
- The VEP VCF to ProteinDB ready for tests. [Yasset Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Add nextflow script for database search. [Yafeng Zhu]
- Changed test file name. [husensofteng]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Minor change in the library. [Yasset Perez-Riverol]


1.1.0 (2019-05-10)
------------------
- Added tests dir. [husensofteng]
- Updated params. [husensofteng]
- Updated ensembl.py. [husensofteng]
- Updated comment. [husensofteng]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk.
  [husensofteng]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk. [Yasset
  Perez-Riverol]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk :merge
  some small changes. [Yasset Perez-Riverol]
- Minor changes in the main vcf to proteindb tool. [Yasset Perez-
  Riverol]
- Minor changes to move VCF to proteinDB. [Yasset Perez-Riverol]
- Translate transcripts. [husensofteng]
- Removed genome_fasta param. [husensofteng]
- Added param for fasta input file. [husensofteng]
- Return none in case of no overlap. [husensofteng]
- Changed default to consider all CDSs. [husensofteng]
- Parse VCF and generate mutated protein sequences based on annotated
  transcripts. [husensofteng]
- Minor changes in the 3'frame translation script. [Yasset Perez-
  Riverol]
- Script to map peptides back to genome, generating a GFF file output.
  [Yafeng Zhu]
- Script for translate transcriptome sequences in three reading frames.
  [Yafeng Zhu]
- Minor changes in the README. [Yasset Perez-Riverol]
- Minor changes in the README. [Yasset Perez-Riverol]
- Minor changes in the README and the command names. [Yasset Perez-
  Riverol]
- Minor changes in the README. [Yasset Perez-Riverol]
- Minor changes in the README. [Yasset Perez-Riverol]
- Improvements in the README. [Yasset Perez-Riverol]
- Improvements in the README. [Yasset Perez-Riverol]
- Remove an script. [Yasset Perez-Riverol]
- Minor changes in the README. [Yasset Perez-Riverol]
- Script move into temp db folder. [Yasset Perez-Riverol]
- Translate variants from VEP-annotated VCF based on the corresponding
  GTF. [husensofteng]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk. [Yasset
  Perez-Riverol]
- Minor changes reviewed and corrected. [Yafeng Zhu]
- Improve README. [Yasset Perez-Riverol]
- Minor changes in the library, script to translate mutation into
  proteindb. [Yasset Perez-Riverol]
- Minor changes in the code. [Yasset Perez-Riverol]
- Cosmic2proteindb has been implemented as a pipeline. [Yasset Perez-
  Riverol]
- Command cli updated. [Yasset Perez-Riverol]
- Minor changes. [Yasset Perez-Riverol]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk. [Yafeng
  Zhu]
- Minor changes enabling download of COSMIC files. [Yasset Perez-
  Riverol]
- Minor changes enabling download of COSMIC files. [Yasset Perez-
  Riverol]
- Minor changes in wrong commit. [Yasset Perez-Riverol]
- Cds and ncrna included. [Yasset Perez-Riverol]
- Minor changes in the pipelines. [Yasset Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Minor changes in the pipeline. [Yasset Perez-Riverol]
- Config for cbioportal. [Yasset Perez-Riverol]
- Some minor changes. [Yafeng Zhu]
- Rename. [Yafeng Zhu]
- Script for convert COSMIC annotated mutations to mutant protein
  sequences. [Yafeng Zhu]
- Enable command line options. [Yafeng Zhu]
- Minor changes in the code style. [Yasset Perez-Riverol]
- Merge branch 'master' of https://github.com/bigbio/py-pgatk. [Yasset
  Perez-Riverol]
- Updated TRYPSIN func. [Husen M. Umer]

  To allow for any number of missed cleavages for the digestions.
- Minor changes, converted into a single tool app. [Yasset Perez-
  Riverol]
- Minor changes in previous scripts. [Yasset Perez-Riverol]
- Remove intellIJ. [Yasset Perez-Riverol]
- Remove previous libraries. [Yasset Perez-Riverol]
- Minor changes in the library to download all taxonomies. [Yasset
  Perez-Riverol]
- Some refactoring of the code to allow taxonomy id download. [Yasset
  Perez-Riverol]
- Initial commit. [Yafeng Zhu]
- Expect KeyError on ENSEMBL dict. [husensofteng]
- Added func to get gtf files. [husensofteng]
- Ignore eclipse project files. [husensofteng]
- Get canonical peptide sequences for all species. [husensofteng]
- Get canoncial peptide files for all species. [husensofteng]
- Update README.md. [Yasset Perez-Riverol]
- Initial commit. [Yasset Perez-Riverol]


