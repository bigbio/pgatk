from click.testing import CliRunner

from pypgatk.pypgatk_cli import cli


def vcf_to_proteindb():
    """
    Test the default behaviour of the vcf-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['vcf-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--vep_annotated_vcf', 'testdata/test.vcf',
                            '--input_fasta', 'testdata/test.fa',
                            '--gene_annotations_gtf', 'testdata/test.gtf',
                            '--var_prefix', 'ensvar',
                            '--af_field', 'MAF',
                            '--output_proteindb', 'testdata/proteindb_from_ENSEMBL_VCF.fa'])
    assert result.exit_code == 0


def vcf_gnomad_to_proteindb():
    """
    Test the default behaviour of the vcf-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['vcf-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--vep_annotated_vcf', 'testdata/test_gnomad.vcf',
                            '--input_fasta', 'testdata/test_gencode.fa',
                            '--gene_annotations_gtf', 'testdata/test_gencode.gtf',
                            '--output_proteindb', 'testdata/proteindb_from_gnomad_VCF.fa',
                            '--af_field', 'controls_AF',
                            '--var_prefix', 'gnomvar',
                            '--transcript_index', 6,
                            '--biotype_str', 'transcript_type',
                            '--annotation_field_name', 'vep'])
    assert result.exit_code == 0


def dnaseq_to_proteindb():
    """
    Test the default behaviour (translate CDSs) of the dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_CDSs_DNAseq.fa'])
    assert result.exit_code == 0

def dnaseq_ncrnas_to_proteindb():
    """
    Test generation of proteinDB from short noncoding RNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_ncRNAs_DNAseq.fa',
                            '--var_prefix', 'ncRNA_',
                            '--include_biotypes',
                            'lncRNA,retained_intron,Mt_rRNA,Mt_tRNA,miRNA,misc_RNA,rRNA,ribozyme,sRNA,scRNA,scaRNA,snRNA,snoRNA,vaultRNA',
                            '--skip_including_all_cds'])
    assert result.exit_code == 0

def dnaseq_lncrnas_to_proteindb():
    """
    Test generation of proteinDB from long noncoding RNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_lncRNAs_DNAseq.fa',
                            '--var_prefix', 'lncRNA_',
                            '--include_biotypes',
                            'lncRNA,retained_intron',
                            '--skip_including_all_cds'])
    assert result.exit_code == 0


def dnaseq_sncrnas_to_proteindb():
    """
    Test generation of proteinDB from short noncoding RNAs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_sncRNAs_DNAseq.fa',
                            '--var_prefix', 'sncRNA_',
                            '--include_biotypes',
                            'Mt_rRNA,Mt_tRNA,miRNA,misc_RNA,rRNA,ribozyme,sRNA,scRNA,scaRNA,snRNA,snoRNA,vaultRNA',
                            '--skip_including_all_cds'])
    assert result.exit_code == 0


def dnaseq_pseudogenes_to_proteindb():
    """
    Test generation of proteinDB from pseudogenes using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_pseudogenes_DNAseq.fa',
                            '--var_prefix', 'pseudo_',
                            '--include_biotypes',
                            'disrupted_domain, IGC_pseudogene, IGJ_pseudogene, IG_pseudogene, IGV_pseudogene, processed_pseudogene, transcribed_processed_pseudogene, transcribed_unitary_pseudogene, transcribed_unprocessed_pseudogene, translated_processed_pseudogene, TRJ_pseudogene, unprocessed_pseudogene',
                            '--skip_including_all_cds'])
    assert result.exit_code == 0


def dnaseq_altorfs_to_proteindb():
    """
    Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                            '--input_fasta', 'testdata/test.fa',
                            '--output_proteindb', 'testdata/proteindb_from_altORFs_DNAseq.fa',
                            '--var_prefix', 'altorf_',
                            '--include_biotypes', 'altORFs', '--skip_including_all_cds'])
    assert result.exit_code == 0


def cbioportal_to_proteindb():
    """
    Test generation proteinDB from cBioportal mutations using cbioportal-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['cbioportal-to-proteindb', '--config_file', 'config/cbioportal_config.yaml',
                            '--input_mutation', 'testdata/test_cbioportal_data_mutations_mskcc.txt',
                            '--input_cds', 'testdata/test_cbioportal_genes.fa',
                            '--output_db', 'testdata/test_cbioportal_data_mutations_mskcc_proteindb.fa',
                            '--clinical_sample_file', 'testdata/test_cbioportal_data_clinical_sample.txt',
                            '--filter_column', 'Tumor_Sample_Barcode'
                            '--split_by_filter_column', '--accepted_values', 'all'])
    assert result.exit_code == 0


def cosmic_to_proteindb():
    """
    Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
    :return:
    """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['cosmic-to-proteindb', '--config_file', 'config/cosmic_config.yaml',
                            '--input_mutation', 'testdata/test_cosmic_mutations.tsv',
                            '--input_genes', 'testdata/test_cosmic_genes.fa',
                            '--output_db', 'testdata/test_cosmic_mutations_proteindb.fa',
                            '--filter_column', 'Primary site',
                            '--split_by_filter_column', '--accepted_values', 'all'])
    assert result.exit_code == 0

def generate_decoy_database():
    """
        Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
        :return:
        """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['generate-decoy', '--config_file', 'config/protein_decoy.yaml',
                            '--input', 'testdata/test_db.fa', '--output', 'testdata/output_decoy.fa'])
    assert result.exit_code == 0

def download_ensembl_data():
    """
        Test downloading ensembl data for species with taxonomy identifer = 9103
        :return:
        """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['ensembl-downloader', '--config_file', 'config/ensembl_downloader_config.yaml',
                            '--taxonomy', '9103', '--output_directory', 'testdata'])
    assert result.exit_code == 0


def download_ensembl_data_37():
  """
      Test downloading ensembl data for species with taxonomy identifier = 9606
      :return:
      """
  runner = CliRunner()
  result = runner.invoke(cli,
                         ['ensembl-downloader', '--config_file', 'config/ensembl_downloader_config.yaml',
                          '--taxonomy', '9606', '--output_directory', 'testdata', '--grch37'])
  assert result.exit_code == 0

def download_cbioportal_data():
    """
        Test downloading cbioportal data for study id: paac_jhu_2014
        :return:
        """
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['cbioportal-downloader', '--config_file', 'config/cbioportal_config.yaml',
                            '--download_study', 'paac_jhu_2014', '--output_directory', 'testdata'])
    assert result.exit_code == 0


def check_ensembl_database():
  runner = CliRunner()
  result = runner.invoke(cli,
                         ['ensembl-check', '--config_file', 'config/ensembl_config.yaml',
                          '--input_fasta', 'testdata/proteindb_from_ENSEMBL_VCF.fa', '--output', 'testdata/proteindb_from_ENSEMBL_VCF-clean.fa', '--add_stop_codons', '--num_aa', '6'])
  assert result.exit_code == 0



if __name__ == '__main__':
    vcf_to_proteindb()
    vcf_gnomad_to_proteindb()
    dnaseq_to_proteindb()
    dnaseq_ncrnas_to_proteindb()
    dnaseq_lncrnas_to_proteindb()
    dnaseq_sncrnas_to_proteindb()
    dnaseq_pseudogenes_to_proteindb()
    dnaseq_altorfs_to_proteindb()
    cbioportal_to_proteindb()
    generate_decoy_database()
    cosmic_to_proteindb()
    download_ensembl_data()
    download_ensembl_data_37()
    download_cbioportal_data()
    check_ensembl_database()
