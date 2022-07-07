import unittest

from click.testing import CliRunner

from pypgatk.pypgatk_cli import cli


class PypgatkRunnerTests(unittest.TestCase):

    def test_peptide_class_group_fdr(self):
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['peptide-class-fdr',
                                '-in', 'testdata/20151020_QE3_UPLC8_DBJ_SA_HCT116_Rep2_46frac_10_consensus.idxml',
                                '-out',
                                'testdata/20151020_QE3_UPLC8_DBJ_SA_HCT116_Rep2_46frac_10_consensus_filter.idxml',
                                '--peptide-groups-prefix',
                                '"{non_canonical:[altorf,pseudo,ncRNA];mutations:[COSMIC,cbiomut];variants:[var_mut,var_rs]}"'])
        self.assertEqual(result.exit_code, 0)

    def test_peptide_classes_fdr(self):
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['peptide-class-fdr',
                                '-in', 'testdata/20151020_QE3_UPLC8_DBJ_SA_HCT116_Rep2_46frac_10_consensus.idxml',
                                '-out',
                                'testdata/20151020_QE3_UPLC8_DBJ_SA_HCT116_Rep2_46frac_10_consensus_filter.idxml',
                                '--peptide-classes-prefix', '"altorf,pseudo,ncRNA,COSMIC,cbiomut,var_mut,var_rs"'])
        self.assertEqual(result.exit_code, 0)

    def test_vcf_to_proteindb(self):
        """
        Test the default behaviour of the vcf-to-proteindb tool
        :return:
        """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['vcf-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                                '--vcf', 'testdata/test.vcf',
                                '--input_fasta', 'testdata/test.fa',
                                '--gene_annotations_gtf', 'testdata/test.gtf',
                                '--var_prefix', 'ensvar',
                                '--af_field', 'MAF',
                                '--output_proteindb', 'testdata/proteindb_from_ENSEMBL_VCF.fa',
                                '--annotation_field_name', 'CSQ',
                                '--biotype_str', 'feature_type',
                                '--include_biotypes', 'mRNA,ncRNA'])
        self.assertEqual(result.exit_code, 0)

    def test_vcf_to_proteindb_notannotated(self):
        """
          Test the default behaviour of the vcf-to-proteindb tool using not-annotated vcf
          :return:
          """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['vcf-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                                '--vcf', 'testdata/test.vcf',
                                '--input_fasta', 'testdata/test.fa',
                                '--gene_annotations_gtf', 'testdata/test.gtf',
                                '--var_prefix', 'varsample',
                                '--output_proteindb', 'testdata/proteindb_from_custom_VCF.fa',
                                '--annotation_field_name', ""])
        if result.exit_code != 0:
            print(result)
        self.assertEqual(result.exit_code, 0)

    def test_vcf_gnomad_to_proteindb(self):
        """
          Test the default behaviour of the vcf-to-proteindb tool
          :return:
          """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['vcf-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                                '--vcf', 'testdata/test_gnomad.vcf',
                                '--input_fasta', 'testdata/test_gencode.fa',
                                '--gene_annotations_gtf', 'testdata/test_gencode.gtf',
                                '--output_proteindb', 'testdata/proteindb_from_gnomad_VCF.fa',
                                '--af_field', 'controls_AF',
                                '--af_threshold', '0.001',
                                '--var_prefix', 'gnomvar',
                                '--annotation_field_name', 'vep', ])
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_to_proteindb(self):
        """
          Test the default behaviour (translate CDSs) of the dnaseq-to-proteindb tool
          :return:
          """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['dnaseq-to-proteindb', '--config_file', 'config/ensembl_config.yaml',
                                '--input_fasta', 'testdata/test.fa',
                                '--output_proteindb', 'testdata/proteindb_from_CDSs_DNAseq.fa'])
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_ncrnas_to_proteindb(self):
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
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_lncrnas_to_proteindb(self):
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
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_sncrnas_to_proteindb(self):
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
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_pseudogenes_to_proteindb(self):
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
        self.assertEqual(result.exit_code, 0)

    def test_dnaseq_altorfs_to_proteindb(self):
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
        self.assertEqual(result.exit_code, 0)

    def test_cbioportal_to_proteindb(self):
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
                                '--filter_column', 'CANCER_TYPE'
                                                   '--split_by_filter_column', '--accepted_values', 'all'])
        self.assertEqual(result.exit_code, 0)

    def test_cosmic_to_proteindb(self):
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
        if result.exit_code != 0:
            print(result.exception)
        self.assertEqual(result.exit_code, 0)

    def test_generate_decoy_database(self):
        """
        Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
        :return:
        """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['generate-decoy', '--config_file', 'config/protein_decoy.yaml',
                                '-in', 'testdata/test_db.fa', '-out', 'testdata/output_decoy.fa', '--method',
                                'protein-reverse'])
        self.assertEqual(result.exit_code, 0)

    def test_generate_decoy_database_noconfig(self):
        """
              Test generation proteinDB from altORFs using dnaseq-to-proteindb tool
              :return:
              """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['generate-decoy', '-in', 'testdata/test_db.fa', '-out', 'testdata/output_decoy.fa',
                                '--method',
                                'protein-reverse'])
        self.assertEqual(result.exit_code, 0)

    def test_download_ensembl_data(self):
        """
              Test downloading ensembl data for species with taxonomy identifer = 9103
              :return:
              """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['ensembl-downloader', '--config_file', 'config/ensembl_downloader_config.yaml',
                                '--skip_dna',
                                '--taxonomy', '9103', '--output_directory', 'testdata'])
        self.assertEqual(result.exit_code, 0)

    def test_download_ensembl_data_37(self):
        """
              Test downloading ensembl data for species with taxonomy identifier = 9606
              :return:
              """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['ensembl-downloader', '--taxonomy', '9103', '--output_directory', 'testdata',
                                '--grch37'])
        self.assertEqual(result.exit_code, 0)

    def test_download_cbioportal_data(self):
        """
              Test downloading cbioportal data for study id: paac_jhu_2014
              :return:
              """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['cbioportal-downloader', '--config_file', 'config/cbioportal_config.yaml',
                                '--download_study', 'paac_jhu_2014', '--output_directory', 'testdata'])
        self.assertEqual(result.exit_code, 0)

    def test_download_cbioportal_data_noconfig(self):
        """
              Test downloading cbioportal data for study id: paac_jhu_2014
              :return:
              """
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['cbioportal-downloader', '--download_study', 'paac_jhu_2014', '--output_directory',
                                'testdata'])
        self.assertEqual(result.exit_code, 0)

    def test_check_ensembl_database(self):
        runner = CliRunner()
        result = runner.invoke(cli,
                               ['ensembl-check', '--config_file', 'config/ensembl_config.yaml',
                                '--input_fasta', 'testdata/proteindb_from_ENSEMBL_VCF.fa', '--output',
                                'testdata/proteindb_from_ENSEMBL_VCF-clean.fa', '--add_stop_codons', '--num_aa', '6'])
        self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()
