import re

from Bio import SeqIO

from pypgatk.cgenomes.models import SNP
from pypgatk.toolbox.general import ParameterConfiguration


class CancerGenomesService(ParameterConfiguration):
    CONFIG_CANCER_GENOMES_MUTATION_FILE = 'mutation_file'
    CONFIG_COMPLETE_GENES_FILE = "genes_file"
    CONFIG_OUTPUT_FILE = "output_file"
    CONFIG_COSMIC_DATA = "cosmic_data"
    CONFIG_KEY_DATA = 'proteindb'
    CONFIG_TISSUE_INFO = 'tissue_info'
    TISSUE_TYPE = "tissue_type"
    SPLIT_BY_TISSUE_TYPE = "split_by_tissue_type"
    CLINICAL_SAMPLE_FILE = 'clinical_sample_file'
    
    def __init__(self, config_file, pipeline_arguments):
        """
        Init the class with the specific parameters.
        :param config_file configuration file
        :param pipeline_arguments pipelines arguments
        """
        super(CancerGenomesService, self).__init__(self.CONFIG_COSMIC_DATA, config_file, pipeline_arguments)

        if self.CONFIG_CANCER_GENOMES_MUTATION_FILE in self.get_pipeline_parameters():
            self._local_mutation_file = self.get_pipeline_parameters()[self.CONFIG_CANCER_GENOMES_MUTATION_FILE]

        if self.CONFIG_COMPLETE_GENES_FILE in self.get_pipeline_parameters():
            self._local_complete_genes = self.get_pipeline_parameters()[self.CONFIG_COMPLETE_GENES_FILE]

        if self.CONFIG_OUTPUT_FILE in self.get_pipeline_parameters():
            self._local_output_file = self.get_pipeline_parameters()[self.CONFIG_OUTPUT_FILE]
        
        self._tissue_type = self.get_multiple_options(
            self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_TISSUE_INFO][self.TISSUE_TYPE])
        if self.TISSUE_TYPE in self.get_pipeline_parameters():
            self._tissue_type = self.get_multiple_options(self.get_pipeline_parameters()[self.TISSUE_TYPE])

        self._split_by_tissue_type = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_TISSUE_INFO][self.SPLIT_BY_TISSUE_TYPE]
        if self.SPLIT_BY_TISSUE_TYPE in self.get_pipeline_parameters():
            self._split_by_tissue_type = self.get_pipeline_parameters()[self.SPLIT_BY_TISSUE_TYPE]
        
        self._local_clinical_sample_file = self.get_default_parameters()[self.CONFIG_KEY_DATA][self.CONFIG_TISSUE_INFO][self.CLINICAL_SAMPLE_FILE]
        if self.CLINICAL_SAMPLE_FILE in self.get_pipeline_parameters():
            self._local_clinical_sample_file = self.get_pipeline_parameters()[self.CLINICAL_SAMPLE_FILE]
    
    @staticmethod
    def get_multiple_options(options_str: str):
        """
        This method takes an String like option1, option2, ... and produce and array [option1, option2,... ]
        :param options_str:
        :return: Array
        """
        return list(map(lambda x: x.strip(), options_str.split(",")))

    def cosmic_to_proteindb(self):
        """
        This function translate the mutation file + COSMIC genes into a protein Fasta database. The
        method write into the file system the output Fasta.
        :return:
        """
        self.get_logger().debug("Starting reading the All cosmic genes")
        COSMIC_CDS_DB = SeqIO.index(self._local_complete_genes, 'fasta')  # All_COSMIC_Genes.fasta

        cosmic_input = open(self._local_mutation_file, 'r')  # CosmicMutantExport.tsv

        header = cosmic_input.readline().split("\t")
        regex = re.compile('[^a-zA-Z]')
        gene_col = header.index("Gene name")
        enst_col = header.index("Accession Number")
        cds_col = header.index("Mutation CDS")
        aa_col = header.index("Mutation AA")
        muttype_col = header.index("Mutation Description")
        tissue_col = header.index('Primary site')
        
        output = open(self._local_output_file, 'w')

        mutation_dic = {}
        tissue_mutations_dict = {}
        nucleotide = ["A", "T", "C", "G"]
        self.get_logger().debug("Reading input CosmicMutantExport.tsv ...")
        line_counter = 1
        for line in cosmic_input:
            if line_counter % 10000 == 0:
                msg = "Number of lines finished -- '{}'".format(line_counter)
                self.get_logger().debug(msg)
            line_counter += 1
            row = line.strip().split("\t")
            #filter out mutations from unspecified tissues
            if row[tissue_col] not in self._tissue_type and self._tissue_type!=['all']:
                continue
            
            if "coding silent" in row[muttype_col]:
                continue

            snp = SNP(gene=row[gene_col], mrna=row[enst_col], dna_mut=row[cds_col], aa_mut=row[aa_col],
                      type=row[muttype_col])
            header = "COSMIC:%s:%s:%s" % (snp.gene, snp.aa_mut, snp.type.replace(" ", ""))
            try:
                seq = COSMIC_CDS_DB[snp.gene].seq
            except KeyError:  # geneID is not in All_COSMIC_Genes.fasta
                continue

            mut_pro_seq = ""
            if "?" not in row[cds_col]:  # unambiguous DNA change known in CDS sequence
                positions = re.findall(r'\d+', snp.dna_mut)
                if ">" in snp.dna_mut and len(positions) == 1:  # Substitution
                    tmplist = snp.dna_mut.split(">")
                    ref_dna = re.sub("[^A-Z]+", "", tmplist[0])
                    mut_dna = re.sub("[^A-Z]+", "", tmplist[1])
                    index = int(positions[0]) - 1

                    if ref_dna == str(seq[index]).upper() and mut_dna in nucleotide:  #
                        seq_mut = seq[:index] + mut_dna + seq[index + 1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
                elif "ins" in snp.dna_mut:
                    index = snp.dna_mut.index("ins")
                    insert_dna = snp.dna_mut[index + 3:]
                    if insert_dna.isalpha():
                        #print(insert_dna, snp.dna_mut, snp.aa_mut)
                        ins_index1 = int(positions[0])
                        seq_mut = seq[:ins_index1] + insert_dna + seq[ins_index1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)

                elif "del" in snp.dna_mut:
                    if len(positions) == 2:
                        del_index1 = int(positions[0]) - 1
                        del_index2 = int(positions[1])
                        seq_mut = seq[:del_index1] + seq[del_index2:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
                    elif len(positions) == 1:
                        del_index1 = int(positions[0]) - 1
                        seq_mut = seq[:del_index1] + seq[del_index1 + 1:]
                        mut_pro_seq = seq_mut.translate(to_stop=True)
            else:
                if "?" not in row[aa_col]:  # unambiguous aa change known in protein sequence
                    positions = re.findall(r'\d+', snp.aa_mut)
                    protein_seq = str(seq.translate(to_stop=True))

                    if "Missense" in snp.type:
                        mut_aa = snp.aa_mut[-1]
                        index = int(positions[0]) - 1
                        mut_pro_seq = protein_seq[:index] + mut_aa + protein_seq[index + 1:]
                    elif "Nonsense" in snp.type:
                        index = int(positions[0]) - 1
                        mut_pro_seq = protein_seq[:index]
                    elif "Insertion - In frame" in snp.type:
                        index = snp.aa_mut.index("ins")
                        insert_aa = snp.aa_mut[index + 3:]
                        if insert_aa.isalpha():
                            ins_index1 = int(positions[0])
                            mut_pro_seq = protein_seq[:ins_index1] + insert_aa + protein_seq[ins_index1:]
                    elif "Deletion - In frame" in snp.type:
                        if len(positions) == 2:
                            del_index1 = int(positions[0]) - 1
                            del_index2 = int(positions[1])
                            mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index2:]
                        elif len(positions) == 1:
                            del_index1 = int(positions[0]) - 1
                            mut_pro_seq = protein_seq[:del_index1] + protein_seq[del_index1 + 1:]
                    elif "Complex" in snp.type and "frameshift" not in snp.type:
                        try:
                            index = snp.aa_mut.index(">")
                        except ValueError:
                            # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                            continue
                        mut_aa = snp.aa_mut[index + 1:]
                        if "deletion" in snp.type:
                            try:
                                del_index1 = int(positions[0]) - 1
                                del_index2 = int(positions[1])
                                mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                            except IndexError:
                                # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                continue
                        elif "insertion" in snp.type:
                            try:
                                ins_index1 = int(positions[0]) - 1
                            except IndexError:
                                # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                continue
                            mut_pro_seq = protein_seq[:ins_index1] + mut_aa + protein_seq[ins_index1 + 1:]
                        elif "compound substitution" in snp.type:
                            if "*" not in mut_aa:
                                try:
                                    del_index1 = int(positions[0]) - 1
                                    del_index2 = int(positions[1])
                                    mut_pro_seq = protein_seq[:del_index1] + mut_aa + protein_seq[del_index2:]
                                except IndexError:
                                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                    continue
                            else:
                                try:
                                    del_index1 = int(positions[0]) - 1
                                    mut_pro_seq = protein_seq[:del_index1] + mut_aa.replace("*", "")
                                except IndexError:
                                    # print (snp.gene,snp.mRNA,snp.dna_mut,snp.aa_mut,snp.type)
                                    continue
            if mut_pro_seq != "":
                entry = ">%s\n%s\n" % (header, mut_pro_seq)
                if header not in mutation_dic:
                    output.write(entry)
                    mutation_dic[header] = 1
                
                if self._split_by_tissue_type:
                    try:
                        tissue_mutations_dict[row[tissue_col]][header] = entry
                    except KeyError:
                        tissue_mutations_dict[row[tissue_col]] = {header: entry}
            
        for tissue_type in tissue_mutations_dict.keys():
            with open(self._local_output_file.replace('.fa', '')+ '_' + regex.sub('', tissue_type) +'.fa', 'w') as fn:
                for header in tissue_mutations_dict[tissue_type].keys():
                    fn.write(tissue_mutations_dict[tissue_type][header])
            
        self.get_logger().debug("COSMIC contains in total {} non redundant mutations".format(len(mutation_dic)))
        cosmic_input.close()
        output.close()
        
    @staticmethod
    def get_sample_headers(header_line):
        try:
            tissue_type_col = header_line.index('CANCER_TYPE')
        except ValueError:
            print('CANCER_TYPE was not found in the header row:', header_line)
            return None, None
        try:
            sample_id_col = header_line.index('SAMPLE_ID')
        except ValueError:
            print('SAMPLE_ID was not found in the header row:', header_line)
            return None, None
        return tissue_type_col, sample_id_col
    
    def get_tissue_type_per_sample(self, local_clinical_sample_file):
        sample_tissue_type = {}
        if local_clinical_sample_file:
            with open(local_clinical_sample_file, 'r') as clin_fn:
                tissue_type_col, sample_id_col = None, None
                for line in clin_fn.readlines():
                    if line.startswith('#'):
                        continue
                    sl = line.strip().split('\t')
                    #check for header and re-assign columns
                    if 'SAMPLE_ID' in sl and 'CANCER_TYPE' in sl:
                        tissue_type_col, sample_id_col = self.get_sample_headers(sl)
                    if tissue_type_col and sample_id_col:
                        sample_tissue_type[sl[sample_id_col]] = sl[tissue_type_col].strip().replace(' ','_')
        return sample_tissue_type
    
    @staticmethod
    def get_mut_header_cols(header_cols, row, tissue_type, split_by_tissue_type):
        for col in header_cols.keys():
            header_cols[col] = row.index(col)
        
        #check if tissue type should be considered
        if tissue_type!=['all'] or split_by_tissue_type:
            try:
                header_cols["Tumor_Sample_Barcode"] = row.index("Tumor_Sample_Barcode")
            except ValueError:
                print("Tumor_Sample_Barcode was not found in the header {} of mutations file".format(row))
                header_cols["Tumor_Sample_Barcode"] = None
        
        return header_cols
    
    def cbioportal_to_proteindb(self):
        """cBioportal studies have a data_clinical_sample.txt file that shows the Primary Tumor Site per Sample Identifie
        it matches the the Tumor_Sample_Barcode column in the mutations file.
        """
        regex = re.compile('[^a-zA-Z]')
        mutfile = open(self._local_mutation_file, "r")
        fafile = SeqIO.parse(self._local_complete_genes, "fasta")
        output = open(self._local_output_file, "w")
        sample_tissue_type_dict = {}
        tissue_mutations_dict = {}
        
        seq_dic = {}
        for record in fafile:
            newacc = record.id.split(".")[0]
            if newacc not in seq_dic:
                seq_dic[newacc] = record.seq
        
        header_cols = {"HGVSc": None, "Transcript_ID": None, "Variant_Classification": None, 
                       "Variant_Type": None, "HGVSp_Short": None}
        nucleotide = ["A", "T", "C", "G"]
        mutclass = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                    "Nonsense_Mutation"]
        
        # check if sample id and clinical files are given, if not and tissue type is required then exit
        if self._tissue_type!=['all'] or self._split_by_tissue_type:
            if self._local_clinical_sample_file:
                sample_tissue_type_dict = self.get_tissue_type_per_sample(self._local_clinical_sample_file)
                if sample_tissue_type_dict=={}:
                    return
            else:
                print('No clinical sample file is given therefore no tissue type can be detected.')
                return
            
        for i,line in enumerate(mutfile):
            row = line.strip().split("\t")
            if row[0]=='#':
                print("skipping line ({}): {}".format(i, row))
                continue
            #check for header in the mutations file and get column indecis
            if set(header_cols.keys()).issubset(set(row)):
                header_cols = self.get_mut_header_cols(header_cols, row, self._tissue_type, self._split_by_tissue_type)
            #check if any is none in header_cols then continue
            if None in header_cols.values():
                print("Incorrect header column is given")
                continue
            #get tissue type and check it
            tissue_type = None
            if self._tissue_type!=['all'] or self._split_by_tissue_type:
                try:
                    tissue_type = sample_tissue_type_dict[row[header_cols["Tumor_Sample_Barcode"]]]
                except KeyError:
                    if self._tissue_type!=['all'] or self._split_by_tissue_type:
                        print("No clinical info was found for sample {}. Skipping (line {}): {}".format(row[header_cols["Tumor_Sample_Barcode"]], i, line))
                        continue
                except IndexError:
                    print("No sampleID was found in (line {}): {}".format(i, row))
            if tissue_type not in self._tissue_type and self._tissue_type != ['all']:
                continue
            
            gene = row[0]
            try:
                pos = row[header_cols["HGVSc"]]
                enst = row[header_cols["Transcript_ID"]]

                seq_mut = ""
                aa_mut = row[header_cols["HGVSp_Short"]]

                vartype = row[header_cols["Variant_Type"]]
                varclass = row[header_cols["Variant_Classification"]]
            except IndexError:
                print("Incorrect line (i):", row)
                continue
            if varclass not in mutclass:
                continue
            
            try:
                seq = seq_dic[enst]
            except KeyError:
                print("%s not found:" % enst)
                continue
            
            if ":" in pos:
                cdna_pos = pos.split(":")[1]
            else:
                cdna_pos = pos

            if vartype == "SNP":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos)[0])
                except IndexError:
                    print("Incorrect SNP format or record", i, pos, line)
                    continue
                idx = pos.index(">")
                ref_dna = pos[idx - 1]
                mut_dna = pos[idx + 1]

                if mut_dna not in nucleotide:
                    print(mut_dna, "is not a nucleotide base", pos)
                    continue
                try:
                    if ref_dna == seq[enst_pos - 1]:
                        seq_mut = seq[:enst_pos - 1] + mut_dna + seq[enst_pos:]
                    else:
                        print("incorrect substitution, unmatched nucleotide", pos, enst)
                except IndexError:
                    print("incorrect substitution, out of index", pos)
            elif vartype == "DEL":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                except IndexError:
                    print("incorrect del format or record", i, pos, line)
                    continue
                del_dna = pos.split("del")[1]
                if del_dna == seq[enst_pos - 1:enst_pos - 1 + len(del_dna)]:
                    seq_mut = seq[:enst_pos - 1] + seq[enst_pos - 1 + len(del_dna):]
                else:
                    print("incorrect deletion, unmatched nucleotide", pos)

            elif vartype == "INS":
                try:
                    enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[0])[0])
                except IndexError:
                    print("incorrect ins/dup format or record", i, pos, line)
                    continue
                if "ins" in pos:
                    ins_dna = pos.split("ins")[1]
                elif "dup" in pos:
                    ins_dna = pos.split("dup")[1]
                    if len(ins_dna) > 1:
                        enst_pos = int(re.findall(r'\d+', cdna_pos.split("_")[1])[0])
                else:
                    print("unexpected insertion format")
                    continue;

                seq_mut = seq[:enst_pos] + ins_dna + seq[enst_pos:]

            if seq_mut == "":
                continue;

            mut_pro_seq = seq_mut.translate(to_stop=True)
            if len(mut_pro_seq) > 6:
                header = "cbiomut:%s:%s:%s:%s" % (enst, gene, aa_mut, varclass)
                output.write(">%s\n%s\n" % (header, mut_pro_seq))
                
                if self._split_by_tissue_type:
                    try:
                        tissue_mutations_dict[tissue_type][header] = mut_pro_seq
                    except KeyError:
                        tissue_mutations_dict[tissue_type] = {header: mut_pro_seq}
                    
        output.close()
        mutfile.close()
        fafile.close()
        
        for tissue_type in tissue_mutations_dict.keys():
            with open(self._local_output_file.replace('.fa', '')+ '_' + regex.sub('', tissue_type) +'.fa', 'w') as fn:
                for header in tissue_mutations_dict[tissue_type].keys():
                    fn.write(">{}\n{}\n".format(header, tissue_mutations_dict[tissue_type][header]))
         
