package org.bigbio.pgatk.pepgenome.gui;

import com.sun.org.apache.xpath.internal.operations.Bool;
import javafx.concurrent.Task;
import org.apache.log4j.Logger;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.PepGenomeTool;
import org.bigbio.pgatk.pepgenome.common.Assembly;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.io.GTFParser;
import org.bigbio.pgatk.pepgenome.io.GenomeFastaParser;
import org.bigbio.pgatk.pepgenome.io.MzTabInputPeptideFileParser;
import org.bigbio.pgatk.pepgenome.io.TabInputPeptideFileParser;
import org.bigbio.pgatk.pepgenome.io.custom.PeptideAtlasPeptideParser;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import org.bigbio.pgatk.pepgenome.kmer.inmemory.KmerSortedMap;
import org.bigbio.pgatk.pepgenome.kmer.inmemory.KmerTreeMap;

public class PepGenomeTask extends Task {

    private static final org.apache.log4j.Logger log = Logger.getLogger(PepGenomeTask.class);

    private String fileStr;
    private String fastaStr;
    private String gtfStr;
    private Boolean mmMode;
    private int numberMisMatches;
    private final Boolean chrIncluded;
    private Boolean mergeOuput;
    private Boolean gtfOutFlag;
    private Boolean gctOutFlag;
    private Boolean bedOutFlag;
    private Boolean ptmbedOutFlag;


    public PepGenomeTask(String fileStr, String fastaStr, String gtfStr,
                         Boolean mmMode, int numberMisMatches,
                         Boolean chrIncluded,
                         Boolean mergeOuput,
                         Boolean gtfOutFlag, Boolean bedOutFlag, Boolean gctOutFlag, Boolean ptmBedOutFlag) {
        this.fileStr = fileStr;
        this.fastaStr = fastaStr;
        this.gtfStr = gtfStr;
        this.mmMode = mmMode;
        this.numberMisMatches = numberMisMatches;
        this.chrIncluded = chrIncluded;
        this.mergeOuput = mergeOuput;
        this.gtfOutFlag = gtfOutFlag;
        this.gctOutFlag = gctOutFlag;
        this.bedOutFlag = bedOutFlag;
        this.ptmbedOutFlag = ptmBedOutFlag;
    }

    @Override
    protected Object call() throws Exception {

        String[] peptideInputFilePaths = Utils.tokenize(this.fileStr, ",", true);


        boolean inMemory = true;

        if (mergeOuput && peptideInputFilePaths.length == 1) {
            System.err.println("cannot merge output files for one input file, default (-merge false) assumed");
            mergeOuput = false;
        }

        String pluralString = GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES == 1 ? " mismatch" : " mismatches";
        log.info("Start: allowing " + GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES + pluralString);

        try {
            if (fastaStr != null) {
                log.info("reading genome FASTA: " + fastaStr);
                GenomeFastaParser.readGenomeFASTA(fastaStr);
            }

            log.info("reading FASTA: " + fastaStr);
            CoordinateWrapper coordinate_wrapper = new CoordinateWrapper();
            coordinate_wrapper.read_fasta_file(fastaStr);
            log.info("Fasta done: " + coordinate_wrapper.size() + " proteins read.");
            log.info("building KmerTreeMap...");
            int kmerSize = coordinate_wrapper.getTotalAACount() / coordinate_wrapper.getNumberOfProteins();
            kmerSize = kmerSize / GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH * coordinate_wrapper.getNumberOfProteins();
            Object kmer_map;
            if (inMemory) {
                kmer_map = new KmerTreeMap();
            } else {
                kmer_map = new KmerSortedMap(kmerSize);
            }

            coordinate_wrapper.add_all_proteins_to_kmer_map((IKmerMap)kmer_map);
            if (!inMemory) {
                ((KmerSortedMap)kmer_map).sortKmer();
            }

            MappedPeptides mapped_peptides = new MappedPeptides();
            Assembly assem = GTFParser.get_instance().read(gtfStr, coordinate_wrapper, mapped_peptides);

            String filename_mm_postfix = "";
            if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES > 0) {
                StringBuilder ss = new StringBuilder();
                ss.append("_").append(GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES).append("MM");
                if (GenomeMapper.PEPTIDE_MAPPER.ONE_IN_FIVE_MODE) {
                    ss.append("-(1in5)");
                }

                filename_mm_postfix = ss.toString();
            }

            String[] var52 = peptideInputFilePaths;
            int var27 = peptideInputFilePaths.length;

            String source = "PoGo";

            String path7;
            String path8;
            String path81;
            String path4;
            String path5;

            PepGenomeTool.INPUT_FILE_FORMAT fileFormat = PepGenomeTool.INPUT_FILE_FORMAT.MZTAB;

            for(int var28 = 0; var28 < var27; ++var28) {
                path7 = var52[var28];
                log.info(path7);
                path8 = path7;
                if (path7.endsWith(".txt")) {
                    path8 = Utils.removeExtension(path7, ".txt");
                } else if (path7.endsWith(".tsv")) {
                    path8 = Utils.removeExtension(path7, ".tsv");
                } else if (path7.endsWith(".pogo")) {
                    path8 = Utils.removeExtension(path7, ".pogo");
                } else if (path7.endsWith(".mztab")) {
                    path8 = Utils.removeExtension(path7, ".mztab");
                } else if (path7.endsWith(".mzid")) {
                    path8 = Utils.removeExtension(path7, ".mzid");
                }

                path81 = path8 + "_unmapped.txt";
                if (fileFormat == PepGenomeTool.INPUT_FILE_FORMAT.MZTAB) {
                    (new MzTabInputPeptideFileParser()).read(path7, coordinate_wrapper, mapped_peptides, path81, (IKmerMap)kmer_map);
                } else if (fileFormat == PepGenomeTool.INPUT_FILE_FORMAT.PEPTIDEATLAS) {
                    (new PeptideAtlasPeptideParser()).read(path7, coordinate_wrapper, mapped_peptides, path81, (IKmerMap)kmer_map);
                } else {
                    (new TabInputPeptideFileParser()).read(path7, coordinate_wrapper, mapped_peptides, path81, (IKmerMap)kmer_map);
                }

                log.info("Results done! (" + path7 + ")");
                log.info("writing output files");
                if (!mergeOuput) {
                    path4 = path8 + filename_mm_postfix + "_out.gtf";
                    path5 = path8 + filename_mm_postfix + ".bed";
                    path7 = path8 + filename_mm_postfix + ".gct";
                    path8 = path8 + filename_mm_postfix + "_ptm.bed";
                    path81 = path8 + filename_mm_postfix + "_no-ptm.bed";
                    String path9 = path8 + filename_mm_postfix + "_" + assem.toString() + "_out.gtf";
                    String path10 = path8 + filename_mm_postfix + "_" + assem.toString() + ".bed";
                    String path11 = path8 + filename_mm_postfix + "_" + assem.toString() + ".gct";
                    String path12 = path8 + filename_mm_postfix + "_" + assem.toString() + "_ptm.bed";
                    String path121 = path8 + filename_mm_postfix + "_" + assem.toString() + "_no-ptm.bed";
                    if (assem == Assembly.patchhaploscaff) {
                        path9 = path8 + filename_mm_postfix + "_patch_hapl_scaff_out.gtf";
                        path10 = path8 + filename_mm_postfix + "_patch_hapl_scaff.bed";
                        path11 = path8 + filename_mm_postfix + "_patch_hapl_scaff.gct";
                        path12 = path8 + filename_mm_postfix + "_patch_hapl_scaff_ptm.bed";
                        path121 = path8 + filename_mm_postfix + "_patch_hapl_scaff_no-ptm.bed";
                    }

                    if (gtfOutFlag) {
                        mapped_peptides.to_gtf(path4, source);
                        mapped_peptides.to_gtf(path9, source, assem);
                    }

                    if (bedOutFlag) {
                        mapped_peptides.to_bed(path5, Assembly.primary, this.chrIncluded);
                        mapped_peptides.to_bed(path10, assem, this.chrIncluded);
                    }

                    if (gctOutFlag) {
                        mapped_peptides.to_gct(path7, Assembly.primary, this.chrIncluded);
                        mapped_peptides.to_gct(path11, assem);
                    }

                    if (ptmbedOutFlag) {
                        mapped_peptides.to_ptmbed(path8, path81);
                        mapped_peptides.to_ptmbed(path12, path121, assem);
                    }

                    mapped_peptides.remove_all_peptides();
                }
            }

            if (mergeOuput) {
                String final_peptide_path_results = peptideInputFilePaths[0];
                if (final_peptide_path_results.endsWith(".txt")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".txt");
                } else if (final_peptide_path_results.endsWith(".tsv")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".tsv");
                } else if (final_peptide_path_results.endsWith(".pogo")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".pogo");
                } else if (final_peptide_path_results.endsWith(".mztab")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".mztab");
                } else if (final_peptide_path_results.endsWith(".mzid")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".mzid");
                }

                path4 = final_peptide_path_results + filename_mm_postfix + "_merged.gtf";
                path5 = final_peptide_path_results + filename_mm_postfix + "_merged.bed";
                path7 = final_peptide_path_results + filename_mm_postfix + "_merged.gct";
                path8 = final_peptide_path_results + filename_mm_postfix + "_merged_ptm.bed";
                path81 = final_peptide_path_results + filename_mm_postfix + "_merged_no-ptm.bed";
                path4 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.gtf";
                path5 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.bed";
                path7 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.gct";
                path8 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged_ptm.bed";
                path81 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged_no-ptm.bed";
                if (assem == Assembly.patchhaploscaff) {
                    path4 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.gtf";
                    path5 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.bed";
                    path7 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.gct";
                    path8 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged_ptm.bed";
                    path81 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged_no-ptm.bed";
                }

                if (gtfOutFlag) {
                    mapped_peptides.to_gtf(path4, source);
                    mapped_peptides.to_gtf(path4, source, assem);
                }

                if (bedOutFlag) {
                    mapped_peptides.to_bed(path5);
                    mapped_peptides.to_bed(path5, assem);
                }

                if (gctOutFlag) {
                    mapped_peptides.to_gct(path7);
                    mapped_peptides.to_gct(path7, assem);
                }

                if (ptmbedOutFlag) {
                    mapped_peptides.to_ptmbed(path8, path81);
                    mapped_peptides.to_ptmbed(path8, path81, assem);
                }
            }
        } catch (Exception var45) {
            log.info(var45.getMessage());
            System.err.println(var45.getMessage());
            var45.printStackTrace();
        }

        log.info("DONE..");
        return null;
    }

    @Override
    protected void succeeded() {
        super.succeeded();
    }
}
