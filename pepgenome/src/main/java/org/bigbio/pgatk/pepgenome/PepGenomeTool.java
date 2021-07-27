package org.bigbio.pgatk.pepgenome;


import lombok.extern.slf4j.Slf4j;
import org.apache.commons.cli.*;
import org.apache.commons.io.FilenameUtils;
import org.bigbio.pgatk.pepgenome.common.Assembly;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.io.*;
import org.bigbio.pgatk.pepgenome.io.custom.PeptideAtlasPeptideParser;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import org.bigbio.pgatk.pepgenome.kmer.inmemory.KmerSortedMap;
import org.bigbio.pgatk.pepgenome.kmer.inmemory.KmerTreeMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.stream.Stream;

@Slf4j
public class PepGenomeTool {

    private static Logger log = LoggerFactory.getLogger(PepGenomeTool.class);

    public enum INPUT_FILE_FORMAT {
        TAB("tab", "Tab delimited input (.pogo, .tsv, .txt)", TabInputPeptideFileParser.class),
        MZTAB("mztab", "mzTab file format (.mztab)", MzTabInputPeptideFileParser.class),
        PRIDE("pavro", "Peptide binary avro file", PeptideAvroFileParser.class),
        PEPTIDEATLAS("peptideatlas", "PeptideAtlas PeptideBuild (.tsv)", PeptideAtlasPeptideParser.class),
        MZIDENML("mzid", "MzIndetML file format (.mzid)", MzTabInputPeptideFileParser.class);

        private String name;
        private String description;
        private String Class;

        INPUT_FILE_FORMAT(String name, String description, Class classType) {
            this.name = name;
            this.description = description;
        }

        public String getName() {
            return name;
        }

        public String getDescription() {
            return description;
        }

        public static INPUT_FILE_FORMAT findByString(String key) {
            for (INPUT_FILE_FORMAT value : values())
                if (value.getName().equalsIgnoreCase(key))
                    return value;
            return INPUT_FILE_FORMAT.TAB;
        }
    }


    //exit codes --------------------------------
    private static final int GENOME_MAPPER_EXIT_HELP = 1;
    private static final int GENOME_MAPPER_EXIT_TOO_FEW_ARGS = 2;
    private static final int GENOME_MAPPER_EXIT_INVALID_ARG = 3;

    //-----------------input args--------------------------
    private static final String ARG_FASTA = "fasta";
    private static final String ARG_GTF = "gtf";
    private static final String ARG_GFF = "gff"; // Addition - GFF3 (only) accepted
    private static final String ARG_ANN = "ann"; // Addition - Annotation, GTF or GFF accepted
    private static final String ARG_EXON_COORDS = "exco"; // Addition - Use genomic coordinates of exons rather than CDS annotation features.  (Peptides unannotated).
    private static final String ARG_PEPTIDE_FILTER = "variant_filter"; // Addition - Peptide Filter to mark variants.
    private static final String ARG_IN = "in";
    private static final String ARG_MERGE = "merge";
    private static final String ARG_FORMAT = "format";
    private static final String ARG_SOURCE = "source";
    private static final String ARG_MM = "mm";
    private static final String ARG_MMMODE = "mmmode";
    private static final String ARG_CHR = "chr";
    private static final String ARG_HELP = "h";
    private static final String ARG_INMEMORY = "inm";
    private static final String ARG_INPUT_FORMAT = "inf";
    private static final String ARG_GENOME_FASTA = "genome";

    //DEFAULT values
    private static boolean mergeFlag = false;
    private static boolean gtfOutFlag = true;
    private static boolean gctOutFlag = true;
    private static boolean bedOutFlag = true;
    private static boolean ptmbedOutFlag = true;
    private static String source = "PoGo";
    private static boolean chrincluded = false;
    private static boolean inMemory = true;
    private static INPUT_FILE_FORMAT fileFormat = INPUT_FILE_FORMAT.TAB;
    public static boolean useExonCoords = false; // ||Exco mode||
    public static boolean usePeptideFilter = false; // ||Peptide filter mode||

    // Map of transcript IDs to CDS offsets for non-CDS-annotated use.
    public static HashMap<String, Integer> m_translation_offset_map = new HashMap<>();

    public static void main(String[] args) {

        System.setProperty("hadoop.home.dir", "/");

        long startTime = System.nanoTime();

        Options options = new Options();
        options.addOption(Option.builder(ARG_FASTA).hasArg(true)
                .desc("Filepath for file containing protein sequences in FASTA format").build())
                .addOption(Option.builder(ARG_GTF).hasArg(true)
                        .desc("Filepath for file containing genome annotation in GTF format").build())
                .addOption(Option.builder(ARG_GFF).hasArg(true).desc("Filepath for file containing genome annotation in GFF3 format").build())
                .addOption(Option.builder(ARG_ANN).hasArg(true).desc("Filepath for file containing genome annotation in GTF or GFF3 format").build())
                .addOption(Option.builder(ARG_EXON_COORDS).hasArg(true).desc("Use exon coordinates rather than CDS (Unannotated peptides)").build())
                .addOption(Option.builder(ARG_PEPTIDE_FILTER).hasArg(true).desc("Peptide filter mode.").build())
                .addOption(Option.builder(ARG_IN).hasArg(true)
                        .desc("Comma(,) separated file paths for files containing peptide identifications " +
                                "(Contents of the file can tab separated format. i.e., " +
                                "File format: four columns: SampleName\t\tPeptideSequence\t\tPSMs\tQuant; " +
                                "or mzTab, and mzIdentML)").build())
                .addOption(Option.builder(ARG_MERGE).hasArg(true)
                        .desc("Set 'true' to merge mappings from all files from input (default 'false')").build())
                .addOption(Option.builder(ARG_FORMAT).hasArg(true)
                        .desc("Select the output formats from gtf, gct, bed, ptmbed, all or combinations" +
                                " thereof separated by ',' (default all)").build())
                .addOption(Option.builder(ARG_SOURCE).hasArg(true)
                        .desc("Please give a source name which will be used in the second column" +
                                " in the output gtf file (default: PoGo)").build())
                .addOption(Option.builder(ARG_MM).hasArg(true)
                        .desc("Allowed mismatches (0, 1 or 2; default: 0)").build())
                .addOption(Option.builder(ARG_MMMODE).hasArg(true)
                        .desc("Mismatch mode (true or false): if true mismatching with two mismatches will" +
                                " only allow 1 mismatch every kmersize (default: 5) positions." +
                                " (default: false)").build())
                .addOption(Option.builder(ARG_GENOME_FASTA).hasArg(true)
                        .desc("Filepath for file containing genome sequence in FASTA format " +
                                "used to extract chromosome names and order and differenciate between" +
                                " assembly and scaffolds. If not set chromosome and scaffold names and" +
                                " order is extracted from GTF input.").build())
                .addOption(Option.builder(ARG_CHR).hasArg(true)
                        .desc("Export chr prefix Allowed 0, 1  (default: 0)").build())
                .addOption(Option.builder(ARG_INMEMORY).hasArg(true)
                        .desc("Compute the kmer algorithm in memory or using database" +
                                " algorithm (default 0, database 1)").build())
                .addOption(Option.builder(ARG_INPUT_FORMAT).hasArg(true)
                        .desc("Format of the input file (mztab, mzid, pavro, or tsv). (default tsv) ").build())
                .addOption(Option.builder(ARG_HELP).hasArg(false)
                        .desc("Print this help & exit").build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
            if (cmd.hasOption(ARG_HELP)) {
                Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_HELP);
            }
        } catch (ParseException e) {
            log.info(" *** Error in parsing the input arguments. Please check the arguments ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_TOO_FEW_ARGS);
        }

        // Check for protein sequence input and peptide input.
        if (!cmd.hasOption(ARG_FASTA) || !cmd.hasOption(ARG_IN)) {
            log.info("*** Missing mandatory parameters: -fasta and/or -in ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_TOO_FEW_ARGS);
        }

        // Check for genome annotation input.  Can be ANN (GTF/GFF3), GTF or GFF3.
        if (!cmd.hasOption(ARG_GTF) && !cmd.hasOption(ARG_ANN) && !cmd.hasOption(ARG_GFF)) {
            log.info("*** Missing mandatory parameters: -gtf/-gff3/-ann ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_TOO_FEW_ARGS);
        }

        if (cmd.hasOption(ARG_INPUT_FORMAT))
            fileFormat = INPUT_FILE_FORMAT.findByString(cmd.getOptionValue(ARG_INPUT_FORMAT));

        if (cmd.hasOption(ARG_INMEMORY) && cmd.getOptionValue(ARG_INMEMORY).equalsIgnoreCase("1")) {
            inMemory = false;
        }

        // Exon coords argument check
        if (cmd.hasOption(ARG_EXON_COORDS)) {
            useExonCoords = true;
        }

        // Peptide filter argument check
        if (cmd.hasOption(ARG_PEPTIDE_FILTER)) {
            usePeptideFilter = true;
        }

        String fastaFilePath = cmd.getOptionValue(ARG_FASTA);
        String peptideInputFilePathsParam = cmd.getOptionValue(ARG_IN);
        String fastaGenomeFilePath = cmd.getOptionValue(ARG_GENOME_FASTA);

        String annFilePath = "";
        if (cmd.getOptionValue(ARG_GFF) != null && cmd.getOptionValue(ARG_GFF).endsWith(".gff3")) {
            annFilePath = cmd.getOptionValue(ARG_GFF);

        }
        else if (cmd.getOptionValue(ARG_GTF) != null && cmd.getOptionValue(ARG_GTF).endsWith(".gtf")) {
            annFilePath = cmd.getOptionValue(ARG_GTF);

        }
        else if (cmd.getOptionValue(ARG_ANN) != null && (cmd.getOptionValue(ARG_ANN).endsWith(".gff3") || cmd.getOptionValue(ARG_ANN).endsWith(".gtf"))) {
            annFilePath = cmd.getOptionValue(ARG_ANN);
        }
        else {
            log.info(" *** Annotation Input Error: Please provide valid input for -ann, -gtf or -gff. Allowed file extensions are .gtf or .gff3 ***");
            log.info(" *** Note: -ann accepts both .gtf and .gff3 ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        if (fastaFilePath == null || !(fastaFilePath.endsWith(".fasta") || fastaFilePath.endsWith(".fa"))) {
            log.info(" *** Please provide valid input for -fasta. The input filename has to end with .fa or .fasta *** ");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        if (fastaGenomeFilePath != null && !(fastaGenomeFilePath.endsWith(".fasta") || fastaGenomeFilePath.endsWith(".fa"))) {
            log.info(" *** Please provide valid input for -genome. The input filename has to end with .fa or .fasta *** ");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        if (peptideInputFilePathsParam == null) {
            log.info(" *** Please provide valid input for -in. Allowed file extensions are .txt, .tsv or .pogo (e.g. filename.txt or filename1.txt,filename2.txt)  ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        String[] peptideInputFilePaths = Utils.tokenize(peptideInputFilePathsParam, ",", true);
        String[] validpeptideInputFileExts = {".txt", ".tsv", ".pogo", ".mztab", ".mzid", ".avro"};
        if (Stream.of(peptideInputFilePaths)
                .filter(filePath -> Stream.of(validpeptideInputFileExts).anyMatch(filePath::endsWith)).count() != peptideInputFilePaths.length) {
            log.info(" *** Please provide valid input for -in. Allowed file extensions are .mztab, .mzid, .txt, .tsv, .avro or" +
                    " .pogo (e.g. filename.txt or filename1.txt,filename2.txt) ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        String mergeParam = cmd.getOptionValue(ARG_MERGE);
        if (mergeParam != null && mergeParam.toLowerCase().startsWith("t")) {
            mergeFlag = true;
        }

        String formatParam = cmd.getOptionValue(ARG_FORMAT);
        if (formatParam == null || formatParam.toLowerCase().contains("all")) {
            gtfOutFlag = true;
            gctOutFlag = true;
            bedOutFlag = true;
            ptmbedOutFlag = true;
        } else {
            String[] formats = Utils.tokenize(formatParam.toLowerCase(), ",", true);
            ArrayList<String> formatsList = new ArrayList<>(Arrays.asList(formats));
            if (!formatsList.contains("gtf")) {
                gtfOutFlag = false;
            }
            if (!formatsList.contains("gct")) {
                gctOutFlag = false;
            }
            if (!formatsList.contains("bed")) {
                bedOutFlag = false;
            }
            if (!formatsList.contains("ptmbed")) {
                ptmbedOutFlag = false;
            }
        }

        String sourceParam = cmd.getOptionValue(ARG_SOURCE);
        if (sourceParam != null) {
            source = sourceParam;
        }

        String mmParam = cmd.getOptionValue(ARG_MM);
        if (mmParam != null) {
            int par = -1;
            try {
                par = Integer.parseInt(mmParam);
            } catch (Exception e) {
                System.err.println("ERROR: -mm param: invalid input received : " + mmParam);
            }
            if (par >= 0 && par <= 2) {
                GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES = par;
            } else {
                System.err.println("-mm: allowed mismatches need to be between 0 and 2. default (0) assumed");
            }
        }

        String mmModeParam = cmd.getOptionValue(ARG_MMMODE);
        if (mmModeParam != null) {
            if (mmModeParam.toLowerCase().startsWith("t")) {
                if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES > 1) {
                    GenomeMapper.PEPTIDE_MAPPER.ONE_IN_FIVE_MODE = true;
                } else {
                    System.err.println("-mmmode: cannot use mode with less than 2 mismatches. default (false) assumed.");
                }
            } else {
                System.err.println("-mmmode: invalid input. default (false) assumed.");
            }
        }

        String chrParam = cmd.getOptionValue(ARG_CHR);
        if (chrParam != null) {
            int par = -1;
            try {
                par = Integer.parseInt(chrParam);
            } catch (Exception e) {
                System.err.println("ERROR: -chr param: invalid input received : " + chrParam);
            }
            if (par == 1) {
                chrincluded = true;
            } else if (par != 0) {
                System.err.println("ERROR: The Chromosome prefix is not valid " + par + ". Export chr prefix Allowed (0, 1)  default (0) assumed");
            }
        }

        //files cannot be merged if there is only one input file
        if (mergeFlag && peptideInputFilePaths.length == 1) {
            System.err.println("cannot merge output files for one input file, default (-merge false) assumed");
            mergeFlag = false;
        }

        String pluralString = (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES == 1) ? " mismatch" : " mismatches";
        log.info("Start: allowing " + GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES + pluralString);


        try {

            if (fastaGenomeFilePath != null) {
                log.info("reading genome FASTA: " + fastaGenomeFilePath);
                GenomeFastaParser.readGenomeFASTA(fastaGenomeFilePath);
            }

            log.info("reading FASTA: " + fastaFilePath);
            CoordinateWrapper coordinate_wrapper = new CoordinateWrapper();
            coordinate_wrapper.read_fasta_file(fastaFilePath);

            log.info("Fasta done: " + coordinate_wrapper.size() + " proteins read.");
            log.info("building KmerTreeMap...");

            int kmerSize = (coordinate_wrapper.getTotalAACount() / coordinate_wrapper.getNumberOfProteins());
            kmerSize = (kmerSize / GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH) * coordinate_wrapper.getNumberOfProteins();

            IKmerMap kmer_map;
            if (inMemory)
                kmer_map = new KmerTreeMap();
            else
                kmer_map = new KmerSortedMap(kmerSize);

            coordinate_wrapper.add_all_proteins_to_kmer_map(kmer_map);

            if (!inMemory)
                ((KmerSortedMap) kmer_map).sortKmer();

            log.info("KmerTreeMap done: " + kmer_map.size() + " unique " + GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH + "-mers created.");
            log.info("reading GTF: " + annFilePath);

            MappedPeptides mapped_peptides = new MappedPeptides();
            // Check the file extension.  If GTF, run the GTFParser, else if GFF3, run the GFF3Parser.

            Assembly assem = null;
            if (annFilePath.endsWith(".gtf")) {
                assem = GTFParser.get_instance().read(annFilePath, coordinate_wrapper, mapped_peptides);
            } else if (annFilePath.endsWith(".gff3")) {
                assem = GFF3Parser.get_instance().read(annFilePath, coordinate_wrapper, mapped_peptides);
            }
            log.info("Annotation (GTF/GFF3) done!");

            // Creating a post-fix for file name specifying mode of mapping using mismatches
            String filename_mm_postfix = "";
            if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES > 0) {
                StringBuilder ss = new StringBuilder();
                ss.append("_").append(GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES).append("MM");

                if (GenomeMapper.PEPTIDE_MAPPER.ONE_IN_FIVE_MODE) {
                    ss.append("-(1in5)");
                }

                filename_mm_postfix = ss.toString();
            }

            for (String peptideInputFilePath : peptideInputFilePaths) {
                log.info("Computing genomic coordinates for: " + peptideInputFilePath);
                String final_peptide_path_results = FilenameUtils.removeExtension(peptideInputFilePath);

//                ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(curr_input_file_path, ".")));

                String path6 = final_peptide_path_results + "_unmapped.txt";

                if (fileFormat == INPUT_FILE_FORMAT.MZTAB)
                    new MzTabInputPeptideFileParser().read(peptideInputFilePath, coordinate_wrapper, mapped_peptides, path6, kmer_map);
                else if (fileFormat == INPUT_FILE_FORMAT.PEPTIDEATLAS)
                    new PeptideAtlasPeptideParser().read(peptideInputFilePath, coordinate_wrapper, mapped_peptides, path6, kmer_map);
                else if(fileFormat == INPUT_FILE_FORMAT.PRIDE){
                    new PeptideAvroFileParser(final_peptide_path_results + "_gene_coordinates.avro").read(peptideInputFilePath,  coordinate_wrapper, mapped_peptides, path6, kmer_map);
                }else
                    new TabInputPeptideFileParser().read(peptideInputFilePath, coordinate_wrapper, mapped_peptides, path6, kmer_map);

                log.info("Results done! (" + peptideInputFilePath + ")");
                log.info("writing output files");

                if (!mergeFlag) {
                    //the gtf overwrites the input gtf if they are in the same folder

                    String path4 = final_peptide_path_results + filename_mm_postfix + "_out.gtf";
                    String path5 = final_peptide_path_results + filename_mm_postfix + ".bed";
                    String path7 = final_peptide_path_results + filename_mm_postfix + ".gct";
                    String path8 = final_peptide_path_results + filename_mm_postfix + "_ptm.bed";
                    String path81 = final_peptide_path_results + filename_mm_postfix + "_no-ptm.bed";
                    String path9 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_out.gtf";
                    String path10 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + ".bed";
                    String path11 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + ".gct";
                    String path12 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_ptm.bed";
                    String path121 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_no-ptm.bed";

                    if (assem == Assembly.patchhaploscaff) {
                        path9 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_out.gtf";
                        path10 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff.bed";
                        path11 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff.gct";
                        path12 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_ptm.bed";
                        path121 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_no-ptm.bed";
                    }
                    if (gtfOutFlag) {
                        mapped_peptides.to_gtf(path4, source);
                        mapped_peptides.to_gtf(path9, source, assem);
                    }
                    if (bedOutFlag) {
                        mapped_peptides.to_bed(path5, Assembly.primary, chrincluded);
                        mapped_peptides.to_bed(path10, assem, chrincluded);
                    }
                    if (gctOutFlag) {
                        mapped_peptides.to_gct(path7, Assembly.primary, chrincluded);
                        mapped_peptides.to_gct(path11, assem);
                    }
                    if (ptmbedOutFlag) {
                        mapped_peptides.to_ptmbed(path8, path81);
                        mapped_peptides.to_ptmbed(path12, path121, assem);
                    }
                    mapped_peptides.remove_all_peptides();
                }
            }
            if (mergeFlag) {
//                ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(peptideInputFilePaths[0], ".")));

                String final_peptide_path_results = FilenameUtils.removeExtension(peptideInputFilePaths[0]);

                String path4 = final_peptide_path_results + filename_mm_postfix + "_merged.gtf";
                String path5 = final_peptide_path_results + filename_mm_postfix + "_merged.bed";
                String path7 = final_peptide_path_results + filename_mm_postfix + "_merged.gct";
                String path8 = final_peptide_path_results + filename_mm_postfix + "_merged_ptm.bed";
                String path81 = final_peptide_path_results + filename_mm_postfix + "_merged_no-ptm.bed";
                String path9 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.gtf";
                String path10 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.bed";
                String path11 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged.gct";
                String path12 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged_ptm.bed";
                String path121 = final_peptide_path_results + filename_mm_postfix + "_" + assem.toString() + "_merged_no-ptm.bed";

                if (assem == Assembly.patchhaploscaff) {

                    path9 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.gtf";
                    path10 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.bed";
                    path11 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged.gct";
                    path12 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged_ptm.bed";
                    path121 = final_peptide_path_results + filename_mm_postfix + "_patch_hapl_scaff_merged_no-ptm.bed";
                }

                if (gtfOutFlag) {
                    mapped_peptides.to_gtf(path4, source);
                    mapped_peptides.to_gtf(path9, source, assem);
                }
                if (bedOutFlag) {
                    mapped_peptides.to_bed(path5);
                    mapped_peptides.to_bed(path10, assem);
                }
                if (gctOutFlag) {
                    mapped_peptides.to_gct(path7);
                    mapped_peptides.to_gct(path11, assem);
                }
                if (ptmbedOutFlag) {
                    mapped_peptides.to_ptmbed(path8, path81);
                    mapped_peptides.to_ptmbed(path12, path121, assem);
                }
            }
            //if there is a problem with the reading of crucial files the program will end prematurely.
        } catch (Exception e) {
            log.info(e.getMessage());
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
        log.info("DONE..");
        long endTime = System.nanoTime();
        long totalTime = (long) ((endTime - startTime) / 1000000000.0);
        log.debug("Running time -- " + totalTime + " Min");
    }
}
