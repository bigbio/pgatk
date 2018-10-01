package org.bigbio.pgatk.pepgenome;

import org.bigbio.pgatk.pepgenome.common.*;
import org.apache.commons.cli.*;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import org.bigbio.pgatk.pepgenome.kmer.inmemory.KmerMap;
import org.ehcache.sizeof.SizeOf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class PepGenomeTool {
    //exit codes --------------------------------
    private static final int GENOME_MAPPER_EXIT_HELP = 1;
    private static final int GENOME_MAPPER_EXIT_TOO_FEW_ARGS = 2;
    private static final int GENOME_MAPPER_EXIT_INVALID_ARG = 3;

    //-----------------input args--------------------------
    private static final String ARG_FASTA = "fasta";
    private static final String ARG_GTF = "gtf";
    private static final String ARG_IN = "in";
    private static final String ARG_MERGE = "merge";
    private static final String ARG_FORMAT = "format";
    private static final String ARG_SOURCE = "source";
    private static final String ARG_MM = "mm";
    private static final String ARG_MMMODE = "mmmode";
    private static final String ARG_SPECIES = "species";
    private static final String ARG_CHR = "chr";
    private static final String ARG_HELP = "h";

    //DEFAULT values
    private static boolean mergeFlag = false;
    private static boolean gtfOutFlag = true;
    private static boolean gctOutFlag = true;
    private static boolean bedOutFlag = true;
    private static boolean ptmbedOutFlag = true;
    private static String source = "PoGo";
    private static boolean chrincluded = false;


    public static void main(String[] args) {
        Options options = new Options();
        options.addOption(Option.builder(ARG_FASTA).hasArg(true).desc("Filepath for file containing protein sequences in FASTA format").build())
                .addOption(Option.builder(ARG_GTF).hasArg(true).desc("Filepath for file containing genome annotation in GTF format").build())
                .addOption(Option.builder(ARG_IN).hasArg(true).desc("Comma(,) separated filepaths for files containing peptide identifications (Contents of the file should be in tab seperated format. i.e., File format: four columns: SampleName\t\tPeptideSequence\t\tPSMs\tQuant)").build())
                .addOption(Option.builder(ARG_MERGE).hasArg(true).desc("Set 'true' to merge mappings from all files from input (default 'false')").build())
                .addOption(Option.builder(ARG_FORMAT).hasArg(true).desc("Select the output formats from gtf, gct, bed, ptmbed, all or combinations thereof separated by ',' (default all)").build())
                .addOption(Option.builder(ARG_SOURCE).hasArg(true).desc("Please give a source name which will be used in the second column in the output gtf file (default: PoGo)").build())
                .addOption(Option.builder(ARG_MM).hasArg(true).desc("Allowed mismatches (0, 1 or 2; default: 0)").build())
                .addOption(Option.builder(ARG_MMMODE).hasArg(true).desc("Mismatch mode (true or false): if true mismatching with two mismaches will only allow 1 mismatch every kmersize (default: 5) positions. (default: false)").build())
                .addOption(Option.builder(ARG_SPECIES).hasArg(true).desc("Please give species using common or scientific name (default human). For a full list of supported species please go to https://github.com/bigbio/pgatk/tree/master/PepGenome").build())
                .addOption(Option.builder(ARG_CHR).hasArg(true).desc("Export chr prefix Allowed 0, 1  (default: 0)").build())
                .addOption(Option.builder(ARG_HELP).hasArg(false).desc("Print this help & exit").build());

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
            if (cmd.hasOption(ARG_HELP)) {
                Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_HELP);
            }
        } catch (ParseException e) {
            System.out.println(" *** Error in parsing the input arguments. Please check the arguments ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_TOO_FEW_ARGS);
        }

        if (!cmd.hasOption(ARG_FASTA) || !cmd.hasOption(ARG_GTF) || !cmd.hasOption(ARG_IN)) {
            System.out.println("*** Missing mandatory parameters: -fasta, -gtf and -in ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_TOO_FEW_ARGS);
        }

        String fastaFilePath = cmd.getOptionValue(ARG_FASTA);
        String gtfFilePath = cmd.getOptionValue(ARG_GTF);
        String peptideInputFilePathsParam = cmd.getOptionValue(ARG_IN);

        if (fastaFilePath == null || !(fastaFilePath.endsWith(".fasta") || fastaFilePath.endsWith(".fa"))) {
            System.out.println(" *** Please provide valid input for -fasta. The input filename has to end with .fa or .fasta *** ");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        if (gtfFilePath == null || !gtfFilePath.endsWith(".gtf")) {
            System.out.println(" *** Please provide valid input for -gtf. The input filename has to end with .gtf  ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        if (peptideInputFilePathsParam == null) {
            System.out.println(" *** Please provide valid input for -in. Allowed file extensions are .txt, .tsv or .pogo (e.g. filename.txt or filename1.txt,filename2.txt)  ***");
            Utils.printHelpAndExitProgram(options, true, GENOME_MAPPER_EXIT_INVALID_ARG);
        }

        String[] peptideInputFilePaths = Utils.tokenize(peptideInputFilePathsParam, ",", true);
        String[] validpeptideInputFileExts = {".txt", ".tsv", ".pogo"};
        if (Stream.of(peptideInputFilePaths)
                .filter(filePath -> Stream.of(validpeptideInputFileExts).anyMatch(filePath::endsWith))
                .collect(Collectors.toList())
                .size() != peptideInputFilePaths.length) {
            System.out.println(" *** Please provide valid input for -in. Allowed file extensions are .txt, .tsv or .pogo (e.g. filename.txt or filename1.txt,filename2.txt) ***");
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

        String speciesParam = cmd.getOptionValue(ARG_SPECIES);
        if (speciesParam != null) {
            speciesParam = speciesParam.toLowerCase();
            if (GenomeMapper.TAX.containsKey(speciesParam)) {
                GenomeMapper.ID.GENE_ID = GenomeMapper.TAX.get(speciesParam).getGeneId();
                GenomeMapper.ID.TRANSCRIPT_ID = GenomeMapper.TAX.get(speciesParam).getTranscriptId();
                GenomeMapper.ID.EXON_ID = GenomeMapper.TAX.get(speciesParam).getExonId();
                GenomeMapper.ID.LENGTH = GenomeMapper.TAX.get(speciesParam).getLength();
            } else {
                System.err.println("ERROR: Species/Taxonomy: " + speciesParam + "is not supported. For a full list of supported species please go to https://github.com/bigbio/pgatk/tree/master/PepGenome");
                System.exit(GENOME_MAPPER_EXIT_HELP);
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
        System.out.println("Start: allowing " + GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES + pluralString);
        System.out.println("reading FASTA: " + fastaFilePath);

        try {
            CoordinateWrapper coordinate_wrapper = new CoordinateWrapper();
            coordinate_wrapper.read_fasta_file(fastaFilePath);

            System.out.println("Fasta done: " + coordinate_wrapper.size() + " proteins read.");
            System.out.println("building KmerMap...");
            IKmerMap kmer_map = new KmerMap();
            coordinate_wrapper.add_all_proteins_to_kmer_map(kmer_map);

            SizeOf sizeOf = SizeOf.newInstance();
            System.out.println((int) (sizeOf.deepSizeOf(kmer_map)/1000000) + " MB");

            System.out.println("KmerMap done: " + kmer_map.size() + " unique " + GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH + "-mers created.");
            System.out.println("reading GTF: " + gtfFilePath);

            MappedPeptides mapped_peptides = new MappedPeptides();
            Assembly assem = GTFParser.get_instance().read(gtfFilePath, coordinate_wrapper, mapped_peptides);
            System.out.println("GTF done!\nComputing genomic coordinates for: ");

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
                System.out.println(peptideInputFilePath);

                String final_peptide_path_results = peptideInputFilePath;
                if (peptideInputFilePath.endsWith(".txt")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".txt");
                } else if (peptideInputFilePath.endsWith(".tsv")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".tsv");
                }

//                ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(curr_input_file_path, ".")));

                String path6 = final_peptide_path_results + "_unmapped.txt";

                ResultParser.read(peptideInputFilePath, coordinate_wrapper, mapped_peptides, path6, kmer_map);
                System.out.println("Results done! (" + peptideInputFilePath + ")");
                System.out.println("writing output files");

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

                String final_peptide_path_results = peptideInputFilePaths[0];
                if (final_peptide_path_results.endsWith(".txt")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".txt");
                } else if (final_peptide_path_results.endsWith(".tsv")) {
                    final_peptide_path_results = Utils.removeExtension(final_peptide_path_results, ".tsv");
                }

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
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
        System.out.println("DONE..");
    }
}