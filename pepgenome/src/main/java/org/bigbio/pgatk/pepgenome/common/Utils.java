package org.bigbio.pgatk.pepgenome.common;


import io.github.bigbio.pgatk.io.pride.ExonInfo;
import io.github.bigbio.pgatk.io.pride.GeneCoordinates;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;

import java.util.*;

public class Utils {

    public static void printHelpAndExitProgram(final Options options, boolean shouldExit, int exitCode) {
        new HelpFormatter().printHelp("Arguments: -fasta TRANSL -gtf ANNO -in *.tsv[,*.tsv] [-format OUTF] [-merge TRUE/FALSE] [-source SRC] [-mm NUM] [-mmmode TRUE/FALSE] [-species SPECIES] [-chr 0/1]", options);
        if (shouldExit) {
            System.exit(exitCode);
        }
    }

    public static String[] tokenize(String string2split) {
        return tokenize(string2split, " ", false);
    }

    public static String[] tokenize(String string2split, String regex) {
        return tokenize(string2split, regex, false);
    }

    public static String[] tokenize(String string2split, boolean trimEmpty) {
        return tokenize(string2split, " ", trimEmpty);
    }

    public static String[] tokenize(String string2split, String regex, boolean trimEmpty) {
        String[] splits = string2split.split(regex, -1);
        if (!trimEmpty) {
            return splits;
        }
        return Arrays.stream(splits)
                .filter(value -> value != null && value.length() > 0)
                .toArray(String[]::new);
    }

    //removes all ptms from a sequence. these must be delimited by '(' and ')', respectively
    //leaves the original string unchanged.
    public static String remove_ptms(String sequence) {
        String tmp = sequence;
        int start = tmp.indexOf("(");
        while (start != -1) {
            int end = tmp.indexOf(")");
            if (end == -1) {
                break;
            }
            String substr = tmp.substring(start, end + 1);
            tmp = StringUtils.remove(tmp, substr);
            start = tmp.indexOf("(");
        }

        return tmp.toUpperCase();
    }




    /**
     * Converts a sequence into an isosequence. this means replacing all 'I' and 'L' chars with 'J'
     * leaves the original string unchanged.
     *
     * @param sequence Original sequence
     * @return new sequence with all I and J replaced.
     */
    public static String make_iso_sequence(String sequence) {
        return sequence.replaceAll("I", "J").replaceAll("L", "J");
    }

    //this function creates the an entry for a CoordinateMapType.
    public static Tuple<Coordinates, GenomeCoordinates> get_coordinates(Coordinates proteinCoords, GenomeCoordinates genomeCoords, Coordinates peptideCoords) {
        int pep_start;
        int pep_end;

        if (peptideCoords.start >= proteinCoords.start && peptideCoords.start <= proteinCoords.end) {
            //if the start of the peptide lies within the protein use this as the protein start
            pep_start = peptideCoords.start;
        } else {
            //otherwise the peptide has to start where the protein starts.
            pep_start = proteinCoords.start;
        }

        if (peptideCoords.end >= proteinCoords.start && peptideCoords.end <= proteinCoords.end) {
            //if the end of the also lies within the protein use this as peptide end
            pep_end = peptideCoords.end;
        } else {
            //otherwise the peptide has to end where the protein ends
            pep_end = proteinCoords.end;
        }
        Coordinates partial_peptide_coords = new Coordinates();
        partial_peptide_coords.start = pep_start;
        partial_peptide_coords.end = pep_end;
        int start_genomic_coord;
        int end_genomic_coord;

        if (genomeCoords.end >= genomeCoords.start) {
            //coding from forward strand
            if (pep_start == proteinCoords.start) {
                //peptide starting at protein start position
                partial_peptide_coords.Nterm = proteinCoords.Nterm;
                start_genomic_coord = genomeCoords.start;
            } else {
                //peptide starting within the protein
                partial_peptide_coords.Nterm = Offset.off3;
                start_genomic_coord = ((pep_start - 1 - proteinCoords.start) * 3) + genomeCoords.start + proteinCoords.Nterm.getValue();
            }
            if (pep_end == proteinCoords.end) {
                //peptide ending at protein end position
                partial_peptide_coords.Cterm = proteinCoords.Cterm;
                end_genomic_coord = genomeCoords.end;
            } else {
                //peptide ending within protein
                partial_peptide_coords.Cterm = Offset.off3;
                end_genomic_coord = ((pep_end - proteinCoords.start) * 3) + genomeCoords.start + (proteinCoords.Nterm.getValue() - 1);
            }
        } else {
            // coding from reverse strand!
            if (pep_start == proteinCoords.start) {
                //peptide starting at protein start position
                partial_peptide_coords.Nterm = proteinCoords.Nterm;
                start_genomic_coord = genomeCoords.start;
            } else {
                //peptide starting within the protein
                partial_peptide_coords.Nterm = Offset.off3;
                start_genomic_coord = genomeCoords.start - (((pep_start - proteinCoords.start - 1) * 3) + proteinCoords.Nterm.getValue());
            }
            if (pep_end == proteinCoords.end) {
                //peptide ending at protein end position
                partial_peptide_coords.Cterm = proteinCoords.Cterm;
                end_genomic_coord = genomeCoords.end;
            } else {
                //peptide ending within protein
                partial_peptide_coords.Cterm = Offset.off3;
                end_genomic_coord = genomeCoords.start - (((pep_end - proteinCoords.start) * 3) + proteinCoords.Nterm.getValue() - 1);
            }
        }

        GenomeCoordinates genome_coordinates = new GenomeCoordinates();
        genome_coordinates.setChr(genomeCoords.getChr());
        genome_coordinates.setChrscaf(genomeCoords.getChrscaf());
        genome_coordinates.setStrand(genomeCoords.getStrand());
        genome_coordinates.setFrame(genomeCoords.getFrame());
        genome_coordinates.start = start_genomic_coord;
        genome_coordinates.end = end_genomic_coord;
        genome_coordinates.setTranscriptid(genomeCoords.getTranscriptid());
        genome_coordinates.setExonid(genomeCoords.getExonid());

        //flip start and end if coding from reverse strand
        if (genome_coordinates.getStrand() == Strand.rev) {
            int start = genome_coordinates.end;
            int end = genome_coordinates.start;
            genome_coordinates.start = start;
            genome_coordinates.end = end;
        }

        return new Tuple<>(partial_peptide_coords, genome_coordinates);
    }

    //a more sophisticated operator== for common.GenomeCoordinates.
    public static boolean same_coordinates(GenomeCoordinates lhs, GenomeCoordinates rhs) {
        return ((!lhs.getChr().isScaffold() && lhs.getChr() == rhs.getChr())
                || (lhs.getChr().isScaffold() && lhs.getChrscaf().equals(rhs.getChrscaf())))
                && (lhs.start == rhs.start) && (lhs.end == rhs.end)
                && (lhs.getFrame() == rhs.getFrame())
                && (lhs.getStrand() == rhs.getStrand());
    }

    //this compare_function returns true if lhs.start < rhs.start and false otherwise.
    //it should ONLY be used in sort funtions, where only the rough order matters
    //for comparisions if two GenomicCoordinates are equal, use compare_genome_coordinate_sets_ascending or compare_coordinates_ascending_whole
    public static boolean compare_coordinates_ascending(GenomeCoordinates lhs, GenomeCoordinates rhs) {
        return lhs.start < rhs.start;
    }

    //returns true if lhs.start < rhs.start
    //otherwise returns true if lhs.end < rhs.end
    //otherwise returns false.
    public static boolean compare_coordinates_ascending_whole(GenomeCoordinates lhs, GenomeCoordinates rhs) {
        return lhs.start < rhs.start || (lhs.start == rhs.start && lhs.end < rhs.end);
    }

    public static String coordinates_to_string(GenomeCoordinates coords) {
        return coordinates_to_string(coords, true);
    }

    private static StringBuilder coord2StrCommon(GenomeCoordinates coords, boolean chrincluded, StringBuilder ss) {
        if (coords.getChr().isScaffold() && !chrincluded) {
            ss.append(coords.getChrscaf());
        } else if (coords.getChr().isScaffold() && chrincluded) {
            ss.append("scaffold").append(coords.getChrscaf());
        } else if (!chrincluded) {
            ss.append(EnumStringMapper.enumToString(coords.getChr()));
        } else if (chrincluded) {
            ss.append(EnumStringMapper.enumToChrString(coords.getChr()));
        }
        return ss;
    }

    //to_string functions for Genomic coordinates.
    //as simple string
    public static String coordinates_to_string(GenomeCoordinates coords, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        ss = coord2StrCommon(coords, chrincluded, ss);
        ss.append(":").append(coords.start).append("-").append(coords.end).append(" ").append(EnumStringMapper.enumToString(coords.getStrand(), false));
        return ss.toString();
    }

    public static String coordinates_to_short_string(GenomeCoordinates coords, int offset) {
        return coordinates_to_short_string(coords, offset, true);
    }

    public static String coordinates_to_short_string(GenomeCoordinates coords, boolean chrincluded) {
        return coordinates_to_short_string(coords, 0, chrincluded);
    }

    public static String coordinates_to_short_string(GenomeCoordinates coords) {
        return coordinates_to_short_string(coords, 0, true);
    }

    //as a shortend version.
    public static String coordinates_to_short_string(GenomeCoordinates coords, int offset, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        ss = coord2StrCommon(coords, chrincluded, ss);
        ss.append(":").append(coords.start - offset).append("-").append(coords.end);
        return ss.toString();
    }

    public static String coordinates_to_gtf_string(GenomeCoordinates coords, String type, boolean frameinclude, String source) {
        return coordinates_to_gtf_string(coords, type, frameinclude, source, true);
    }

    //as a line in a gtffile
    public static String coordinates_to_gtf_string(GenomeCoordinates coords, String type, boolean frameinclude, String source, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        ss = coord2StrCommon(coords, chrincluded, ss);
        ss.append("\t").append(source).append("\t").append(type).append("\t").append(coords.start).append("\t").append(coords.end).append("\t.\t").append(EnumStringMapper.enumToString(coords.getStrand(), false));
        if (frameinclude) {
            ss.append("\t").append(coords.getFrame()).append("\t");
        } else {
            ss.append("\t.\t");
        }
        return ss.toString();
    }

    public static String coordinates_to_bed_string(GenomeCoordinates coords, String name, int score) {
        return coordinates_to_bed_string(coords, name, score, true);
    }

    public static String coordinates_to_bed_string(GenomeCoordinates coords, String name) {
        return coordinates_to_bed_string(coords, name, 1000, true);
    }

    public static List<GeneCoordinates> coordinates_info(PeptideEntry peptide, String pSequence, boolean noptm, boolean chrincluded){

        Set<GeneCoordinates> summaryCoordinates = new HashSet<>();

        for (PeptideCoordinates coord : peptide.getPepCoordinates()) {
            //os.write(Utils.coordinates_to_bed_string(coord.get_transcript_coordinates(), pSequence).getBytes());
            //std::cout << coordinates_to_bed_string((*it)->get_transcript_coordinates(), pSequence) << std::endl;

            ArrayList<GenomeCoordinates> exon_coordinates = coord.get_exon_coordinates();
            int exon_count = 0;
            List<ExonInfo> exonInfoList = new ArrayList<>();
            for (GenomeCoordinates exon_coordinate : exon_coordinates) {
                exon_count += 1;

                int exon_start = exon_coordinate.getStart() - coord.get_transcript_coordinates().getStart();
                int exon_length = exon_coordinate.getEnd() - exon_coordinate.getStart() + 1;
                exonInfoList.add(ExonInfo
                        .builder()
                        .exonStarts(exon_start)
                        .exonLengths(exon_length)
                        .exonCount(exon_count)
                        .exonAccession(exon_coordinate.getExonid())
                        .build());
            }

            GeneCoordinates coords = GeneCoordinates.builder()
                    .start(coord.get_transcript_coordinates().start -1)
                    .end(coord.get_transcript_coordinates().end)
                    .chromosome(coord2StrCommon(coord.get_transcript_coordinates(), chrincluded, new StringBuilder("")).toString())
                    .transcriptAccession(coord.get_transcript_coordinates().getTranscriptid())
                    .geneUnique(peptide.isGeneUnique())
                    .transcriptUnique(peptide.isTranscriptUnique())
                    .exonInfoList(exonInfoList)
                    .geneAccession(peptide.getAssociatedGene().get_id())
                    .geneName(peptide.getAssociatedGene().get_name())
                    .geneType(peptide.getAssociatedGene().get_type())
                    .strand(coord.get_transcript_coordinates().getStrand().toString())
                    .build();
            summaryCoordinates.add(coords);
        }
        return new ArrayList<>(summaryCoordinates);
    }

    public static String coordinates_to_bed_string(GenomeCoordinates coords, String name, boolean chrincluded) {
        return coordinates_to_bed_string(coords, name, 1000, chrincluded);
    }

    //as a line in a bed file
    public static String coordinates_to_bed_string(GenomeCoordinates coords, String name, int score, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        ss = coord2StrCommon(coords, chrincluded, ss);
        ss.append("\t").append((coords.start - 1)).append("\t")
                .append(coords.end).append("\t").append(name).append("\t")
                .append(score).append("\t").append(EnumStringMapper.enumToString(coords.getStrand(), false)).append("\t").append((coords.start - 1)).append("\t")
                .append(coords.start - 1).append("\t");
        return ss.toString();
    }

    public static String coordinates_to_short_bed_string(GenomeCoordinates coords, String name, int score) {
        return coordinates_to_short_bed_string(coords, name, score, true);
    }

    public static String coordinates_to_short_bed_string(GenomeCoordinates coords, String name) {
        return coordinates_to_short_bed_string(coords, name, 1000, true);
    }

    public static String coordinates_to_short_bed_string(GenomeCoordinates coords, String name, boolean chrincluded) {
        return coordinates_to_short_bed_string(coords, name, 1000, chrincluded);
    }

    //as a short bed string
    public static String coordinates_to_short_bed_string(GenomeCoordinates coords, String name, int score, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        ss = coord2StrCommon(coords, chrincluded, ss);
        ss.append("\t").append((coords.start - 1)).append("\t").append(coords.end).append("\t").append(name).append("\t").append(score).append("\t").append(EnumStringMapper.enumToString(coords.getStrand(), false)).append("\t");
        return ss.toString();
    }

    public static String coordinates_to_gct_string(ArrayList<GenomeCoordinates> coords) {
        return coordinates_to_gct_string(coords, true);
    }

    //as a line in a gct file
    public static String coordinates_to_gct_string(ArrayList<GenomeCoordinates> coords, boolean chrincluded) {
        StringBuilder ss = new StringBuilder();
        for (int i = 0; i < coords.size(); ++i) {
            if (i > 0) {
                ss.append(",");
            }
            ss.append(coordinates_to_short_string(coords.get(i), 1, chrincluded));
        }
        return ss.toString();
    }

    //compares two genomic coordinates in a way that they can be sorted ascendingly.
    //this is important for map/set <common.GenomeCoordinates> because otherwise they would overwrite
    //other common.Coordinates even if they arent the same.
    public static boolean compare_genome_coordinate_sets_ascending(ArrayList<GenomeCoordinates> lhs, ArrayList<GenomeCoordinates> rhs) {
        boolean same = false;
        int lhs_it = 0;
        int rhs_it = 0;
        if (lhs.size() > 1 && (lhs.size() == rhs.size())) {
            //if both sets have the same amount of entries
            while (lhs_it < lhs.size() && rhs_it < rhs.size() && (!same)) {
                //compare the entries.
                same = compare_coordinates_ascending_whole(lhs.get(lhs_it), rhs.get(rhs_it));
                ++lhs_it;
                ++rhs_it;
            }
        } else if (lhs.size() == 1 && rhs.size() == 1) {
            //if both sets only contain one entry
            same = compare_coordinates_ascending_whole(lhs.get(lhs_it), rhs.get(rhs_it));
        } else {
            //if one contains less entries than another
            int it_end = (lhs.size() < rhs.size()) ? lhs.size() : rhs.size();
            for (int i = 0; (i < it_end) && (!same); ++i) {
                same = compare_coordinates_ascending_whole(lhs.get(i), rhs.get(i));
            }
        }
        return same;
    }

    //given a tokenized string this function will generate the resulting genomic coordinates.
    public static GenomeCoordinates extract_coordinates_from_gtf_line(List<String> tokens) {
        GenomeCoordinates coord = new GenomeCoordinates();
        coord.setChr(EnumStringMapper.string_to_chromosome(tokens.get(0)));
        if (coord.getChr().isScaffold()) {
            coord.setChrscaf(tokens.get(0));
        } else {
            coord.setChrscaf("");
        }
        coord.setStrand(EnumStringMapper.string_to_strand(tokens.get(6)));
        coord.setFrame(EnumStringMapper.string_to_frame(tokens.get(7)));
        if (coord.getStrand() == Strand.fwd) {
            coord.start = Integer.parseInt(tokens.get(3));
            coord.end = Integer.parseInt(tokens.get(4));
        } else if (coord.getStrand() == Strand.rev) {
            coord.end = Integer.parseInt(tokens.get(3));
            coord.start = Integer.parseInt(tokens.get(4));
        }
        return coord;
    }

    /**
     * This function remove the extension of the file if is found, if the extension is not found, it will keep
     * the original value.
     *
     * @param fileName  Name of the file.
     * @param extension Extension to be removed
     * @return new File Name
     */
    public static String removeExtension(String fileName, String extension) {
        if (fileName.endsWith(extension)) {
            int indx = fileName.lastIndexOf(extension);
            return fileName.substring(0, indx);
        }
        return fileName;
    }

    /**
     * Equivalent of "string.c_str() + index";
     */
    public static String getCppStyleSubStringByShift(String str, int from) {
        try {
            return str.substring(from);
        } catch (Exception ex) {
            return str;
        }
    }

    /**
     * We use the folowing function to retrieve the specific sub-string in a protein sequence.
     * It is equivalent to use {@link String} str.substr(from,len)
     * @param str String
     * @param from starting point in the String
     * @param len number of chars to retrieve.
     * @return SubString
     */
    public static String getCppStyleSubString(String str, int from, int len) {
        int end = ((from + len) > str.length()) ? str.length() : ((from + len));
        return str.substring(from, end);
    }

}