package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.*;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GTFParser extends GenomeAnnotationParser {

    private static Logger log = LoggerFactory.getLogger(GTFParser.class);

    public static GTFParser instance;

    private GTFParser() {
    }

    //pattern for gene ID
    private static final Pattern GTFGENEPATTERN = Pattern.compile("gene_id \"([^\"\\.]*)[^\"]*\"");
    //pattern for transcript ID
    private static final Pattern GTFTRANSCRIPTPATTERN = Pattern.compile("transcript_id \"([^\"\\.]*)[^\"]*\"");
    //pattern for exon ID
    private static final Pattern GTFEXONPATTERN = Pattern.compile("exon_id \"([^\"\\.]*)[^\"]*\"");

    //singleton get_instance method.
    public static GTFParser get_instance() {
        if (instance == null) {
            instance = new GTFParser();
        }
        return instance;
    }

    //reads a gtf file and parses it into CoordinateWrapper and MappedPeptides.
    public final Assembly read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping) throws Exception {
        if (!open(file)) {
            throw new IllegalStateException("Problem in reading GTF file");
        }

        String exonID = "";
        ProteinEntry proteinEntry = null;
        ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinatesMap = new ArrayList<>();
//        Coordinates protein_coordinates = new Coordinates();
        Coordinates prevProteinCoordinates = new Coordinates();
        Assembly assem = Assembly.none;
        ArrayList<String> tokens;
        while ((line = reader.readLine()) != null) {
            if ((line.startsWith("#"))) {
                continue;
            }
            tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t")));

            if (is_next_gene(tokens)) {
                Assembly assemtemp = mapping.add_gene_from_annotation(line);
                if (assem == Assembly.none) {
                    if (assemtemp == Assembly.patchhaploscaff) {
                        assem = assemtemp;
                    }
                }
            }
            String transcriptId = extract_transcript_id(line);
            if (is_next_transcript(tokens)) {
                exonID = "";
                mapping.add_transcript_id_to_gene(line);
                if (proteinEntry != null) {
                    proteinEntry.set_coordinate_map(coordinatesMap);
                }
                proteinEntry = coordwrapper.lookup_entry(transcriptId);
                if (proteinEntry == null) {
                    log.info("ERROR: No entry for transcript ID: " + transcriptId);
                    continue;
                }
//                protein_coordinates = new Coordinates();
                prevProteinCoordinates = new Coordinates();
                prevProteinCoordinates.setCterm(Offset.off3);
                prevProteinCoordinates.setNterm(Offset.off3);
                prevProteinCoordinates.setStart(0);
                prevProteinCoordinates.setEnd(0);
                coordinatesMap = new ArrayList<>();
            } else if (is_exon(tokens)) {
                exonID = extract_exon_id(line);
            } else if (is_cds(tokens)) {
                GenomeCoordinates genCoord = Utils.extract_coordinates_from_gtf_line(tokens);
                genCoord.setTranscriptid(transcriptId);
                String tmp_exonID = extract_exon_id(line);
                if(tmp_exonID.equals("")) {
                    tmp_exonID = exonID;
                }
                genCoord.setExonid(tmp_exonID);
                Coordinates proteinCoordinates = new Coordinates();
                // get nterm from prev exon
                if (genCoord.getFrame() != Frame.unknown) {
                    proteinCoordinates.setNterm(Offset.forValue(genCoord.getFrame().getValue()));
                } else {
                    if (prevProteinCoordinates.getCterm() != Offset.off3) {
                        proteinCoordinates.setNterm(Offset.forValue(3 - prevProteinCoordinates.getCterm().getValue()));
                    } else {
                        proteinCoordinates.setNterm(Offset.off3);
                    }
                }

                int length = 0;

                if (is_first_strand(tokens)) {
                    length = genCoord.getEnd() - genCoord.getStart() + 1;
                } else if (!is_first_strand(tokens)) {
                    length = genCoord.getStart() - genCoord.getEnd() + 1;
                }

                // calc cterm
                if (length % 3 == 0) {
                    if (proteinCoordinates.getNterm() != Offset.off3) {
                        proteinCoordinates.setCterm(Offset.forValue(3 - proteinCoordinates.getNterm().getValue()));
                    } else {
                        proteinCoordinates.setCterm(Offset.off3);
                    }
                } else if (length % 3 == 2) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off2);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off1);
                    }
                } else if (length % 3 == 1) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off1);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off2);
                    }
                }

                // calc protein coordinates
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    proteinCoordinates.setStart(prevProteinCoordinates.getEnd());
                } else {
                    if (prevProteinCoordinates.getEnd() == 0 && coordinatesMap.isEmpty()) {
                        proteinCoordinates.setStart(0);
                    } else {
                        proteinCoordinates.setStart(prevProteinCoordinates.getEnd() + 1);
                    }
                }

                int offsets = 0;
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    offsets = offsets + proteinCoordinates.getNterm().getValue();
                }

                if (is_first_strand(tokens)) {
                    length = genCoord.getEnd() - genCoord.getStart() + 1 - offsets;
                } else if (!is_first_strand(tokens)) {
                    length = genCoord.getStart() - genCoord.getEnd() + 1 - offsets;
                }

                int peplength = length / 3;

                int pepend = proteinCoordinates.getStart() + peplength - 1;
                if (proteinCoordinates.getCterm() != Offset.off3) {
                    pepend = pepend + 1;
                }
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    pepend = pepend + 1;
                }

                proteinCoordinates.setEnd(pepend);

                prevProteinCoordinates = proteinCoordinates;

                coordinatesMap.add(new Tuple<>(proteinCoordinates, genCoord));
            }
        }
        if (proteinEntry != null) {
            proteinEntry.set_coordinate_map(coordinatesMap);
        }
        close();
        return assem;
    }

    //looks for the text specified GENEPATTERN and returns the ID.
    public static String extract_gene_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GTFGENEPATTERN);
    }

    //looks for the text specified in TRNASCRIPTPATTERN and returns the ID.
    public static String extract_transcript_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GTFTRANSCRIPTPATTERN);
    }

    //looks for the text specified in EXONPATTERN and returns the ID.
    public static String extract_exon_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GTFEXONPATTERN);
    }

    public static String extract_id(String gtfGeneLine, Pattern pattern) {
        String value = "";
        Matcher matcher = pattern.matcher(gtfGeneLine);
        if (matcher.find()) {
            value = matcher.group(1);
        }
        return value;
    }

    // after tokenizing the gene line these functions can be used to extract information.
    // extracts the gene type (protein_coding,...)
    public static String extract_type(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = GeneEntry.extract_by_tag("gene_type", tokens.get(8));
            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }

    //extracts the gene status (KNOWN,...)
    public static String extract_status(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = GeneEntry.extract_by_tag("gene_status", tokens.get(8));
            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }

    //extracts the gene symbol
    public static String extract_gene_name(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = GeneEntry.extract_by_tag("gene_name", tokens.get(8));

            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }
}