package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.*;
import javafx.util.Pair;
import org.bigbio.pgatk.pepgenome.common.*;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GTFParser {

    private static Logger log = LoggerFactory.getLogger(GTFParser.class);

    //inputstream
    private BufferedReader m_reader;
    private FileInputStream m_ifs;
    //current line
    private String m_line;
    //meyers singleton instance.
    private static GTFParser m_instance;

    private GTFParser() {
    }

    //singleton get_instance method.
    public static GTFParser get_instance() {
        if (m_instance == null) {
            m_instance = new GTFParser();
        }
        return m_instance;
    }

    //opens filestream, returns true if sucessful
    private boolean open(String file) throws Exception {
        if (m_reader == null) {
            m_line = "";
            m_ifs = new FileInputStream(file);
            m_reader = new BufferedReader(new InputStreamReader(m_ifs));
        }
        boolean status = true;
        try{
            status = m_reader.ready();
        }catch (IOException ex){
            log.debug("The gtf file stream is closed -- " + file);
            m_ifs = new FileInputStream(file);
            m_reader = new BufferedReader(new InputStreamReader(m_ifs));
        }
        return status;
    }

    //closes the filestream
    private void close() throws Exception {
        m_line = "";
        m_reader.close();
        m_ifs.close();
    }

    //returns true if in the GTF at position 6 there is a + (plus strand)
    private static boolean is_first_strand(List<String> tokens) {
        return tokens.get(6).equals("+");
    }

    //returns true if position 2 in the GTF says "CDS"
    private static boolean is_cds(List<String> tokens) {
        return tokens.get(2).equals("CDS");
    }

    //returns true if position 2 in the GTF says "transcript"
    private static boolean is_next_transcript(List<String> tokens) {
        return tokens.get(2).equals("transcript");
    }

    //returns true if position 2 in the GTF says "gene"
    private static boolean is_next_gene(List<String> tokens) {
        return tokens.get(2).equals("gene");
    }

    //reads a gtf file and parses it into CoordinateWrapper and MappedPeptides.
    public final Assembly read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping) throws Exception {
        if (!open(file)) {
            throw new IllegalStateException("Problem in reading GTF file");
        }

        ProteinEntry p_protein_entry = null;
        ArrayList<Pair<Coordinates, GenomeCoordinates>> coordinatesMap = new ArrayList<>();
//        Coordinates protein_coordinates = new Coordinates();
        Coordinates prevProteinCoordinates = new Coordinates();
        Assembly assem = Assembly.none;
        ArrayList<String> tokens;
        while ((m_line = m_reader.readLine()) != null) {
            if ((m_line.startsWith("#"))) {
                continue;
            }
            tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(m_line, "\t")));

            if (is_next_gene(tokens)) {
                Assembly assemtemp = mapping.add_gene_from_gtf(m_line);
                if (assem == Assembly.none) {
                    if (assemtemp == Assembly.patchhaploscaff) {
                        assem = assemtemp;
                    }
                }
            }
            String transcriptId = GeneEntry.extract_transcript_id(m_line);
            if (is_next_transcript(tokens)) {
                mapping.add_transcript_id_to_gene(m_line);
                if (p_protein_entry != null) {
                    p_protein_entry.set_coordinate_map(coordinatesMap);
                }
                p_protein_entry = coordwrapper.lookup_entry(transcriptId);
                if (p_protein_entry == null) {
                    log.info("ERROR: No entry for with transcript ID: " + transcriptId);
                    continue;
                }
//                protein_coordinates = new Coordinates();
                prevProteinCoordinates = new Coordinates();
                prevProteinCoordinates.setCterm(Offset.off3);
                prevProteinCoordinates.setNterm(Offset.off3);
                prevProteinCoordinates.setStart(0);
                prevProteinCoordinates.setEnd(0);
                coordinatesMap = new ArrayList<>();
            } else if (is_cds(tokens)) {
                GenomeCoordinates genCoord = Utils.extract_coordinates_from_gtf_line(tokens);
                genCoord.setTranscriptid(transcriptId);
                genCoord.setExonid(GeneEntry.extract_exon_id(m_line));
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

                coordinatesMap.add(new Pair<>(proteinCoordinates, genCoord));
            }
        }
        if (p_protein_entry != null) {
            p_protein_entry.set_coordinate_map(coordinatesMap);
        }
        close();
        return assem;
    }
}