package bigbio.pgatk.jpogo.common;

import bigbio.pgatk.jpogo.CoordinateWrapper;
import bigbio.pgatk.jpogo.GeneEntry;
import bigbio.pgatk.jpogo.MappedPeptides;
import bigbio.pgatk.jpogo.ProteinEntry;
import javafx.util.Pair;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GTFParser {

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
        return m_reader.ready();
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
        ArrayList<Pair<Coordinates, GenomeCoordinates>> coordinates_map = new ArrayList<>();
//        Coordinates protein_coordinates = new Coordinates();
        Coordinates prev_proteint_coordinates = new Coordinates();
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
                    p_protein_entry.set_coordinate_map(coordinates_map);
                }
                p_protein_entry = coordwrapper.lookup_entry(transcriptId);
                if (p_protein_entry == null) {
                    System.out.println("ERROR: No entry for with transcript ID: " + transcriptId);
                    continue;
                }
//                protein_coordinates = new Coordinates();
                prev_proteint_coordinates = new Coordinates();
                prev_proteint_coordinates.Cterm = Offset.off3;
                prev_proteint_coordinates.Nterm = Offset.off3;
                prev_proteint_coordinates.start = 0;
                prev_proteint_coordinates.end = 0;
                coordinates_map = new ArrayList<>();
            } else if (is_cds(tokens)) {
                GenomeCoordinates genCoord = Utils.extract_coordinates_from_gtf_line(tokens);
                genCoord.setTranscriptid(transcriptId);
                genCoord.setExonid(GeneEntry.extract_exon_id(m_line));
                Coordinates protein_coordinates = new Coordinates();
                // get nterm from prev exon
                if (genCoord.getFrame() != Frame.unknown) {
                    protein_coordinates.Nterm = Offset.forValue(genCoord.getFrame().getValue());
                } else {
                    if (prev_proteint_coordinates.Cterm != Offset.off3) {
                        protein_coordinates.Nterm = Offset.forValue(3 - prev_proteint_coordinates.Cterm.getValue());
                    } else {
                        protein_coordinates.Nterm = Offset.off3;
                    }
                }

                int length = 0;

                if (is_first_strand(tokens)) {
                    length = genCoord.end - genCoord.start + 1;
                } else if (!is_first_strand(tokens)) {
                    length = genCoord.start - genCoord.end + 1;
                }

                // calc cterm
                if (length % 3 == 0) {
                    if (protein_coordinates.Nterm != Offset.off3) {
                        protein_coordinates.Cterm = Offset.forValue(3 - protein_coordinates.Nterm.getValue());
                    } else {
                        protein_coordinates.Cterm = Offset.off3;
                    }
                } else if (length % 3 == 2) {
                    if (protein_coordinates.Nterm == Offset.off3) {
                        protein_coordinates.Cterm = Offset.off2;
                    } else if (protein_coordinates.Nterm == Offset.off2) {
                        protein_coordinates.Cterm = Offset.off3;
                    } else if (protein_coordinates.Nterm == Offset.off1) {
                        protein_coordinates.Cterm = Offset.off1;
                    }
                } else if (length % 3 == 1) {
                    if (protein_coordinates.Nterm == Offset.off3) {
                        protein_coordinates.Cterm = Offset.off1;
                    } else if (protein_coordinates.Nterm == Offset.off1) {
                        protein_coordinates.Cterm = Offset.off3;
                    } else if (protein_coordinates.Nterm == Offset.off2) {
                        protein_coordinates.Cterm = Offset.off2;
                    }
                }

                // calc protein coordinates
                if (protein_coordinates.Nterm != Offset.off3) {
                    protein_coordinates.start = prev_proteint_coordinates.end;
                } else {
                    if (prev_proteint_coordinates.end == 0 && coordinates_map.isEmpty()) {
                        protein_coordinates.start = 0;
                    } else {
                        protein_coordinates.start = prev_proteint_coordinates.end + 1;
                    }
                }

                int offsets = 0;
                if (protein_coordinates.Nterm != Offset.off3) {
                    offsets = offsets + protein_coordinates.Nterm.getValue();
                }

                if (is_first_strand(tokens)) {
                    length = genCoord.end - genCoord.start + 1 - offsets;
                } else if (!is_first_strand(tokens)) {
                    length = genCoord.start - genCoord.end + 1 - offsets;
                }

                int peplength = length / 3;

                int pepend = protein_coordinates.start + peplength - 1;
                if (protein_coordinates.Cterm != Offset.off3) {
                    pepend = pepend + 1;
                }
                if (protein_coordinates.Nterm != Offset.off3) {
                    pepend = pepend + 1;
                }

                protein_coordinates.end = pepend;

                prev_proteint_coordinates = protein_coordinates;

                coordinates_map.add(new Pair<>(protein_coordinates, genCoord));
            }
        }
        if (p_protein_entry != null) {
            p_protein_entry.set_coordinate_map(coordinates_map);
        }
        close();
        return assem;
    }
}