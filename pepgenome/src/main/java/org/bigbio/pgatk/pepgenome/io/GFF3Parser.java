package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.PepGenomeTool;
import org.bigbio.pgatk.pepgenome.common.*;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.mortbay.log.Log;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// GFF3 Parser
public class GFF3Parser extends GenomeAnnotationParser {

    public static GFF3Parser instance;

    private GFF3Parser() {
    }

    // For use with exco mode
    private static int translationOffset = 0;

    // GFF3 Patterns
    private static final Pattern GFFIDPATTERN = Pattern.compile("ID=([^;]*)"); // ID tag
    private static final Pattern GFFPARENTPATTERN = Pattern.compile("Parent=([^;]*)"); // Parent tag

    // Singleton get_instance method.
    public static GFF3Parser get_instance() {
        if (instance == null) {
            instance = new GFF3Parser();
        }
        return instance;
    }

    // Reads a gff3 file and parses it into CoordinateWrapper and MappedPeptides.
    public final Assembly read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping) throws Exception {
        if (!open(file)) {
            throw new IllegalStateException("Problem in reading GFF3 file");
        }

        String exonID = "";
        ProteinEntry proteinEntry = null;
        ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinatesMap = new ArrayList<>();
        Coordinates prevProteinCoordinates = new Coordinates();
        Assembly assem = Assembly.none;
        ArrayList<String> tokens;
        int remainingProteinLength = 0;

        // Hash map of transcript IDs to gene IDs so that gene ID may be retrieved for exons via transcripts.
        HashMap<String, String> idMap = new HashMap<>();
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }

            // Convert GFF3 line into 9 tokens
            tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t")));

            // GENE
            if (is_next_gene(tokens)) {
                Assembly assemtemp = mapping.add_gene_from_annotation(line);
                if (assem == Assembly.none) {
                    if (assemtemp == Assembly.patchhaploscaff) {
                        assem = assemtemp;
                    }
                }
            }

            // TRANSCRIPT
            String transcriptId = "";
            if (is_next_transcript(tokens)) {

                String geneId = extract_id(line, GFFPARENTPATTERN);
                transcriptId = extract_id(line, GFFIDPATTERN);
                idMap.put(transcriptId,geneId); // Places transcript id and gene id into hash map.

                exonID = "";
                mapping.add_transcript_id_to_gene(line);
                if (proteinEntry != null) {
                    proteinEntry.set_coordinate_map(coordinatesMap);
                }
                proteinEntry = coordwrapper.lookup_entry(transcriptId);
                if (proteinEntry == null) {
                    Log.info("ERROR: No entry for transcript ID: " + transcriptId);
                    continue;
                }

                remainingProteinLength = proteinEntry.get_sequence().length();
                prevProteinCoordinates = new Coordinates();
                prevProteinCoordinates.setCterm(Offset.off3);
                prevProteinCoordinates.setNterm(Offset.off3);
                prevProteinCoordinates.setStart(0);
                prevProteinCoordinates.setEnd(0);
                coordinatesMap = new ArrayList<>();

                // Retrieve transcript's translation offset based on its ID.
                if (PepGenomeTool.m_translation_offset_map.get(transcriptId) == null) {
                    translationOffset = 0;
                }
                else {
                    translationOffset = PepGenomeTool.m_translation_offset_map.get(transcriptId);
                }

                // EXON
            } else if (is_exon(tokens)) {

                if (PepGenomeTool.useExonCoords) {

                    // CDS
                    // Check exon has parent ID matching the last transcript ID - Ensure offset is being applied correctly.  Possibly unnecessary.
                    if (extract_id(line, GFFPARENTPATTERN).equals(transcriptId)) {

                        // Extract genomic coordinates from exon line
                        GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens);

                        //  Determining length of exon
                        int length = 0;
                        if (is_first_strand(tokens)) {
                            length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                        } else if (!is_first_strand(tokens)) {
                            length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                        }

                        // Check length of exon against offset.  Adjust coords accordingly.
                        // Check remaining protein length, 0 or below indicates all following exons are untranslated.

                        // UNTRANSLATED EXON
                        if (translationOffset > length || remainingProteinLength <= 0 ) {
                            // Offset value greater than total exon length - Exon  is completely untranslated

                            // Adjust offset
                            translationOffset = translationOffset - length;
                            System.out.println(translationOffset);

                            continue;


                            // PARTIALLY TRANSLATED EXON
                        } else if (translationOffset > 0) {
                            // Offset value not greater than exon length, but still greater than 0 - Exon is partially translated.
                            if (is_first_strand(tokens)) {
                                genomeCoordinates.setStart(genomeCoordinates.getStart() + translationOffset);
                            } else {
                                genomeCoordinates.setEnd(genomeCoordinates.getEnd() - translationOffset);
                            }

                            // Adjust offset
                            translationOffset = 0;
                        }

                        // PARTIALLY/FULLY TRANSLATED EXON
                        genomeCoordinates.setTranscriptid(transcriptId);  // Using previous transcript's ID
                        exonID = extract_exon_id(line);  // Extracted from current exon line
                        genomeCoordinates.setExonid(exonID);
                        Coordinates proteinCoordinates = new Coordinates();

                        // Get N term from previous exon
                        if (genomeCoordinates.getFrame() != Frame.unknown) {
                            proteinCoordinates.setNterm(Offset.forValue(genomeCoordinates.getFrame().getValue()));
                        } else {
                            if (prevProteinCoordinates.getCterm() != Offset.off3) {
                                proteinCoordinates.setNterm(Offset.forValue(3 - prevProteinCoordinates.getCterm().getValue()));
                            } else {
                                proteinCoordinates.setNterm(Offset.off3);
                            }
                        }

                        // Calculating C term
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

                        // Calculate protein coordinates
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
                            length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1 - offsets;
                        } else if (!is_first_strand(tokens)) {
                            length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1 - offsets;
                        }

                        int peplength = length / 3;
                        remainingProteinLength = remainingProteinLength - peplength;

                        int pepend = proteinCoordinates.getStart() + peplength - 1;
                        if (proteinCoordinates.getCterm() != Offset.off3) {
                            ++pepend;
                        }
                        if (proteinCoordinates.getNterm() != Offset.off3) {
                            ++pepend;
                        }

                        proteinCoordinates.setEnd(pepend);
                        prevProteinCoordinates = proteinCoordinates;
                        coordinatesMap.add(new Tuple<>(proteinCoordinates, genomeCoordinates));

                    }

                }

            } else if (is_cds(tokens)) {


                GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens);
                genomeCoordinates.setTranscriptid(transcriptId);
                String tmp_exonID = extract_id(line,GFFPARENTPATTERN);

                if (tmp_exonID.equals("")){
                    tmp_exonID = exonID;
                }

                genomeCoordinates.setExonid(tmp_exonID);
                Coordinates proteinCoordinates = new Coordinates();

                // Get N term from previous exon
                if (genomeCoordinates.getFrame() != Frame.unknown) {
                    proteinCoordinates.setNterm(Offset.forValue(genomeCoordinates.getFrame().getValue()));
                } else {
                    if (prevProteinCoordinates.getCterm() != Offset.off3) {
                        proteinCoordinates.setNterm(Offset.forValue(3 - prevProteinCoordinates.getCterm().getValue()));
                    } else {
                        proteinCoordinates.setNterm(Offset.off3);
                    }
                }

                //  Determining length
                int length = 0;
                if (is_first_strand(tokens)) {
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                }

                // Calculating C term
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

                // Calculate protein coordinates
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
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1 - offsets;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1 - offsets;
                }

                int peplength = length / 3;

                int pepend = proteinCoordinates.getStart() + peplength - 1;
                if (proteinCoordinates.getCterm() != Offset.off3) {
                    ++pepend;
                }
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    ++pepend;
                }

                proteinCoordinates.setEnd(pepend);

                prevProteinCoordinates = proteinCoordinates;

                coordinatesMap.add(new Tuple<>(proteinCoordinates, genomeCoordinates));

            }

        }

        if (proteinEntry != null) {
            proteinEntry.set_coordinate_map(coordinatesMap);
        }
        close();
        return assem;
    }

    public static String extract_gene_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    public static String extract_transcript_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    public static String extract_exon_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    public static String extract_id(String gtfGeneLine, Pattern pattern) {
        String value = "";
        Matcher matcher = pattern.matcher(gtfGeneLine);
        if (matcher.find()) {
            value = matcher.group(1);
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
