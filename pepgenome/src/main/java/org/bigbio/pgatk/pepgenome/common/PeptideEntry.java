package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;

import org.apache.commons.lang3.StringUtils;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.*;

/**
 * A peptide entry contains information about a single peptide:
 *  - all its occurrences in the genome,
 *  - the associated gene, if it is unique to transcript or gene,
 *  - the tissues it appears in,
 *  - as well as the first and last appearance in the genome.
 *
 * @author ypriverol
 */

public class PeptideEntry implements Comparable<PeptideEntry>, Serializable {

    private static final long serialVersionUID = -3518477125981077637L;

    //holds the peptide sequence.
    private String pSequence;

    //number of found transcripts.
    private int numTranscripts;

    //is the peptide unique to one gene?
    private boolean geneUnique;

    //is the peptide unique to one transcript?
    private boolean transcriptUnique;

    //set of coordinates.
    private TreeSet<PeptideCoordinates> pepCoordinates = new TreeSet<>(new PeptidecoordsPCompare());

    //set of tissues.
    private Map<String, Tuple<ArrayList<Integer>, ArrayList<Double>>> tissueTags = new TreeMap<>();

    //lowest startcoord of peptides that share the same sequence.
    private int startCoord;

    //highest endcoord of peptides that share the same sequence.
    private int endCoord;

    //pointer to the associated gene
    private GeneEntry associatedGene;

    //this map holds the different PTMs for a peptide.
    private Map<String, Map<String, PTMEntry>> pepForms = new TreeMap<>();

    private Set<String> transcriptIds = new TreeSet<>();

    private Set<String> exonIds = new TreeSet<>();

    @Override
    public int compareTo(PeptideEntry o) {
        if (lessThan(o)) {
            return -1;
        } else if (o.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    public PeptideEntry(GeneEntry associatedgene) {
        this.pSequence = "";
        this.numTranscripts = 0;
        this.geneUnique = true;
        this.transcriptUnique = true;
        this.startCoord = 0;
        this.endCoord = 0;
        this.associatedGene = associatedgene;
    }

    //generates a map that holds all PTMs in the given sequence.
    private static Map<String, PTMEntry> ptmSet(String sequence) {
        Map<String, PTMEntry> map = new TreeMap<>();

        int startB = sequence.indexOf("(");
        int endB = sequence.indexOf(")");
        String name;
        while (startB != -1 && endB != -1) {
            name = sequence.substring(startB + 1, endB);
            String substr = "(" + name + ")";
            sequence = StringUtils.replaceOnce(sequence, substr, "");
            if (startB != 0) {
                startB = startB - 1;
            }
            if (!name.equals("")) {
                if (!map.containsKey(name)) {
                    map.put(name, new PTMEntry(name, startB, startB));
                } else {
                    map.get(name).add_coord(startB);
                }
            }
            //}
            startB = sequence.indexOf("(");
            endB = sequence.indexOf(")");
        }
        return map;
    }

    //comparator to check if one PeptideEntry is smaller than another one.
    //returns true if the startcoordinate are smaller than rhs' startcoordinate
    //otherwise returns true if endcoordinate is smaller than rhs' endcoordinate
    //otherwise returns true if lhs' sequence is lexicographically lesser than rhs'
    //otherwise returns false.
    public boolean lessThan(PeptideEntry rhs) {
        return (startCoord < rhs.startCoord) || (startCoord == rhs.startCoord && endCoord < rhs.endCoord)
                || (startCoord == rhs.startCoord && endCoord == rhs.endCoord && pSequence.compareTo(rhs.pSequence) < 0);
    }

    public final OutputStream to_gtf(String source, OutputStream os) throws Exception {
        return to_gtf(source, os, true);
    }

    public final OutputStream to_gtf(String source, boolean chrincluded) throws Exception {
        return to_gtf(source, System.out, chrincluded);
    }

    public final OutputStream to_gtf(String source) throws Exception {
        return to_gtf(source, System.out, true);
    }

    //generates a string in the gtf format and writes it to the specified ostream.
    public final OutputStream to_gtf(String source, OutputStream os, boolean chrincluded) throws Exception {
        String sequence_add = "";
        int count = 0;
        for (PeptideCoordinates coord : pepCoordinates) {
            count += 1;
            if (pepCoordinates.size() > 1) {
                sequence_add = "." + count;
            }
            if (count > 1) {
                os.write("\n".getBytes());
            }

            os.write(Utils.coordinates_to_gtf_string(coord.get_transcript_coordinates(), "transcript", false, source).getBytes());

            os.write(("gene_id \"" + associatedGene.get_id() + "\"; transcript_id \"" +
                    associatedGene.get_id() + "." + pSequence + sequence_add + "\"; gene_type \"" +
                    associatedGene.get_type() + "\"; gene_status \"" + associatedGene.get_status() +
                    "\"; gene_name \"" + associatedGene.get_name()).getBytes());

            os.write(("\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"" +
                    associatedGene.get_name() + "." + pSequence + sequence_add + "\";").getBytes());

            os.write((" tag \"Transcripts:" + numTranscripts + "\";").getBytes());

            os.write((" tag \"TranscriptIDs:" + transcriptids_to_string() + "\";"
                    + " tag \"ExonIDs:" + exonids_to_string() + "\";").getBytes());

            for (Map.Entry<String, Tuple<ArrayList<Integer>, ArrayList<Double>>> current : tissueTags.entrySet()) {
                os.write((" tag \"" + current.getKey() + ":").getBytes());
                StringBuilder ss_tissue = new StringBuilder();
                for (int i_tissue = 0; i_tissue < current.getValue().getKey().size(); ++i_tissue) {
                    if (i_tissue > 0) {
                        os.write("/".getBytes());
                        ss_tissue.append("/");
                    }
                    os.write(String.valueOf(current.getValue().getKey().get(i_tissue)).getBytes());
                    ss_tissue.append(current.getValue().getValue().get(i_tissue));
                }
                os.write((" sig PSMs" + " " + ss_tissue + " Quant" + "\";").getBytes());
            }

            ArrayList<GenomeCoordinates> exon_coordinates = coord.get_exon_coordinates();
            int exon_count = 0;
            String exon_add = "";
            for (int exit = 0; exit < exon_coordinates.size(); ++exit) {
                exon_count += 1;
                if (exon_coordinates.size() > 1) {
                    exon_add = "." + exon_count;
                }

                os.write("\n".getBytes());
                os.write(Utils.coordinates_to_gtf_string(exon_coordinates.get(exit), "exon", false, source).getBytes());

                os.write(("gene_id \"" + associatedGene.get_id() + "\"; transcript_id \"" + associatedGene.get_id()
                        + "." + pSequence + sequence_add + "\"; gene_type \"" + associatedGene.get_type()
                        + "\"; gene_status \"" + associatedGene.get_status() + "\"; gene_name \"" + associatedGene.get_name()).getBytes());

                os.write(("\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \""
                        + associatedGene.get_name() + "." + pSequence + sequence_add + "\";").getBytes());

                os.write((" exon_number " + exon_count + "; exon_id \"" + associatedGene.get_id() +
                        "." + pSequence + sequence_add + exon_add + "\";").getBytes());

                os.write((" tag \"TranscriptIDs:" + transcriptids_to_string() + "\";" + " tag \"ExonIDs:" + exonids_to_string() + "\";").getBytes());
            }
        }
        return os;
    }

    public final OutputStream to_bed(boolean noptm, boolean chrincluded) throws Exception {
        return to_bed(System.out, noptm, chrincluded);
    }

    public final OutputStream to_bed(OutputStream os, boolean noptm) throws Exception {
        return to_bed(os, noptm, true);
    }

    public final OutputStream to_bed(OutputStream os) throws Exception {
        return to_bed(os, false, true);
    }

    public final OutputStream to_bed() throws Exception {
        return to_bed(System.out, false, true);
    }

    //generates a bed line and writes it to the specified ostream.
    public final OutputStream to_bed(OutputStream os, boolean noptm, boolean chrincluded) throws Exception {
        for (PeptideCoordinates coord : pepCoordinates) {
            os.write(Utils.coordinates_to_bed_string(coord.get_transcript_coordinates(), pSequence).getBytes());
            //std::cout << coordinates_to_bed_string((*it)->get_transcript_coordinates(), pSequence) << std::endl;

            StringBuilder exon_starts = new StringBuilder();
            StringBuilder exon_lengths = new StringBuilder();
            ArrayList<GenomeCoordinates> exon_coordinates = coord.get_exon_coordinates();
            int exon_count = 0;
            for (GenomeCoordinates exon_coordinate : exon_coordinates) {
                exon_count += 1;
                if (exon_count > 1) {
                    exon_starts.append(",");
                    exon_lengths.append(",");
                }
                int exon_start = exon_coordinate.getStart() - coord.get_transcript_coordinates().getStart();
                int exon_length = exon_coordinate.getEnd() - exon_coordinate.getStart() + 1;
                exon_starts.append(exon_start);
                exon_lengths.append(exon_length);
            }
            String colour = "128,128,128";

            if (!noptm) {
                if (geneUnique && !transcriptUnique) {
                    colour = "0,0,0";
                } else if (geneUnique && transcriptUnique) {
                    colour = "204,0,0";
                }
            }

            os.write((colour + "\t" + exon_count + "\t" + exon_lengths + "\t" + exon_starts + "\n").getBytes());
        }
        return os;
    }

    public final OutputStream to_gct(String geneID, ArrayList<String> tissuelist, OutputStream os) throws Exception {
        return to_gct(geneID, tissuelist, os, true);
    }

    public final OutputStream to_gct(String geneID, ArrayList<String> tissuelist) throws Exception {
        return to_gct(geneID, tissuelist, System.out, true);
    }

    public final OutputStream to_gct(String geneID, ArrayList<String> tissuelist, boolean chrincluded) throws Exception {
        return to_gct(geneID, tissuelist, System.out, chrincluded);
    }

    //generates a gct line and writes it to the specified ostream.
    public final OutputStream to_gct(String geneID, ArrayList<String> tissuelist, OutputStream os, boolean chrincluded) throws Exception {
        String sequence_add = "";
        int count = 0;
        for (PeptideCoordinates coord : pepCoordinates) {
            count += 1;
            if (pepCoordinates.size() > 1) {
                sequence_add = "." + count;
            }

            os.write((geneID + "." + pSequence + sequence_add + "\t\"" + geneID + "|@").getBytes());

            ArrayList<GenomeCoordinates> exoncoords = coord.get_exon_coordinates();
            os.write((Utils.coordinates_to_gct_string(exoncoords) + "|\"").getBytes());

            os.write((tissue_quant_to_string(tissuelist) + "\n").getBytes());
        }
        return os;
    }

    //this function is used to generate a date for a gct line.
    private String tissue_quant_to_string(ArrayList<String> tissuelist) {
        StringBuilder ss = new StringBuilder();
        for (int i = 1; i < tissuelist.size(); ++i) {
            ss.append("\t");
            Tuple<ArrayList<Integer>, ArrayList<Double>> pair = tissueTags.get(tissuelist.get(i));
            if (pair != null) {
                ArrayList<Double> doubles = pair.getValue();
                double sum = doubles.stream().mapToDouble(Double::doubleValue).sum();
                double mean = sum / doubles.size();
                ss.append(mean);
            }
        }
        return ss.toString();
    }

    public final OutputStream to_ptmbed(OutputStream os) throws Exception {
        return to_ptmbed(os, true);
    }

    public final OutputStream to_ptmbed(boolean chrincluded) throws Exception {
        return to_ptmbed(System.out, chrincluded);
    }

    public final OutputStream to_ptmbed() throws Exception {
        return to_ptmbed(System.out, true);
    }

    //generates a bed line with ptms and writes it to the specified ostream.
    public final OutputStream to_ptmbed(OutputStream os, boolean chrincluded) throws Exception {
        for (Map.Entry<String, Map<String, PTMEntry>> ptm_it : pepForms.entrySet()) {
            for (Map.Entry<String, PTMEntry> ptm_single_it : ptm_it.getValue().entrySet()) {
                List<Tuple<PeptideCoordinates, GenomeCoordinates>> coord = ptm_single_it.getValue().get_genome_coordinates();
                for (Tuple<PeptideCoordinates, GenomeCoordinates> coord_it : coord) {
                    String bed_string = Utils.coordinates_to_short_bed_string(coord_it.getKey().get_transcript_coordinates(), ptm_it.getKey());
                    ArrayList<GenomeCoordinates> exon_coords = coord_it.getKey().get_exon_coordinates();

                    StringBuilder exon_starts = new StringBuilder();
                    StringBuilder exon_lengths = new StringBuilder();
                    int exon_count = 0;
                    for (GenomeCoordinates exon_coord : exon_coords) {
                        exon_count += 1;
                        if (exon_count > 1) {
                            exon_starts.append(",");
                            exon_lengths.append(",");
                        }
                        int exon_start = exon_coord.getStart() - coord_it.getKey().get_transcript_coordinates().getStart();
                        int exon_length = exon_coord.getEnd() - exon_coord.getStart() + 1;

                        exon_starts.append(exon_start);
                        exon_lengths.append(exon_length);
                    }
                    String colour = EnumStringMapper.ptmToColour(ptm_single_it.getKey());
                    os.write((bed_string + (coord_it.getValue().getStart() - 1) + "\t" + coord_it.getValue().getEnd()
                            + "\t" + "\t" + colour + "\t" + exon_count + "\t" + exon_lengths + "\t"
                            + exon_starts + "\n").getBytes());
                }
            }
        }
        return os;
    }

    //adds new tissuetags if they havent existed before.
    public final void add_tags(String tag, int sigPSMs, double quant) {

        Tuple<ArrayList<Integer>, ArrayList<Double>> pair = tissueTags.computeIfAbsent(tag, k -> new Tuple<>(new ArrayList<>(), new ArrayList<>()));

        pair.getKey().add(sigPSMs);
        pair.getValue().add(quant);
    }

    //adds ptms if a sequence matches another with ptms.
    public final void add_ptm(String ptmsequence) {
        //adds the modfication to the peptide entry
        pepForms.put(ptmsequence, ptmSet(ptmsequence));
        for (Map.Entry<String, PTMEntry> ptm_it : pepForms.get(ptmsequence).entrySet()) {
            for (PeptideCoordinates peptide_coordinates_it : pepCoordinates) {
                ArrayList<GenomeCoordinates> genomic_coordinates = peptide_coordinates_it.find_coordinates(ptm_it.getValue().get_range());
                GenomeCoordinates ptm_coordinates = new GenomeCoordinates();
                for (int j = 0; j < genomic_coordinates.size(); ++j) {
                    GenomeCoordinates genomeCoordinates = new GenomeCoordinates(genomic_coordinates.get(j));
                    if (j == 0) {
                        ptm_coordinates = genomeCoordinates;
                    } else {
                        if (ptm_coordinates.getStart() > genomeCoordinates.getStart()) {
                            ptm_coordinates.setStart(genomeCoordinates.getStart());
                        }
                        if (ptm_coordinates.getEnd() < genomeCoordinates.getEnd()) {
                            ptm_coordinates.setEnd(genomeCoordinates.getEnd());
                        }
                    }
                }
                ptm_it.getValue().add_genome_coordinates(peptide_coordinates_it, ptm_coordinates);
            }
        }
    }

    //this function creates a coordinate_map_type and works similar to a CoordinateMapTypeConstructor.
    private static ArrayList<Tuple<Coordinates, GenomeCoordinates>> create_coordinate_map_type(ArrayList<GenomeCoordinates> genomecoords) {
        ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordMap = new ArrayList<>();
        Coordinates prevProtCoord = new Coordinates();
        prevProtCoord.setCterm(Offset.off3);
        prevProtCoord.setNterm(Offset.off3);
        prevProtCoord.setStart(0);
        prevProtCoord.setEnd(0);

        for (GenomeCoordinates genomeCoordinatesItr : genomecoords) {
            Coordinates protCoord = new Coordinates();
            GenomeCoordinates genomeCoordinates = new GenomeCoordinates(genomeCoordinatesItr);
            if (genomeCoordinates.getStrand() == Strand.rev) {
                int start = genomeCoordinates.getEnd();
                int end = genomeCoordinates.getStart();
                genomeCoordinates.setStart(start);
                genomeCoordinates.setEnd(end);
            }

            if (prevProtCoord.getCterm() != Offset.off3) {
                protCoord.setNterm(Offset.forValue(3 - prevProtCoord.getCterm().getValue()));
            } else {
                protCoord.setNterm(Offset.off3);
            }

            int length;

            if (genomeCoordinates.getStrand() == Strand.fwd) {
                length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
            } else {
                length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
            }
            // calc cterm
            if (length % 3 == 0) {
                if (protCoord.getNterm() != Offset.off3) {
                    protCoord.setCterm(Offset.forValue(3 - protCoord.getNterm().getValue()));
                } else {
                    protCoord.setCterm(Offset.off3);
                }
            } else if (length % 3 == 2) {
                if (protCoord.getNterm() == Offset.off3) {
                    protCoord.setCterm(Offset.off2);
                } else if (protCoord.getNterm() == Offset.off2) {
                    protCoord.setCterm(Offset.off3);
                } else if (protCoord.getNterm() == Offset.off1) {
                    protCoord.setCterm(Offset.off1);
                }
            } else if (length % 3 == 1) {
                if (protCoord.getNterm() == Offset.off3) {
                    protCoord.setCterm(Offset.off1);
                } else if (protCoord.getNterm() == Offset.off1) {
                    protCoord.setCterm(Offset.off3);
                } else if (protCoord.getNterm() == Offset.off2) {
                    protCoord.setCterm(Offset.off2);
                }
            }

            // calc protein coordinates
            if (protCoord.getNterm() != Offset.off3) {
                protCoord.setStart(prevProtCoord.getEnd());
            } else {
                if (prevProtCoord.getEnd() == 0 && coordMap.isEmpty()) {
                    protCoord.setStart(0);
                } else {
                    protCoord.setStart(prevProtCoord.getEnd() + 1);
                }
            }

            int offsets = 0;
            if (protCoord.getNterm() != Offset.off3) {
                offsets = offsets + protCoord.getNterm().getValue();
            }

            if (genomeCoordinates.getStrand() == Strand.fwd) {
                length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1 - offsets;
            } else {
                length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1 - offsets;
            }

            int pep_length = length / 3;

            int pep_end = protCoord.getStart() + pep_length - 1;
            if (protCoord.getCterm() != Offset.off3) {
                pep_end = pep_end + 1;
            }
            if (protCoord.getNterm() != Offset.off3) {
                pep_end = pep_end + 1;
            }

            protCoord.setEnd(pep_end);

            prevProtCoord = protCoord;

            coordMap.add(new Tuple<>(protCoord, genomeCoordinates));
        }
        return coordMap;
    }

    //adds a peptide. this function is used if this specific peptide has not yet been found. it will also find the peptides genomic coordinates.
    public final void add_peptide(CoordinateWrapper coordwrapper, String sequence, String ptmSequence, String tag, int sigPSMs, TranscriptsT transcripts, int genes, FileOutputStream ofstream, double quant) {
        if (pSequence.isEmpty()) {
            pSequence = sequence;
        }
        if (numTranscripts == 0) {
            numTranscripts = transcripts.getM_entries().size();
        }

        if (genes > 1) {
            geneUnique = false;
        }
        if (transcripts.getM_entries().size() > 1) {
            transcriptUnique = false;
        }

        add_tags(tag, sigPSMs, quant);
        //iterate all found transcripts.
        for (Map.Entry<String, ArrayList<PositionMismatchT>> it : transcripts.getM_entries().entrySet()) {
            // find all genomic coordinates
            ArrayList<ArrayList<GenomeCoordinates>> genomic_coordinates = coordwrapper.lookup_entry(it.getKey()).find_coordinates(sequence.length(), it.getValue());
            int CDS_annotation_correct = coordwrapper.lookup_entry(it.getKey()).get_cds_annotation_correct();
            //iterate all genomic coordinates.
            for (ArrayList<GenomeCoordinates> genomic_coordinate : genomic_coordinates) {
                //creates PetideCoordinates
                PeptideCoordinates pep_coord = new PeptideCoordinates(create_coordinate_map_type(genomic_coordinate), CDS_annotation_correct);
                TreeSet<String> transcriptids = pep_coord.get_trasncript_ids();
                transcriptIds.addAll(transcriptids);
                TreeSet<String> exonids = pep_coord.get_exon_ids();
                exonIds.addAll(exonids);
                if (genomic_coordinate.size() != 0) {
                    //and saves them.
                    pepCoordinates.add(pep_coord);
                    //sets start and end coord of the peptide entry to min/max values. these are used for comparing PeptideEntry objects.
                    if (startCoord == 0) {
                        startCoord = pep_coord.get_transcript_coordinates().getStart();
                    }
                    if (endCoord == 0) {
                        endCoord = pep_coord.get_transcript_coordinates().getEnd();
                    }
                    if (pep_coord.get_transcript_coordinates().getStart() < startCoord) {
                        startCoord = pep_coord.get_transcript_coordinates().getStart();
                    }
                    if (pep_coord.get_transcript_coordinates().getEnd() > endCoord) {
                        endCoord = pep_coord.get_transcript_coordinates().getEnd();
                    }
                }
            }
        }
        add_ptm(ptmSequence);
    }

    //adds a peptide. this function is used if a peptide has appeared before.
    public final void add_peptide(String ptmsequence, String tag, int sigPSMs, double quant) {
        if (!pSequence.equals(ptmsequence)) {
            add_ptm(ptmsequence);
        }
        add_tags(tag, sigPSMs, quant);
    }

    private String transcriptids_to_string() {
        StringBuilder transcriptid_string = new StringBuilder();
        for (String m_transcriptid : transcriptIds) {
            if (!transcriptid_string.toString().isEmpty()) {
                transcriptid_string.append("|");
            }
            transcriptid_string.append(m_transcriptid);
        }
        return transcriptid_string.toString();
    }

    private String exonids_to_string() {
        StringBuilder exonid_string = new StringBuilder();
        for (String it : exonIds) {
            if (!exonid_string.toString().isEmpty()) {
                exonid_string.append("|");
            }
            exonid_string.append(it);
        }
        return exonid_string.toString();
    }

    //checks whether the peptide entry has a modified form. Returns true if pepForms only contains peptide without PTM (size==1)
    public final boolean noPTM() {
        return pepForms.containsKey(pSequence);
    }
}