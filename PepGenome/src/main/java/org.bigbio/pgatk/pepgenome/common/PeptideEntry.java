package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import javafx.util.Pair;
import org.apache.commons.lang3.StringUtils;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.*;

//a peptide entry contains information about a single peptide:
//all its occurences in the genome, 
//the associated gene, if it is unique to transcript or gene, 
//the tissues it appears in, 
//as well as the first and last appearance in the genome.
public class PeptideEntry implements Comparable<PeptideEntry>, Serializable {

    private static final long serialVersionUID = -3518477125981077637L;
    //holds the peptide sequence.
    private String m_sequence;
    //number of found transcripts.
    private int m_numtranscripts;
    //is the peptide unique to one gene?
    private boolean m_geneunique;
    //is the peptide unique to one transcript?
    private boolean m_transcriptunique;
    //set of coordinates.
    private TreeSet<PeptideCoordinates> m_pepcoordinates = new TreeSet<>(new PeptidecoordsPCompare());
    //set of tissues.
    private Map<String, Pair<ArrayList<Integer>, ArrayList<Double>>> m_tissuetags = new TreeMap<>();
    //lowest startcoord of peptides that share the same sequence.
    private int m_startcoord;
    //highest endcoord of peptides that share the same sequence.
    private int m_endcoord;
    //pointer to the associated gene
    private GeneEntry m_associatedgene;
    //this map holds the different PTMs for a peptide.
    private Map<String, Map<String, PTMEntry>> m_pepforms = new TreeMap<>();

    private Set<String> m_transcriptids = new TreeSet<>();
    private Set<String> m_exonids = new TreeSet<>();

    @Override
    public int compareTo(PeptideEntry otherInstance) {
        if (lessThan(otherInstance)) {
            return -1;
        } else if (otherInstance.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    public PeptideEntry(GeneEntry associatedgene) {
        this.m_sequence = "";
        this.m_numtranscripts = 0;
        this.m_geneunique = true;
        this.m_transcriptunique = true;
//        this.m_pepcoordinates = new TreeSet<>(new PeptidecoordsPCompare());
//        this.m_tissuetags = new TreeMap<>();
        this.m_startcoord = 0;
        this.m_endcoord = 0;
        this.m_associatedgene = associatedgene;
//        this.m_pepforms = = new TreeMap<>();
//        this.m_transcriptids = new TreeSet<>();
//        this.m_exonids = new TreeSet<>();
    }

    //generates a map that holds all PTMs in the given sequence.
    private static Map<String, PTMEntry> ptm_set(String sequence) {
        Map<String, PTMEntry> map = new TreeMap<>();

        int start_b = sequence.indexOf("(");
        int end_b = sequence.indexOf(")");
        String name;
        while (start_b != -1 && end_b != -1) {
            name = sequence.substring(start_b + 1, end_b);
            String substr = "(" + name + ")";
            sequence = StringUtils.replaceOnce(sequence, substr, "");
            if (start_b != 0) {
                start_b = start_b - 1;
            }
            if (!name.equals("")) {
                if (!map.containsKey(name)) {
                    map.put(name, new PTMEntry(name, start_b, start_b));
                } else {
                    map.get(name).add_coord(start_b);
                }
            }
            //}
            start_b = sequence.indexOf("(");
            end_b = sequence.indexOf(")");
        }
        return map;
    }

    //comparator to check if one PeptideEntry is smaller than another one.
    //returns true if the startcoordinate are smaller than rhs' startcoordinate
    //otherwise returns true if endcoordinate is smaller than rhs' endcoordinate
    //otherwise returns true if lhs' sequence is lexicographically lesser than rhs'
    //otherwise returns false.
    public boolean lessThan(PeptideEntry rhs) {
        return (m_startcoord < rhs.m_startcoord) || (m_startcoord == rhs.m_startcoord && m_endcoord < rhs.m_endcoord)
                || (m_startcoord == rhs.m_startcoord && m_endcoord == rhs.m_endcoord && m_sequence.compareTo(rhs.m_sequence) < 0);
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
        for (PeptideCoordinates coord : m_pepcoordinates) {
            count += 1;
            if (m_pepcoordinates.size() > 1) {
                sequence_add = "." + count;
            }
            if (count > 1) {
                os.write("\n".getBytes());
            }

            os.write(Utils.coordinates_to_gtf_string(coord.get_transcript_coordinates(), "transcript", false, source).getBytes());

            os.write(("gene_id \"" + m_associatedgene.get_id() + "\"; transcript_id \"" +
                    m_associatedgene.get_id() + "." + m_sequence + sequence_add + "\"; gene_type \"" +
                    m_associatedgene.get_type() + "\"; gene_status \"" + m_associatedgene.get_status() +
                    "\"; gene_name \"" + m_associatedgene.get_name()).getBytes());

            os.write(("\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"" +
                    m_associatedgene.get_name() + "." + m_sequence + sequence_add + "\";").getBytes());

            os.write((" tag \"Transcripts:" + m_numtranscripts + "\";").getBytes());

            os.write((" tag \"TranscriptIDs:" + transcriptids_to_string() + "\";"
                    + " tag \"ExonIDs:" + exonids_to_string() + "\";").getBytes());

            for (Map.Entry<String, Pair<ArrayList<Integer>, ArrayList<Double>>> current : m_tissuetags.entrySet()) {
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
                os.write((" sig PSMs" + " " + ss_tissue.toString() + " Quant" + "\";").getBytes());
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

                os.write(("gene_id \"" + m_associatedgene.get_id() + "\"; transcript_id \"" + m_associatedgene.get_id()
                        + "." + m_sequence + sequence_add + "\"; gene_type \"" + m_associatedgene.get_type()
                        + "\"; gene_status \"" + m_associatedgene.get_status() + "\"; gene_name \"" + m_associatedgene.get_name()).getBytes());

                os.write(("\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \""
                        + m_associatedgene.get_name() + "." + m_sequence + sequence_add + "\";").getBytes());

                os.write((" exon_number " + exon_count + "; exon_id \"" + m_associatedgene.get_id() +
                        "." + m_sequence + sequence_add + exon_add + "\";").getBytes());

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
        for (PeptideCoordinates coord : m_pepcoordinates) {
            os.write(Utils.coordinates_to_bed_string(coord.get_transcript_coordinates(), m_sequence).getBytes());
            //std::cout << coordinates_to_bed_string((*it)->get_transcript_coordinates(), m_sequence) << std::endl;

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
                if (m_geneunique && !m_transcriptunique) {
                    colour = "0,0,0";
                } else if (m_geneunique && m_transcriptunique) {
                    colour = "204,0,0";
                }
            }

            os.write((colour + "\t" + exon_count + "\t" + exon_lengths.toString() + "\t" + exon_starts.toString() + "\n").getBytes());
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
        for (PeptideCoordinates coord : m_pepcoordinates) {
            count += 1;
            if (m_pepcoordinates.size() > 1) {
                sequence_add = "." + count;
            }

            os.write((geneID + "." + m_sequence + sequence_add + "\t\"" + geneID + "|@").getBytes());

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
            Pair<ArrayList<Integer>, ArrayList<Double>> pair = m_tissuetags.get(tissuelist.get(i));
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
        for (Map.Entry<String, Map<String, PTMEntry>> ptm_it : m_pepforms.entrySet()) {
            for (Map.Entry<String, PTMEntry> ptm_single_it : ptm_it.getValue().entrySet()) {
                List<Pair<PeptideCoordinates, GenomeCoordinates>> coord = ptm_single_it.getValue().get_genome_coordinates();
                for (Pair<PeptideCoordinates, GenomeCoordinates> coord_it : coord) {
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
                            + "\t" + "\t" + colour + "\t" + exon_count + "\t" + exon_lengths.toString() + "\t"
                            + exon_starts.toString() + "\n").getBytes());
                }
            }
        }
        return os;
    }

    //adds new tissuetags if they havent existed before.
    public final void add_tags(String tag, int sigPSMs, double quant) {

        Pair<ArrayList<Integer>, ArrayList<Double>> pair = m_tissuetags.computeIfAbsent(tag, k -> new Pair<>(new ArrayList<>(), new ArrayList<>()));

        pair.getKey().add(sigPSMs);
        pair.getValue().add(quant);
    }

    //adds ptms if a sequence matches another with ptms.
    public final void add_ptm(String ptmsequence) {
        //adds the modfication to the peptide entry
        m_pepforms.put(ptmsequence, ptm_set(ptmsequence));
        for (Map.Entry<String, PTMEntry> ptm_it : m_pepforms.get(ptmsequence).entrySet()) {
            for (PeptideCoordinates peptide_coordinates_it : m_pepcoordinates) {
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
    private static ArrayList<Pair<Coordinates, GenomeCoordinates>> create_coordinate_map_type(ArrayList<GenomeCoordinates> genomecoords) {
        ArrayList<Pair<Coordinates, GenomeCoordinates>> coord_map = new ArrayList<>();
        Coordinates prev_prot_coord = new Coordinates();
        prev_prot_coord.setCterm(Offset.off3);
        prev_prot_coord.setNterm(Offset.off3);
        prev_prot_coord.setStart(0);
        prev_prot_coord.setEnd(0);

        for (GenomeCoordinates genome_coordinates_itr : genomecoords) {
            Coordinates prot_coord = new Coordinates();
            GenomeCoordinates genome_coordinates = new GenomeCoordinates(genome_coordinates_itr);
            if (genome_coordinates.getStrand() == Strand.rev) {
                int start = genome_coordinates.getEnd();
                int end = genome_coordinates.getStart();
                genome_coordinates.setStart(start);
                genome_coordinates.setEnd(end);
            }

            if (prev_prot_coord.getCterm() != Offset.off3) {
                prot_coord.setNterm(Offset.forValue(3 - prev_prot_coord.getCterm().getValue()));
            } else {
                prot_coord.setNterm(Offset.off3);
            }

            int length;

            if (genome_coordinates.getStrand() == Strand.fwd) {
                length = genome_coordinates.getEnd() - genome_coordinates.getStart() + 1;
            } else {
                length = genome_coordinates.getStart() - genome_coordinates.getEnd() + 1;
            }
            // calc cterm
            if (length % 3 == 0) {
                if (prot_coord.getNterm() != Offset.off3) {
                    prot_coord.setCterm(Offset.forValue(3 - prot_coord.getNterm().getValue()));
                } else {
                    prot_coord.setCterm(Offset.off3);
                }
            } else if (length % 3 == 2) {
                if (prot_coord.getNterm() == Offset.off3) {
                    prot_coord.setCterm(Offset.off2);
                } else if (prot_coord.getNterm() == Offset.off2) {
                    prot_coord.setCterm(Offset.off3);
                } else if (prot_coord.getNterm() == Offset.off1) {
                    prot_coord.setCterm(Offset.off1);
                }
            } else if (length % 3 == 1) {
                if (prot_coord.getNterm() == Offset.off3) {
                    prot_coord.setCterm(Offset.off1);
                } else if (prot_coord.getNterm() == Offset.off1) {
                    prot_coord.setCterm(Offset.off3);
                } else if (prot_coord.getNterm() == Offset.off2) {
                    prot_coord.setCterm(Offset.off2);
                }
            }

            // calc protein coordinates
            if (prot_coord.getNterm() != Offset.off3) {
                prot_coord.setStart(prev_prot_coord.getEnd());
            } else {
                if (prev_prot_coord.getEnd() == 0 && coord_map.isEmpty()) {
                    prot_coord.setStart(0);
                } else {
                    prot_coord.setStart(prev_prot_coord.getEnd() + 1);
                }
            }

            int offsets = 0;
            if (prot_coord.getNterm() != Offset.off3) {
                offsets = offsets + prot_coord.getNterm().getValue();
            }

            if (genome_coordinates.getStrand() == Strand.fwd) {
                length = genome_coordinates.getEnd() - genome_coordinates.getStart() + 1 - offsets;
            } else {
                length = genome_coordinates.getStart() - genome_coordinates.getEnd() + 1 - offsets;
            }

            int pep_length = length / 3;

            int pep_end = prot_coord.getStart() + pep_length - 1;
            if (prot_coord.getCterm() != Offset.off3) {
                pep_end = pep_end + 1;
            }
            if (prot_coord.getNterm() != Offset.off3) {
                pep_end = pep_end + 1;
            }

            prot_coord.setEnd(pep_end);

            prev_prot_coord = prot_coord;

            coord_map.add(new Pair<>(prot_coord, genome_coordinates));
        }
        return coord_map;
    }

    //adds a peptide. this funciton is used if this specific peptide has not yet been found. it will also find the peptides genomic coordinates.
    public final void add_peptide(CoordinateWrapper coordwrapper, String sequence, String ptmSequence, String tag, int sigPSMs, TranscriptsT transcripts, int genes, FileOutputStream ofstream, double quant) {
        if (m_sequence.isEmpty()) {
            m_sequence = sequence;
        }
        if (m_numtranscripts == 0) {
            m_numtranscripts = transcripts.getM_entries().size();
        }

        if (genes > 1) {
            m_geneunique = false;
        }
        if (transcripts.getM_entries().size() > 1) {
            m_transcriptunique = false;
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
                m_transcriptids.addAll(transcriptids);
                TreeSet<String> exonids = pep_coord.get_exon_ids();
                m_exonids.addAll(exonids);
                if (genomic_coordinate.size() != 0) {
                    //and saves them.
                    m_pepcoordinates.add(pep_coord);
                    //sets start and end coord of the peptide entry to min/max values. these are used for comparing PeptideEntry objects.
                    if (m_startcoord == 0) {
                        m_startcoord = pep_coord.get_transcript_coordinates().getStart();
                    }
                    if (m_endcoord == 0) {
                        m_endcoord = pep_coord.get_transcript_coordinates().getEnd();
                    }
                    if (pep_coord.get_transcript_coordinates().getStart() < m_startcoord) {
                        m_startcoord = pep_coord.get_transcript_coordinates().getStart();
                    }
                    if (pep_coord.get_transcript_coordinates().getEnd() > m_endcoord) {
                        m_endcoord = pep_coord.get_transcript_coordinates().getEnd();
                    }
                }
            }
        }
        add_ptm(ptmSequence);
    }

    //adds a peptide. this function is used if a peptide has appeared before.
    public final void add_peptide(String ptmsequence, String tag, int sigPSMs, double quant) {
        if (!m_sequence.equals(ptmsequence)) {
            add_ptm(ptmsequence);
        }
        add_tags(tag, sigPSMs, quant);
    }

    private String transcriptids_to_string() {
        StringBuilder transcriptid_string = new StringBuilder();
        for (String m_transcriptid : m_transcriptids) {
            if (!transcriptid_string.toString().isEmpty()) {
                transcriptid_string.append("|");
            }
            transcriptid_string.append(m_transcriptid);
        }
        return transcriptid_string.toString();
    }

    private String exonids_to_string() {
        StringBuilder exonid_string = new StringBuilder();
        for (String it : m_exonids) {
            if (!exonid_string.toString().isEmpty()) {
                exonid_string.append("|");
            }
            exonid_string.append(it);
        }
        return exonid_string.toString();
    }

    //checks whether the peptide entry has a modified form. Returns true if m_pepforms only contains peptide without PTM (size==1)
    public final boolean noPTM() {
        return m_pepforms.containsKey(m_sequence);
    }
}