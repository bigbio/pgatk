package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.common.*;
import javafx.util.Pair;

import java.io.Serializable;
import java.util.ArrayList;

public class ProteinEntry implements Serializable {
    private static final long serialVersionUID = -1732196455282495216L;
    //whole fasta header
    private String m_fasta_header;
    //transcript number
    private String m_transcript_id;
    //gene number
    private String m_gene_id;
    //the AA sequence.
    private String m_aa_sequence;

    //std::multimap <Coordinates (protein coordinates), GenomeCoordinates(corresponding genomic coordinates), Coordinates (passing this as third argument will use the Coordinates::operator() as comparator)>
    //the first Coordinate are the coordinates of exons within the protein and the GenomeCoordinate is its corresponding location in the genome.
    private ArrayList<Pair<Coordinates, GenomeCoordinates>> m_coordinates_map;
    //check if the coding sequence is dividable by 3bp and not offset due to incomplete transcript annotation.
    private int m_cds_annotation_correct;

    public ProteinEntry() {
        this.m_fasta_header = "";
        this.m_transcript_id = "";
        this.m_gene_id = "";
        this.m_aa_sequence = "";
        this.m_coordinates_map = new ArrayList<>();
        this.m_cds_annotation_correct = 0;
    }

    public ProteinEntry(String fastaHeader, String AAsequence) {
        init(fastaHeader, AAsequence);
    }

    public ProteinEntry(FastaEntry fastaEntry) {
        init(fastaEntry);
    }

    //QOL: delegates to void init(std::string fastaHeader, std::string AAsequence);
    private void init(FastaEntry fastaEntry) {
        init(fastaEntry.get_header(), fastaEntry.get_sequence());
    }

    //sets all crucial values.
    private void init(String fastaHeader, String AAsequence) {
        if (fastaHeader.substring(0, 1).equals(">")) {
            m_fasta_header = fastaHeader;
            m_transcript_id = extract_transcript_id_fasta(fastaHeader);
            m_gene_id = extract_gene_id_fasta(fastaHeader);
            m_aa_sequence = AAsequence;
            m_coordinates_map = new ArrayList<>();
            m_cds_annotation_correct = 0;
        }
    }

    //gets the transcriptId from a fasta header
    private String extract_transcript_id_fasta(String str) {
        int index = str.indexOf(GenomeMapper.ID.TRANSCRIPT_ID);
        String value = "";

        if (index != -1) {
            if ((index + (GenomeMapper.ID.LENGTH - 1)) < str.length()) {
                value = str.substring(index, index + GenomeMapper.ID.LENGTH);
            }
        }
        return value;
    }

    //gets the gene id from a fasta header
    private String extract_gene_id_fasta(String str) {
        int start = str.indexOf(GenomeMapper.ID.GENE_ID);
        String value = "";
        if (start != -1) {
            if ((start + GenomeMapper.ID.LENGTH - 1) < str.length()) {
                value = str.substring(start, start + GenomeMapper.ID.LENGTH);
            }
        }
        return value;
    }

    //returns the transcript_id number of the current protein
    public String get_transcript_id() {
        return m_transcript_id;
    }

    //returns the gene_id number of the current protein
    public String get_gene_id() {
        return m_gene_id;
    }

    //returns the sequence. the sequences are iso-sequences (I and L are converted to J)
    public String get_sequence() {
        return m_aa_sequence;
    }

    //setter for the coordinatesMap
    public void set_coordinate_map(ArrayList<Pair<Coordinates, GenomeCoordinates>> coordinatesMap) {
        m_coordinates_map = coordinatesMap;
    }

    //getter for CDS_annotation_correct.
    public int get_cds_annotation_correct() {
        return m_cds_annotation_correct;
    }

    //returns the genomic coordinates.
    //mapping function. takes the positions calculated in the KmereMap and generates the genomic coordinates for all peptides of this protein.
    public ArrayList<ArrayList<GenomeCoordinates>> find_coordinates(int peptideseqSize, ArrayList<PositionMismatchT> positions) {
        ArrayList<ArrayList<GenomeCoordinates>> found_coordinates = new ArrayList<>();
        Coordinates peptide_coordinates = new Coordinates();
        peptide_coordinates.setCterm(Offset.off3);
        peptide_coordinates.setNterm(Offset.off3);

        //iterate all found positions
        for (PositionMismatchT current : positions) {
            peptide_coordinates.setStart(current.position_in_protein());
            peptide_coordinates.setEnd(current.position_in_protein() + (peptideseqSize - 1));
            ArrayList<GenomeCoordinates> single = new ArrayList<>();
            m_coordinates_map.stream().filter(e -> e.getKey().equals(peptide_coordinates))
                    .forEach(fe -> {
                        Pair<Coordinates, GenomeCoordinates> coordinatesPartial = Utils.get_coordinates(fe.getKey(), fe.getValue(), peptide_coordinates);
                        single.add(coordinatesPartial.getValue());
                    });
            //and has to be done several times to find all peptides.
            found_coordinates.add(single);
        }
        return found_coordinates;
    }

    @Override
    public String toString() {
        return "ProteinEntry{" +
                "m_fasta_header='" + m_fasta_header + '\'' +
                ", m_transcript_id='" + m_transcript_id + '\'' +
                ", m_gene_id='" + m_gene_id + '\'' +
                ", m_aa_sequence='" + m_aa_sequence + '\'' +
                ", m_coordinates_map=" + m_coordinates_map +
                ", m_cds_annotation_correct=" + m_cds_annotation_correct +
                '}';
    }
}