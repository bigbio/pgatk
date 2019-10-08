package org.bigbio.pgatk.pepgenome.common;



import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

//in this class all the information will be used to find the
//genomic location of a peptide.
public class PeptideCoordinates implements Comparable<PeptideCoordinates> {

    //holds the found coordinates.
    private ArrayList<GenomeCoordinates> m_coordinate_list;
    //passed to the constuctor. holds the peptides coordinates and the associated genomic coordinates.
    //peptide coordinates sorted; used for computing ptm genomic coordinates.
    private ArrayList<Tuple<Coordinates, GenomeCoordinates>> m_coordinates;
    //the coordinates of the current transcript.
    private GenomeCoordinates m_transcript_coordinates;
    //tests if the annotation of the coding sequence is dividable by 3bp and not offset
    //due to incomplete transcript annotation.
    private int m_cds_annotation_correct;

    private TreeSet<String> m_transcriptids;
    private TreeSet<String> m_exonids;


    public PeptideCoordinates() {
        this.m_coordinate_list = new ArrayList<>();
        this.m_coordinates = new ArrayList<>();
        this.m_transcript_coordinates = new GenomeCoordinates();
        this.m_cds_annotation_correct = 0;
        this.m_transcriptids = new TreeSet<>();
        this.m_exonids = new TreeSet<>();
    }

    public PeptideCoordinates(ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinates, int CDSannotationcorrect) {
        this.m_coordinate_list = new ArrayList<>();
        this.m_coordinates = coordinates;
        this.m_cds_annotation_correct = CDSannotationcorrect;
        get_exon_coordinates();
        m_transcript_coordinates = get_transcript_coordinates();
        m_transcriptids = new TreeSet<>();
        m_exonids = new TreeSet<>();
        add_ids();
    }

    //this function is the actual search function that will find genomic coordinates for a peptide.
    public ArrayList<GenomeCoordinates> find_coordinates(Tuple<Integer, Integer> peptideproteincoords) {
        Coordinates ptm_coordinates = new Coordinates();
        ptm_coordinates.start = peptideproteincoords.getKey();
        ptm_coordinates.end = peptideproteincoords.getValue();
        ptm_coordinates.Cterm = Offset.off3;
        ptm_coordinates.Nterm = Offset.off3;

        ArrayList<GenomeCoordinates> ptm = new ArrayList<>();

        m_coordinates.stream().filter(e -> e.getKey().equals(ptm_coordinates))
                .forEach(fe -> {
                    Tuple<Coordinates, GenomeCoordinates> coordinates_partial = Utils.get_coordinates(fe.getKey(), fe.getValue(), ptm_coordinates);
                    ptm.add(coordinates_partial.getValue());
                });

        //sorted
        ptm.sort(getGenomeCoordinatesComparator());
        return ptm;
    }

    //generates exon coordinates from the coordinate map type.
    public ArrayList<GenomeCoordinates> get_exon_coordinates() {
        //used to generate the exon information in the output
        //this is mainly important to find reverse exons.
        if (m_coordinate_list.isEmpty()) {
            for (Tuple<Coordinates, GenomeCoordinates> gcs : m_coordinates) {
                GenomeCoordinates val = new GenomeCoordinates(gcs.getValue());
                if (val.getStrand() == Strand.rev) {
                    int start = val.end;
                    int end = val.start;
                    val.start = start;
                    val.end = end;
                }
                m_coordinate_list.add(val);
            }
            m_coordinate_list.sort(getGenomeCoordinatesComparator());
        }
        return m_coordinate_list;
    }

    private Comparator<GenomeCoordinates> getGenomeCoordinatesComparator() {
        return (lhs, rhs) -> {
            if (Utils.compare_coordinates_ascending(lhs, rhs)) {
                return -1;
            }
            if (Utils.compare_coordinates_ascending(rhs, lhs)) {
                return 1;
            }
            return 0;
        };
    }

    //generates common.GenomeCoordinates for the associated transcript.
    public GenomeCoordinates get_transcript_coordinates() {
        GenomeCoordinates transcript_coordinates = new GenomeCoordinates();
        for (int i = 0; i < m_coordinate_list.size(); ++i) {
            GenomeCoordinates transcript_coordinates1 = new GenomeCoordinates(m_coordinate_list.get(i));
            if (i == 0) {
                transcript_coordinates = transcript_coordinates1;
            } else {
                if (transcript_coordinates.start > transcript_coordinates1.start) {
                    transcript_coordinates.start = transcript_coordinates1.start;
                }
                if (transcript_coordinates.end < transcript_coordinates1.end) {
                    transcript_coordinates.end = transcript_coordinates1.end;
                }
            }
        }

        if (m_cds_annotation_correct != 0) {
            if (transcript_coordinates.getStrand() == Strand.fwd) {
                transcript_coordinates.start = transcript_coordinates.start - m_cds_annotation_correct;
                transcript_coordinates.end = transcript_coordinates.end - m_cds_annotation_correct;
            } else if (transcript_coordinates.getStrand() == Strand.rev) {
                transcript_coordinates.end = transcript_coordinates.end + m_cds_annotation_correct;
                transcript_coordinates.start = transcript_coordinates.start + m_cds_annotation_correct;
            }
        }
        return transcript_coordinates;
    }

    public final TreeSet<String> get_trasncript_ids() {
        return m_transcriptids;
    }

    public final TreeSet<String> get_exon_ids() {
        return m_exonids;
    }

    private void add_ids() {
        for (Tuple<Coordinates, GenomeCoordinates> gcs : m_coordinates) {
            GenomeCoordinates val = gcs.getValue();
            m_transcriptids.add(val.getTranscriptid());
            m_exonids.add(val.getExonid());
        }
    }

    //lesser than operator. returns compare_coordinates_ascending_whole (this, other)
    //otherwise returns false if compare_coordinates_ascending_whole (other, this)
    //otherwise returns compare_genome_coordinate_sets_ascending(this, other)
    boolean lessThan(PeptideCoordinates rhs) {
        if (Utils.compare_coordinates_ascending_whole(m_transcript_coordinates, rhs.m_transcript_coordinates)) {
            return true;
        } else if (Utils.compare_coordinates_ascending_whole(rhs.m_transcript_coordinates, m_transcript_coordinates)) {
            return false;
        }
        return Utils.compare_genome_coordinate_sets_ascending(m_coordinate_list, rhs.m_coordinate_list);
    }

    @Override
    public int compareTo(PeptideCoordinates o) {
        if (lessThan(o)) {
            return -1;
        }
        if (o.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    @Override
    public boolean equals(Object obj) {
        return compareTo((PeptideCoordinates) obj) == 0;
    }
}