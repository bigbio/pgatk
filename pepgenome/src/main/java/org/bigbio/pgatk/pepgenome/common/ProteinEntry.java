package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.PepGenomeTool;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ProteinEntry implements Serializable {
    private static final long serialVersionUID = -1732196455282495216L;
    //whole fasta header
    private String m_fasta_header;
    //transcript id
    private String m_transcript_id;
    //gene id
    private String m_gene_id;
    // Translation offset
    private int m_translation_offset = 0;
    //the AA sequence.
    private String m_aa_sequence;

    // Default regex patterns.
    private static final Pattern originalGENEPATTERN = Pattern.compile("gene:([^\\s\\.]*)[^\\s]*\\s");
    private static final Pattern originalTRANSCRIPTPATTERN = Pattern.compile("transcript:([^\\s\\.]*)[^\\s]*\\s");

    // Exco regex patterns
    private static final Pattern excoGENEPATTERN = Pattern.compile("geneID=([^\\s]*)[^\\s]*\\s");
    private static final Pattern excoTRANSCRIPTPATTERN = Pattern.compile("transcriptID=([^\\s]*)[^\\s]*");
    private static final Pattern excoOFFSETPATTERN = Pattern.compile(" offset=(\\d*) ");

    //the first Coordinates are the coordinates of exons within the protein and the GenomeCoordinate is its corresponding location in the genome.
    private ArrayList<Tuple<Coordinates, GenomeCoordinates>> m_coordinates_map;
    //check if the coding sequence is dividable by 3bp and not offset due to incomplete transcript annotation.
    private int m_cds_annotation_correct;

    public ProteinEntry() {
        this.m_fasta_header = "";
        this.m_transcript_id = "";
        this.m_gene_id = "";
        this.m_aa_sequence = "";
        this.m_coordinates_map = new ArrayList<>();
        this.m_cds_annotation_correct = 0;
        this.m_translation_offset = 0;
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
            m_aa_sequence = AAsequence;
            m_coordinates_map = new ArrayList<>();
            m_cds_annotation_correct = 0;

            m_gene_id = extract_gene_id_fasta(fastaHeader);
            m_transcript_id = extract_transcript_id_fasta(fastaHeader);

            if (PepGenomeTool.useExonCoords) {
                // Exco mode
                // Using exon coords in place of CDS, offset describes distance in nucleotides from exon start to translation start.
                m_translation_offset = extract_offset_fasta(fastaHeader);

            }
            else {
                // Not using exon coords in place of CDS - Using original format
                m_translation_offset = 0;
            }

            // Put transcript ID and offset value into translation offset map in PepGenomeTool
            PepGenomeTool.m_translation_offset_map.put(m_transcript_id, m_translation_offset);

        }
    }

    private String extract_gene_id_fasta(String str) {
        String value = "";

        Matcher geneMatcher;
        if (PepGenomeTool.useExonCoords) {
            geneMatcher= excoGENEPATTERN.matcher(str);
        } else {
            geneMatcher= originalGENEPATTERN.matcher(str);
        }

        if (geneMatcher.find()) {
            value = geneMatcher.group(1);
            //System.out.println("Gene extracted correctly: "+value);
        }
        else {

            // From original method
            String[] split = str.split("\\|");
            if(split.length==8) {
                String[] dotsplit = split[2].split("\\.");
                value = dotsplit[0];
            }
        }

        return value;
    }

    // Extracts the transcript Id from a fasta header
    private String extract_transcript_id_fasta(String str) {
        String value = "";

        Matcher transcriptMatcher;
        if (PepGenomeTool.useExonCoords) {
            transcriptMatcher = excoTRANSCRIPTPATTERN.matcher(str);
        } else {
            transcriptMatcher = originalTRANSCRIPTPATTERN.matcher(str);

            // From original method
            String[] split = str.split("\\|");
            if(split.length==8) {
                String[] dotsplit = split[1].split("\\.");
                value = dotsplit[0];
            }
        }

        if (transcriptMatcher.find()) {
            value = transcriptMatcher.group(1);
            //System.out.println("Transcript extracted correctly"+value);
        } else {
            //System.out.println("Transcript not extracted correctly: "+value);
        }
        return value;
    }

    // Extracts the translation start point offset from a fasta header
    private Integer extract_offset_fasta(String str) {
        String value = "0"; // Default if no offset given.
        Matcher offsetMatcher = excoOFFSETPATTERN.matcher(str);
        if (offsetMatcher.find()) {
            value = offsetMatcher.group(1);
            //System.out.println("Offset extracted correctly");
            //System.out.println("Offset value: "+value);
        }

        m_translation_offset = Integer.parseInt(value);

        return m_translation_offset;
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
    public void set_coordinate_map(ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinatesMap) {
        m_coordinates_map = coordinatesMap;
    }

    //getter for CDS_annotation_correct.
    public int get_cds_annotation_correct() {
        return m_cds_annotation_correct;
    }

    //returns the genomic coordinates.
    //mapping function. takes the positions calculated in the KmereMap and generates the genomic coordinates for all peptides of this protein.
    public ArrayList<ArrayList<GenomeCoordinates>> find_coordinates(int peptideseqSize, ArrayList<PositionMismatchT> positions) {
        ArrayList<ArrayList<GenomeCoordinates>> foundCoordinates = new ArrayList<>();
        Coordinates peptideCoordinates = new Coordinates();
        peptideCoordinates.setCterm(Offset.off3);
        peptideCoordinates.setNterm(Offset.off3);

        //iterate all found positions
        for (PositionMismatchT current : positions) {
            peptideCoordinates.setStart(current.position_in_protein());
            peptideCoordinates.setEnd(current.position_in_protein() + (peptideseqSize - 1));
            ArrayList<GenomeCoordinates> single = new ArrayList<>();
            m_coordinates_map.stream().filter(e -> e.getKey().equals(peptideCoordinates))
                    .forEach(fe -> {
                        Tuple<Coordinates, GenomeCoordinates> coordinatesPartial = Utils.get_coordinates(fe.getKey(), fe.getValue(), peptideCoordinates);
                        single.add(coordinatesPartial.getValue());
                    });
            //and has to be done several times to find all peptides.
            foundCoordinates.add(single);
        }
        return foundCoordinates;
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