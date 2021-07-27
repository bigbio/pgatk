package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.io.GFF3Parser;
import org.bigbio.pgatk.pepgenome.io.GTFParser;

import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GeneEntry implements Comparable<GeneEntry>, Serializable {

    private static final long serialVersionUID = 1099018767895133035L;
    //the genomic coordinates for a gene are stored here.
    private GenomeCoordinates m_coord = new GenomeCoordinates();
    //gene id as extracted by extract_gene_id();
    private String m_id;
    //type (protein_coding,...)
    private String m_type;
    //status (KNOWN,...)
    private String m_status;
    //gene name (gene symbol)
    private String m_gene_name;
    //tags. (ncRNA_host,...)
    private List<String> m_tags = new ArrayList<>();

    @Override
    public int compareTo(GeneEntry o) {
        if (isLessThan(o)) {
            return -1;
        } else if (o.isLessThan(this)) {
            return 1;
        }
        return 0;
    }

    public GeneEntry() {
    }

    public GeneEntry(String annotationGeneLine) {
        init(annotationGeneLine);
    }

    //returns the gene ID
    public final String get_id() {
        return m_id;
    }

    //returns the type
    public final String get_type() {
        return m_type;
    }

    //returns the genes status
    public final String get_status() {
        return m_status;
    }

    //returns the gene name
    public final String get_name() {
        return m_gene_name;
    }

    /**
     * Compares two genes and returns true if the chromosome number is smaller than rhs' chromosome number
     * otherwise returns true if the start position in the chromosome is smaller otherwise returns
     * true if the end position in the chromosome is smaller otherwise returns false.
     *
     * @param {org.bigbio.pgatk.pepgenome.common.GeneEntry }
     * @return true is lessThan
     */
    public boolean isLessThan(GeneEntry rhs) {
        if (String.valueOf(m_coord.getChr().getValue()).compareTo(String.valueOf(rhs.m_coord.getChr().getValue())) < 0) {
            return true;
        }
        if (m_coord.getChr().getValue() == rhs.m_coord.getChr().getValue()) {
            if (m_coord.getStart() == rhs.m_coord.getStart()) {
                return m_coord.getEnd() < rhs.m_coord.getEnd();
            }
            return m_coord.getStart() < rhs.m_coord.getStart();
        }
        return false;
    }

    //converts a gene into a gtf line and prints it to the given output stream.
    public final OutputStream to_gtf(String source, OutputStream os) throws Exception {
        return to_gtf(source, os, true);
    }

    public final OutputStream to_gtf(String source) throws Exception {
        return to_gtf(source, System.out, true);
    }

    public final OutputStream to_gtf(String source, OutputStream os, boolean chrincluded) throws Exception {
        os.write(Utils.coordinates_to_gtf_string(m_coord, "gene", false, source, chrincluded).getBytes());
        os.write(("gene_id \"" + m_id + "\"; transcript_id \"" + m_id + "\"; gene_type \"" + m_type + "\"; gene_status \""
                + m_status + "\"; gene_name \"" + m_gene_name).getBytes());
        os.write(("\"; transcript_type \"" + m_type + "\"; transcript_status \"" + m_status + "\"; transcript_name \""
                + m_gene_name + "\";").getBytes());

        for (String mTag : m_tags) {
            os.write((" tag \"" + mTag + "\";").getBytes());
        }
        return os;
    }

    //check if gene entry maps to chromosomes (e.g. chr1, chrX, 10, etc.)
    public final boolean is_primary() {
        return !m_coord.getChr().isNA() && !m_coord.getChr().isScaffold();
    }

    //check if gene entry maps to patch, haplotype or scaffold
    public final boolean is_patchhaploscaff() {
        return m_coord.getChr().isScaffold() && !m_coord.getChrscaf().equals("");
    }

    private void init(String annotationGeneLine) {
        ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(annotationGeneLine, "\t")));

        if (GTFParser.instance != null) {
            init(GTFParser.extract_gene_id(annotationGeneLine), Utils.extract_coordinates_from_gtf_line(tokens), GTFParser.extract_type(tokens), GTFParser.extract_status(tokens), GTFParser.extract_gene_name(tokens), extract_tags(tokens));
        }

        else if (GFF3Parser.instance != null) {
            init(GFF3Parser.extract_gene_id(annotationGeneLine), Utils.extract_coordinates_from_gtf_line(tokens), "", "", GFF3Parser.extract_gene_name(tokens), extract_tags(tokens));
        }

    }

    //private member function called in constructor to initialize a gene entry
    private void init(String ID, GenomeCoordinates coordinates, String type, String status, String gene_name) {
        init(ID, coordinates, type, status, gene_name, new ArrayList<>());
    }

    private void init(String ID, GenomeCoordinates coordinates, String type, String status, String gene_name, List<String> tags) {
        m_coord = coordinates;
        m_id = ID;
        m_type = type;
        m_status = status;
        m_gene_name = gene_name;
        m_tags = tags;
    }

    //extracts all possible tags to  be used with extract_by_tag.
    private static List<String> extract_tags(List<String> tokens) {
        List<String> values = new ArrayList<>();
        if (tokens.size() >= 9) {
            values = extract_by_tag("tag", tokens.get(8));
        }
        return values;
    }

    //tag takes the name of the tag to extract, tagList contains all tags separated by \.
    //possible values for tag could be gene_name, gene_type,...
    public static List<String> extract_by_tag(String tag, String tagList) {
        List<String> rValues = new ArrayList<>();
        String[] values = Utils.tokenize(tagList, ";", true);
        for (String value : values) {
            int start = 0;
            if (Utils.getCppStyleSubString(value, 0, 1).equals(" ")) {
                start = 1;
            }
            if (Utils.getCppStyleSubString(value, start, tag.length()).equals(tag)) {
                String[] v = Utils.tokenize(value, "\"");
                if (v.length >= 2) {
                    rValues.add(v[1]);
                }
            }
        }
        return rValues;
    }
}
