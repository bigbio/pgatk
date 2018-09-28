package bigbio.pgatk.jpogo;

import bigbio.pgatk.jpogo.common.Chromosome;
import bigbio.pgatk.jpogo.common.GenomeCoordinates;
import bigbio.pgatk.jpogo.common.GenomeMapper;
import bigbio.pgatk.jpogo.common.Utils;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GeneEntry implements Comparable<GeneEntry> {

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
    public int compareTo(GeneEntry otherInstance) {
        if (lessThan(otherInstance)) {
            return -1;
        } else if (otherInstance.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    public GeneEntry() {
    }

    public GeneEntry(String gtfgeneline) {
        init(gtfgeneline);
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

    //comapres two genes.
    //returns true if the chromosome number is smaller than rhs' chromosome number
    //otherwise returns true if the startposition in the chromosome is smaller
    //otherwise returns true if the endposition in the chromosome is smaller
    //otherwise returns false.
    public boolean lessThan(GeneEntry rhs) {
        if (String.valueOf(m_coord.getChr().getValue()).compareTo(String.valueOf(rhs.m_coord.getChr().getValue())) < 0) {
            return true;
        }
        if (m_coord.getChr() == rhs.m_coord.getChr()) {
            if (m_coord.getStart() == rhs.m_coord.getStart()) {
                return m_coord.getEnd() < rhs.m_coord.getEnd();
            }
            return m_coord.getStart() < rhs.m_coord.getStart();
        }
        return false;
    }

    //converts a gene into a gtf line and prints it to the given outputstream.
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

        for (String m_tag : m_tags) {
            os.write((" tag \"" + m_tag + "\";").getBytes());
        }
        return os;
    }

    //check if gene entry maps to chromosomes (e.g. chr1, chrX, 10, etc.)
    public final boolean is_primary() {
        return m_coord.getChr() != Chromosome.chrNA && m_coord.getChr() != Chromosome.scaffold && m_coord.getChr() != Chromosome.chrXY;
    }

    //check if gene entry maps to patch, haplotype or scaffold
    public final boolean is_patchhaploscaff() {
        return m_coord.getChr() == Chromosome.scaffold && !m_coord.getChrscaf().equals("");
    }

    //looks for the text specified in common.GenomeMapper::ID::geneId and returns the ID (including the trailing number of length common.GenomeMapper::ID::length - common.GenomeMapper::ID::geneId.length()).
    public static String extract_gene_id(String gtfGeneLine) {
        int start = gtfGeneLine.indexOf(GenomeMapper.ID.GENE_ID);
        String value = "";
        if (start != -1) {
            if ((start + GenomeMapper.ID.LENGTH) < gtfGeneLine.length()) {
                value = gtfGeneLine.substring(start, start + GenomeMapper.ID.LENGTH);
            }
        }
        return value;
    }

    //looks for the text specified in common.GenomeMapper::ID::transcriptId and returns the ID (including the trailing number of length common.GenomeMapper::ID::length - common.GenomeMapper::ID::transcriptId.length()).
    public static String extract_transcript_id(String gtfGeneLine) {
        int index = gtfGeneLine.indexOf(GenomeMapper.ID.TRANSCRIPT_ID);
        String value = "";
        if (index != -1) {
            if ((index + GenomeMapper.ID.LENGTH) < gtfGeneLine.length()) {
                value = gtfGeneLine.substring(index, index + GenomeMapper.ID.LENGTH);
            }
        }
        return value;
    }

    //looks for the text specified in GenomeMapper::ID::exonId and returns the ID (including the trailing number of length GenomeMapper::ID::length - GenomeMapper::ID::exonId.length()).
    public static String extract_exon_id(String gtfGeneLine) {
        int index = gtfGeneLine.indexOf(GenomeMapper.ID.EXON_ID);
        String value = "";
        if (index != -1) {
            if ((index + GenomeMapper.ID.LENGTH) < gtfGeneLine.length()) {
                value = gtfGeneLine.substring(index, index + GenomeMapper.ID.LENGTH);
            }
        }
        return value;
    }

    private void init(String gtfgeneline) {
        ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(gtfgeneline, "\t")));
        init(extract_gene_id(gtfgeneline), Utils.extract_coordinates_from_gtf_line(tokens), extract_type(tokens), extract_status(tokens), extract_gene_name(tokens), extract_tags(tokens));
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

    //after tokenizing the geneline these functions can be used to extract information.
    //extracts the gene type (proein_coding,...)
    private static String extract_type(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = extract_by_tag("gene_type", tokens.get(8));
            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }

    //extracts the gene status (KNOWN,...)
    private static String extract_status(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = extract_by_tag("gene_status", tokens.get(8));
            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }

    //extracts the gene symbol
    private static String extract_gene_name(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = extract_by_tag("gene_name", tokens.get(8));
            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
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
    private static List<String> extract_by_tag(String tag, String tagList) {
        List<String> return_values = new ArrayList<>();
        String[] values = Utils.tokenize(tagList, ";", true);
        for (String value : values) {
            int start = 0;
            if (Utils.getCppStyleSubString(value, 0, 1).equals(" ")) {
                start = 1;
            }
            if (Utils.getCppStyleSubString(value, start, tag.length()).equals(tag)) {
                String[] v = Utils.tokenize(value, "\"");
                if (v.length >= 2) {
                    return_values.add(v[1]);
                }
            }
        }
        return return_values;
    }
}