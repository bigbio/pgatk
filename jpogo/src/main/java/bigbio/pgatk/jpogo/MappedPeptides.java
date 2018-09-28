package bigbio.pgatk.jpogo;

import bigbio.pgatk.jpogo.common.*;
import javafx.util.Pair;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.*;

//the mapped peptides class connects the different classes
//and offers functions to wite all found peptides.
public class MappedPeptides {

    //maps the peptide sequence (without PTM) to the corresponding MapEntry (chromosomes).
    private Map<String, MapEntry> m_mapping = new TreeMap<>();
    //maps the peptide sequence (without PTM) to the corresponding MapEntry (patches, haplotypes, scaffolds).
    private Map<String, MapEntry> m_mapping_phs = new TreeMap<>();
    //counts the number of found peptides. similar to m_mapping.size().
    private int m_count_peptides;
    private int m_count_peptides_phs;

    //maps a tissue to an index.
    private Map<String, Integer> m_tissuemap = new TreeMap<>();
    //the index, to be incremented durig the program.
    private int m_tissueindex;

    public MappedPeptides() {
//        this.m_mapping = new TreeMap<>();
//        this.m_mapping_phs = new TreeMap<>();
        this.m_count_peptides = 0;
//        this.m_tissuemap = new TreeMap<>();
        this.m_tissueindex = 0;
    }

    //adds a new gene from a gtfline.
    public final Assembly add_gene_from_gtf(String gtfGeneLine) {
        GeneEntry gene = new GeneEntry(gtfGeneLine);
        if (gene.is_primary()) {
            m_mapping.put(gene.get_id(), new MapEntry(gene));
            return Assembly.primary;
        } else if (gene.is_patchhaploscaff()) {
            m_mapping_phs.put(gene.get_id(), new MapEntry(gene));
            return Assembly.patchhaploscaff;
        }
        return Assembly.none;
    }

    //maps a transcript id to a gene id
    public final void add_transcript_id_to_gene(String gtftranscriptline) {
        String transcript_id = GeneEntry.extract_transcript_id(gtftranscriptline);
        String gene_id = GeneEntry.extract_gene_id(gtftranscriptline);
        if (m_mapping.containsKey(gene_id) && m_mapping_phs.containsKey(gene_id)) {
            m_mapping.get(gene_id).add_transcript_id(transcript_id);
            m_mapping_phs.get(gene_id).add_transcript_id(transcript_id);
        } else if (m_mapping.containsKey(gene_id) && !(m_mapping_phs.containsKey(gene_id))) {
            m_mapping.get(gene_id).add_transcript_id(transcript_id);
        } else if (!(m_mapping.containsKey(gene_id)) && m_mapping_phs.containsKey(gene_id)) {
            m_mapping_phs.get(gene_id).add_transcript_id(transcript_id);
        }
    }

    //print functions
    //converts all peptides to gtf lines
    public final void to_gtf(String filename, String source, Assembly assem) throws Exception {
        to_gtf(filename, source, assem, true);
    }

    public final void to_gtf(String filename, String source) throws Exception {
        to_gtf(filename, source, Assembly.primary, true);
    }

    public final void to_gtf(String filename, String source, Assembly assem, boolean chrincluded) throws Exception {
        FileOutputStream ofs = new FileOutputStream(filename);
        to_gtf(assem, source, ofs);
        ofs.close();
    }

    public final void to_gtf(Assembly assem, String source, OutputStream os) throws Exception {
        to_gtf(assem, source, os, true);
    }

    public final void to_gtf(Assembly assem, String source) throws Exception {
        to_gtf(assem, source, System.out, true);
    }

    public final void to_gtf(Assembly assem, String source, OutputStream os, boolean chrincluded) throws Exception {
        TreeSet<MapEntry> mapping_set = new TreeSet<>(new MapentryPCompare());

        if (assem == Assembly.primary) {
            for (Map.Entry<String, MapEntry> it : m_mapping.entrySet()) {
                mapping_set.add(it.getValue());
            }
        } else if (assem == Assembly.patchhaploscaff) {
            for (Map.Entry<String, MapEntry> it : m_mapping_phs.entrySet()) {
                mapping_set.add(it.getValue());
            }
        }
        for (MapEntry sit : mapping_set) {
            sit.to_gtf(source, os, chrincluded);
        }
    }

    //converts all peptides to bed lines
    public final void to_bed(String filename, Assembly assem) throws Exception {
        to_bed(filename, assem, true);
    }

    public final void to_bed(String filename) throws Exception {
        to_bed(filename, Assembly.primary, true);
    }

    public final void to_bed(String filename, Assembly assem, boolean chrincluded) throws Exception {
        FileOutputStream ofs = new FileOutputStream(filename);
        to_bed(assem, ofs, chrincluded);
        ofs.close();
    }

    public final void to_bed(Assembly assem, OutputStream os) throws Exception {
        to_bed(assem, os, true);
    }

    public final void to_bed(Assembly assem) throws Exception {
        to_bed(assem, System.out, true);
    }

    public final void to_bed(Assembly assem, OutputStream os, boolean chrincluded) throws Exception {
        TreeSet<MapEntry> mapping_set = new TreeSet<>(new MapentryPCompare());

        if (assem == Assembly.primary) {
            for (Map.Entry<String, MapEntry> it : m_mapping.entrySet()) {
                mapping_set.add(it.getValue());
            }
        } else if (assem == Assembly.patchhaploscaff) {
            for (Map.Entry<String, MapEntry> it : m_mapping_phs.entrySet()) {
                mapping_set.add(it.getValue());
            }
        }
        for (MapEntry sit : mapping_set) {
            sit.to_bed(os, chrincluded);
        }
    }

    //converts all peptides to gct lines
    public final void to_gct(String filename, Assembly assem) throws Exception {
        to_gct(filename, assem, true);
    }

    public final void to_gct(String filename) throws Exception {
        to_gct(filename, Assembly.primary, true);
    }

    public final void to_gct(String filename, Assembly assem, boolean chrincluded) throws Exception {
        FileOutputStream ofs = new FileOutputStream(filename);
        to_gct(assem, ofs, chrincluded);
        ofs.close();
    }

    public final void to_gct(Assembly assem, OutputStream os) throws Exception {
        to_gct(assem, os, true);
    }

    public final void to_gct(Assembly assem) throws Exception {
        to_gct(assem, System.out, true);
    }

    public final void to_gct(Assembly assem, OutputStream os, boolean chrincluded) throws Exception {
        TreeSet<MapEntry> mapping_set = new TreeSet<>(new MapentryPCompare());

        os.write("#1.2\t".getBytes());
        for (int i = 0; i < m_tissuemap.size(); ++i) {
            os.write("\t".getBytes());
        }
        if (assem == Assembly.primary) {
            os.write(("\n" + m_count_peptides).getBytes());
        } else if (assem == Assembly.patchhaploscaff) {
            os.write(("\n" + m_count_peptides_phs).getBytes());
        }
        os.write(("\t" + m_tissuemap.size()).getBytes());

        for (int i = 0; i < m_tissuemap.size(); ++i) {
            os.write("\t".getBytes());
        }
        String tissue_string = tissuemap_to_sorted_string("\t");
        os.write(("\nName\tDescription" + tissue_string + "\n").getBytes());

        ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(tissue_string, "\t", false)));

        if (assem == Assembly.primary) {
            for (Map.Entry<String, MapEntry> it : m_mapping.entrySet()) {
                mapping_set.add(it.getValue());
            }
        } else if (assem == Assembly.patchhaploscaff) {
            for (Map.Entry<String, MapEntry> it : m_mapping_phs.entrySet()) {
                mapping_set.add(it.getValue());
            }
        }

        for (MapEntry sit : mapping_set) {
            sit.to_gct(tokens, os, chrincluded);
        }
    }

    //converts all peptides to _ptm.bed lines
    public final void to_ptmbed(String filename, String filename2, Assembly assem) throws Exception {
        to_ptmbed(filename, filename2, assem, true);
    }

    public final void to_ptmbed(String filename, String filename2) throws Exception {
        to_ptmbed(filename, filename2, Assembly.primary, true);
    }

    public final void to_ptmbed(String filename, String filename2, Assembly assem, boolean chrincluded) throws Exception {
        FileOutputStream ofs = new FileOutputStream(filename);
        FileOutputStream ofs2 = new FileOutputStream(filename2);
        to_ptmbed(assem, ofs, ofs2, chrincluded);
        ofs.close();
        ofs2.close();
    }

    public final void to_ptmbed(Assembly assem, OutputStream os, OutputStream os2) throws Exception {
        to_ptmbed(assem, os, os2, true);
    }

    public final void to_ptmbed(Assembly assem, OutputStream os) throws Exception {
        to_ptmbed(assem, os, System.out, true);
    }

    public final void to_ptmbed(Assembly assem) throws Exception {
        to_ptmbed(assem, System.out, System.out, true);
    }

    public final void to_ptmbed(Assembly assem, OutputStream os, OutputStream os2, boolean chrincluded) throws Exception {
        TreeSet<MapEntry> mapping_set = new TreeSet<>(new MapentryPCompare());

        if (assem == Assembly.primary) {
            for (Map.Entry<String, MapEntry> it : m_mapping.entrySet()) {
                mapping_set.add(it.getValue());
            }
        } else if (assem == Assembly.patchhaploscaff) {
            for (Map.Entry<String, MapEntry> it : m_mapping_phs.entrySet()) {
                mapping_set.add(it.getValue());
            }
        }

        for (MapEntry sit : mapping_set) {
            sit.to_ptmbed(os, os2, chrincluded);
        }
    }

    //removes all peptides from the MappedPeptides.
    public final void remove_all_peptides() {
        for (Map.Entry<String, MapEntry> it : m_mapping.entrySet()) {
            it.getValue().remove_peptides();
        }
        for (Map.Entry<String, MapEntry> it : m_mapping_phs.entrySet()) {
            it.getValue().remove_peptides();
        }
        m_tissuemap.clear();
        m_count_peptides = 0;
        m_count_peptides_phs = 0;
        m_tissueindex = 0;
    }

    //converts the above map to a string that contains all found tissues, delimited by '\t' and sorted.
    private String tissuemap_to_sorted_string(String sep) {
        TreeSet<Pair<String, Integer>> tissueset = new TreeSet<>(new ByIntValue());
        for (Map.Entry<String, Integer> it : m_tissuemap.entrySet()) {
            tissueset.add(new Pair<>(it.getKey(), it.getValue()));
        }

        StringBuilder ss = new StringBuilder();
        for (Pair<String, Integer> it : tissueset) {
            ss.append("\t").append(it.getKey());
        }
        return ss.toString();
    }

    //passes the peptide insertion and lookup to the MapEntry.
    public final void add_peptide(CoordinateWrapper coordwrapper, String sequence, String tag, int sigPSMs, int genes, FileOutputStream ofstream, double quant, Map.Entry<String, TranscriptsT> transcriptsEntry) throws Exception {
        if (!m_tissuemap.containsKey(tag)) {
            m_tissuemap.put(tag, m_tissueindex);
            ++m_tissueindex;
        }

        final String geneID = transcriptsEntry.getKey();
        if (!(m_mapping.containsKey(geneID)) && !(m_mapping_phs.containsKey(geneID))) {
            StringJoiner ss = new StringJoiner(",");
            for (Map.Entry<String, ArrayList<PositionMismatchT>> entryIt : transcriptsEntry.getValue().getM_entries().entrySet()) {
                ss.add(entryIt.getKey());
            }
            ofstream.write((geneID + "\t" + sequence + "\t" + ss.toString() + "\t" + genes + "\t" + tag + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
        } else if (m_mapping.containsKey(geneID) && !(m_mapping_phs.containsKey(geneID))) {
            m_count_peptides += (m_mapping.get(geneID).add_peptide(coordwrapper, sequence, tag, sigPSMs, genes, ofstream, quant, transcriptsEntry) != 0) ? 0 : 1;
        } else if (!(m_mapping.containsKey(geneID)) && m_mapping_phs.containsKey(geneID)) {
            m_count_peptides_phs += (m_mapping_phs.get(geneID).add_peptide(coordwrapper, sequence, tag, sigPSMs, genes, ofstream, quant, transcriptsEntry) != 0) ? 0 : 1;
        } else {
            m_count_peptides += (m_mapping.get(geneID).add_peptide(coordwrapper, sequence, tag, sigPSMs, genes, ofstream, quant, transcriptsEntry) != 0) ? 0 : 1;
            m_count_peptides_phs += (m_mapping_phs.get(geneID).add_peptide(coordwrapper, sequence, tag, sigPSMs, genes, ofstream, quant, transcriptsEntry) != 0) ? 0 : 1;
        }
    }

}