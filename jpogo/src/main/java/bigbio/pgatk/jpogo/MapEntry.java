package bigbio.pgatk.jpogo;

import bigbio.pgatk.jpogo.common.TranscriptsT;
import bigbio.pgatk.jpogo.common.Utils;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.*;

//a map entry is used to connect a peptide string to all its variations
//it maps a peptide to a gene, and the transcripts.
public class MapEntry implements Comparable<MapEntry> {

    //pointer to the associated GeneEntry
    private GeneEntry m_p_gene_entry;
    //peptideentries, maps sequence without ptms to the corresponding PeptideEntry
    private Map<String, PeptideEntry> m_peptide_entries = new TreeMap<>();
    //holds all transcripts for this MapEntry
    private Set<String> m_transcripts = new TreeSet<>();

    @Override
    public int compareTo(MapEntry otherInstance) {
        if (lessThan(otherInstance)) {
            return -1;
        } else if (otherInstance.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    public MapEntry(GeneEntry geneentry_p) {
        this.m_p_gene_entry = geneentry_p;
//        this.m_peptide_entries = new TreeMap<>();
//        this.m_transcripts = new TreeSet<>();
    }

    //adds a transcript id to the mapping.
    public final void add_transcript_id(String transcriptID) {
        m_transcripts.add(transcriptID);
    }

    //compares two MapEntry objects. returns true if the lhs 'GeneEntry is lesser than rhs'
    public boolean lessThan(MapEntry rhs) {
        return m_p_gene_entry.lessThan(rhs.m_p_gene_entry);
    }

    //calls the PeptideEntry::to_gtf metod for every peptide.
    public final OutputStream to_gtf(String source, OutputStream os) throws Exception {
        return to_gtf(source, os, true);
    }

    public final OutputStream to_gtf(String source) throws Exception {
        return to_gtf(source, System.out, true);
    }

    public final OutputStream to_gtf(String source, boolean chrincluded) throws Exception {
        return to_gtf(source, System.out, chrincluded);
    }

    public final OutputStream to_gtf(String source, OutputStream os, boolean chrincluded) throws Exception {
        if (m_peptide_entries.size() > 0) {
            m_p_gene_entry.to_gtf(source, os).write("\n".getBytes());

            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : m_peptide_entries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_gtf(source, os).write("\n".getBytes());
            }
        }
        return os;
    }

    //calls the PeptideEntr::to_bed metod for every peptide.
    public final OutputStream to_bed(OutputStream os) throws Exception {
        return to_bed(os, true);
    }

    public final OutputStream to_bed() throws Exception {
        return to_bed(System.out, true);
    }

    public final OutputStream to_bed(OutputStream os, boolean chrincluded) throws Exception {
        if (m_peptide_entries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : m_peptide_entries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_bed(os);
            }
        }
        return os;
    }

    //calls the PeptideEntry::to_gct metod for every peptide.
    public final OutputStream to_gct(ArrayList<String> tissuelist, OutputStream os) throws Exception {
        return to_gct(tissuelist, os, true);
    }

    public final OutputStream to_gct(ArrayList<String> tissuelist) throws Exception {
        return to_gct(tissuelist, System.out, true);
    }

    public final OutputStream to_gct(ArrayList<String> tissuelist, OutputStream os, boolean chrincluded) throws Exception {
        if (m_peptide_entries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : m_peptide_entries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_gct(m_p_gene_entry.get_id(), tissuelist, os);
            }
        }
        return os;
    }

    //calls the PeptideEntry::to_ptmbed metod for every peptide.
    public final OutputStream to_ptmbed(OutputStream os, OutputStream os2) throws Exception {
        return to_ptmbed(os, os2, true);
    }

    public final OutputStream to_ptmbed(OutputStream os) throws Exception {
        return to_ptmbed(os, System.out, true);
    }

    public final OutputStream to_ptmbed() throws Exception {
        return to_ptmbed(System.out, System.out, true);
    }

    public final OutputStream to_ptmbed(OutputStream os, OutputStream os2, boolean chrincluded) throws Exception {
        if (m_peptide_entries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : m_peptide_entries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_ptmbed(os);
                if (pit.noPTM()) {
                    pit.to_bed(os2, true);
                }
            }
        }
        return os;
    }

    //removes all peptides that are associated with a specific sequence.
    public final void remove_peptides() {
        m_peptide_entries.clear();
    }

    //delegates to peptide_entry after checking if that peptide already
    //exists and creating it if it doesn't
    public final int add_peptide(CoordinateWrapper coordwrapper, String sequence, String tag, int sigPSMs, int genes, FileOutputStream ofstream, double quant, Map.Entry<String, TranscriptsT> transcriptsEntry) {
        int newly_added = 0;
        String sequence_wo_ptm = Utils.remove_ptms(sequence);

        if (transcriptsEntry.getValue().getM_entries().size() > 0) {
            if (!m_peptide_entries.containsKey(sequence_wo_ptm)) {
                PeptideEntry new_peptide = new PeptideEntry(m_p_gene_entry);
                m_peptide_entries.put(sequence_wo_ptm, new_peptide);
                coordwrapper.add_to_existing_peptides(sequence_wo_ptm, new_peptide);
                newly_added = 1;
            }
            m_peptide_entries.get(sequence_wo_ptm).add_peptide(coordwrapper, sequence_wo_ptm, sequence, tag, sigPSMs, transcriptsEntry.getValue(), genes, ofstream, quant);
        }
        return newly_added;
    }

}