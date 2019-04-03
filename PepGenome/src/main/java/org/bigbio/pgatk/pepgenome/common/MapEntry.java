package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.comparators.PeptideentryPcompare;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.*;

//a map entry is used to connect a peptide string to all its variations
//it maps a peptide to a gene, and the transcripts.
public class MapEntry implements Comparable<MapEntry>, Serializable {

    private static final long serialVersionUID = 6404015419072151797L;

    //pointer to the associated GeneEntry
    private GeneEntry geneEntry;

    //peptideentries, maps sequence without ptms to the corresponding PeptideEntry
    private Map<String, PeptideEntry> peptideEntries = new TreeMap<>();

    //holds all transcripts for this MapEntry
    private Set<String> transcripts = new TreeSet<>();

    @Override
    public int compareTo(MapEntry o) {
        if (lessThan(o)) {
            return -1;
        } else if (o.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    public MapEntry(GeneEntry geneentry_p) {
        this.geneEntry = geneentry_p;
//        this.peptideEntries = new TreeMap<>();
//        this.transcripts = new TreeSet<>();
    }

    //adds a transcript id to the mapping.
    public final void addTranscriptId(String transcriptID) {
        transcripts.add(transcriptID);
    }

    //compares two MapEntry objects. returns true if the lhs 'GeneEntry is lesser than rhs'
    public boolean lessThan(MapEntry rhs) {
        return geneEntry.isLessThan(rhs.geneEntry);
    }

    //calls the PeptideEntry::toGtf metod for every peptide.
    public final OutputStream toGtf(String source, OutputStream os) throws Exception {
        return toGtf(source, os, true);
    }

    public final OutputStream toGtf(String source) throws Exception {
        return toGtf(source, System.out, true);
    }

    public final OutputStream toGtf(String source, boolean chrincluded) throws Exception {
        return toGtf(source, System.out, chrincluded);
    }

    public final OutputStream toGtf(String source, OutputStream os, boolean chrincluded) throws Exception {
        if (peptideEntries.size() > 0) {
            geneEntry.to_gtf(source, os).write("\n".getBytes());

            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : peptideEntries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_gtf(source, os).write("\n".getBytes());
            }
        }
        return os;
    }

    //calls the PeptideEntr::toBed metod for every peptide.
    public final OutputStream toBed(OutputStream os) throws Exception {
        return toBed(os, true);
    }

    public final OutputStream toBed() throws Exception {
        return toBed(System.out, true);
    }

    public final OutputStream toBed(OutputStream os, boolean chrincluded) throws Exception {
        if (peptideEntries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : peptideEntries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_bed(os);
            }
        }
        return os;
    }

    //calls the PeptideEntry::toGct metod for every peptide.
    public final OutputStream toGct(ArrayList<String> tissuelist, OutputStream os) throws Exception {
        return toGct(tissuelist, os, true);
    }

    public final OutputStream toGct(ArrayList<String> tissuelist) throws Exception {
        return toGct(tissuelist, System.out, true);
    }

    public final OutputStream toGct(ArrayList<String> tissuelist, OutputStream os, boolean chrincluded) throws Exception {
        if (peptideEntries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : peptideEntries.entrySet()) {
                peptide_entries_set.add(it.getValue());
            }

            for (PeptideEntry pit : peptide_entries_set) {
                pit.to_gct(geneEntry.get_id(), tissuelist, os);
            }
        }
        return os;
    }

    //calls the PeptideEntry::toPtmbed metod for every peptide.
    public final OutputStream toPtmbed(OutputStream os, OutputStream os2) throws Exception {
        return toPtmbed(os, os2, true);
    }

    public final OutputStream toPtmbed(OutputStream os) throws Exception {
        return toPtmbed(os, System.out, true);
    }

    public final OutputStream toPtmbed() throws Exception {
        return toPtmbed(System.out, System.out, true);
    }

    public final OutputStream toPtmbed(OutputStream os, OutputStream os2, boolean chrincluded) throws Exception {
        if (peptideEntries.size() > 0) {
            TreeSet<PeptideEntry> peptide_entries_set = new TreeSet<>(new PeptideentryPcompare());
            for (Map.Entry<String, PeptideEntry> it : peptideEntries.entrySet()) {
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
    public final void removePeptides() {
        peptideEntries.clear();
    }

    /**
     * Delegates to peptide_entry after checking if that peptide already
     * exists and creating it if it doesn't
     *
     */

    public final int addPeptide(CoordinateWrapper coordwrapper, String sequence, String tag, int sigPSMs, int genes, FileOutputStream ofstream, double quant, Map.Entry<String, TranscriptsT> transcriptsEntry) {
        int added = 0;
        String sequenceWoPtm = Utils.remove_ptms(sequence);

        if (transcriptsEntry.getValue().getM_entries().size() > 0) {
            if (!peptideEntries.containsKey(sequenceWoPtm)) {
                PeptideEntry newPeptide = new PeptideEntry(geneEntry);
                peptideEntries.put(sequenceWoPtm, newPeptide);
                coordwrapper.add_to_existing_peptides(sequenceWoPtm, newPeptide);
                added = 1;
            }
            peptideEntries.get(sequenceWoPtm).add_peptide(coordwrapper, sequenceWoPtm, sequence, tag, sigPSMs, transcriptsEntry.getValue(), genes, ofstream, quant);
        }
        return added;
    }

}