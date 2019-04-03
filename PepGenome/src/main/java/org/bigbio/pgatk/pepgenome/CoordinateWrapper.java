package org.bigbio.pgatk.pepgenome;

import org.bigbio.pgatk.pepgenome.common.ExistingPeptides;
import org.bigbio.pgatk.pepgenome.common.FastaEntry;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.ProteinEntry;
import org.bigbio.pgatk.pepgenome.io.FastaParser;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class CoordinateWrapper implements Serializable {

    private static final long serialVersionUID = -5604555402952311335L;
    //holds fasta headers and the associated ProteinEntry objects.
    private Map<String, ProteinEntry> m_map;

    //pointer to the existing_peptides
    //common.ExistingPeptides holds information about previously read peptides.
    private ExistingPeptides m_existing_peptides;

    private int totalAACount = 0;

    public CoordinateWrapper() {
        this.m_map = new TreeMap<>();
        this.m_existing_peptides = new ExistingPeptides();
    }

    public final int size() {
        return m_map.size();
    }

    /**
     * Looks up a ProteinEntry given a fasta header and returns a reference
     * to that entry.
     * @param transcriptId Transcript Identifier
     * @return ProteinEntry
     */
    public final ProteinEntry lookup_entry(String transcriptId) {
        return m_map.computeIfAbsent(transcriptId, k -> new ProteinEntry());
    }

    /**
     * Adds a ProteinEntry.
     * @param entry
     */
    public final void add(ProteinEntry entry) {
        m_map.put(entry.get_transcript_id(), entry);
    }

    //reads and parses a fasta file and adds all of them to the CoordinateWrapper.
    public final void read_fasta_file(String file) throws Exception {
        totalAACount = 0;
        FastaParser fastaParserSingleton = FastaParser.get_instance();
        if (!fastaParserSingleton.open(file)) {
            throw new IllegalStateException("Problem while reading Fasta file");
        }

        FastaEntry fastaEntry;
        while (!(fastaEntry = fastaParserSingleton.nextEntry()).is_empty()) {
            add(new ProteinEntry(fastaEntry));
            totalAACount += fastaEntry.get_sequence().length();
        }

        fastaParserSingleton.close();
    }

    //adds all previously added proteins to the given KmerTreeMap.
    public final void add_all_proteins_to_kmer_map(IKmerMap kmerMap) {
        for (ProteinEntry entry : m_map.values()) {
            kmerMap.add_protein(entry);
        }
    }

    /**
     *  Adds a peptide to the existing peptides list. this is used in the TabInputPeptideFileParser so
     *  that already found peptides dont have to be mapped again.
     *
     * @param peptideSequence
     * @param peptideEntry
     */
    public final void add_to_existing_peptides(String peptideSequence, PeptideEntry peptideEntry) {
        m_existing_peptides.add(peptideSequence, peptideEntry);
    }

    //gets a reference to the already existing peptides
    //used for adding PTMs and tags.
    public final ArrayList<PeptideEntry> get_existing_peptides_at(String peptideSequence) {
        return m_existing_peptides.getItem(peptideSequence);
    }

    //returns true if the peptide was found before.
    public final boolean isPeptidePresent(String peptideSequence) {
        return m_existing_peptides.contains(peptideSequence);
    }

    /**
     * Return the protein size sum of all proteins in the fasta file
     * @return total AA count
     */
    public int getTotalAACount() {
        return totalAACount;
    }

    /**
     * Get the number of proteins
     * @return Number of proteins
     */
    public int getNumberOfProteins(){
        return m_map.size();
    }




}