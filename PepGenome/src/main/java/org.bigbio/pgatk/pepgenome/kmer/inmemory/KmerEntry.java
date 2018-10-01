package org.bigbio.pgatk.pepgenome.kmer.inmemory;

import lombok.Data;
import org.bigbio.pgatk.pepgenome.ProteinEntry;
import org.bigbio.pgatk.pepgenome.kmer.IKmerEntry;

/**
 * A kmer Entry holds information on a kmer. It has pointers to the proteinsequence and associated
 * protein to which it belongs, and knows its position therein.
 */

public class KmerEntry implements IKmerEntry {

    // Full Protein entry
    ProteinEntry m_p_protein;

    //the (0 based) index of the first letter of the kmer in the protein string.
    int m_pos_in_protein;

    public KmerEntry(String key, ProteinEntry protein, int pos_in_protein) {
        this.m_p_protein = protein;
        this.m_pos_in_protein = pos_in_protein;
    }

    public KmerEntry() {
        this.m_p_protein = null;
        this.m_pos_in_protein = 0;
    }

    @Override
    public ProteinEntry m_p_protein() {
        return null;
    }

    @Override
    public int m_pos_in_protein() {
        return 0;
    }
}
///#endif
