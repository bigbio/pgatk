package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class ExistingPeptides implements Serializable {

    private Map<String, ArrayList<PeptideEntry>> m_existing_peptides = new TreeMap<>();

    public ExistingPeptides() {
    }

    /**
     * Evaluates if a certain PeptideString is already present
     * @param peptideString
     * @return true if peptide exists
     */
    public final boolean contains(String peptideString) {
        return m_existing_peptides.containsKey(peptideString);
    }

    /**
     *
     */
    public final void add(String peptideString, PeptideEntry peptideEntry) {
        ArrayList<PeptideEntry> peptideEntries = m_existing_peptides.computeIfAbsent(peptideString, k -> new ArrayList<>());
        peptideEntries.add(peptideEntry);
    }

    //access to the elements that are saved in the map.
    public final ArrayList<PeptideEntry> getItem(String peptideString) {
        return m_existing_peptides.get(peptideString);
    }
}
