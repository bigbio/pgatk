package bigbio.pgatk.jpogo.common;

import bigbio.pgatk.jpogo.PeptideEntry;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class ExistingPeptides {

    private Map<String, ArrayList<PeptideEntry>> m_existing_peptides = new TreeMap<>();

    public ExistingPeptides() {
    }

    //evaluates if a certain peptidestring is already present
    public final boolean contains(String peptideString) {
        return m_existing_peptides.containsKey(peptideString);
    }

    // add an entry. creates the entry in the map if it does not already exist
    public final void add(String peptideString, PeptideEntry peptideEntry) {
        ArrayList<PeptideEntry> peptideEntries = m_existing_peptides.computeIfAbsent(peptideString, k -> new ArrayList<>());
        peptideEntries.add(peptideEntry);
    }

    //access to the elements that are saved in the map.
    public final ArrayList<PeptideEntry> getItem(String peptideString) {
        return m_existing_peptides.get(peptideString);
    }
}
