package org.bigbio.pgatk.pepgenome.common;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

public class TranscriptsT {
    //transcript id, mismatch positions
    private Map<String, ArrayList<PositionMismatchT>> m_entries = new TreeMap<>();

    public Map<String, ArrayList<PositionMismatchT>> getM_entries() {
        return m_entries;
    }

    public void setM_entries(Map<String, ArrayList<PositionMismatchT>> m_entries) {
        this.m_entries = m_entries;
    }
}
