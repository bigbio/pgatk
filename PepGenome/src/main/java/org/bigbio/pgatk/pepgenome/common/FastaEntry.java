package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;

public class FastaEntry implements Serializable {

    private static final long serialVersionUID = 2153797062207182214L;
    //holds the fasta header
    private String m_header;
    //holds the fasta entry sequence.
    private String m_aa_sequence;

    public FastaEntry(String header, String AAsequence) {
        this.m_header = header;
        this.m_aa_sequence = AAsequence;
    }

    public FastaEntry() {
        this.m_header = "";
        this.m_aa_sequence = "";
    }

    public final String get_header() {
        return m_header;
    }

    public final String get_sequence() {
        return m_aa_sequence;
    }

    public final boolean is_empty() {
        return m_header.isEmpty() && m_aa_sequence.isEmpty();
    }
}