package org.bigbio.pgatk.pepgenome;

import org.bigbio.pgatk.pepgenome.common.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

//this fastaparser reads fasta files.
//it will convert input sequences into iso-sequences [I, L] will be converted to J
public class FastaParser {

    //inputstream.
    private BufferedReader br = null;
    //current line.
    private String m_line = "";
    //meyers singleton instance.
    private static FastaParser m_instance = null;

    private FastaParser() {
    }

    //singleton get_instance method.
    public static FastaParser get_instance() {
        if (m_instance == null) {
            m_instance = new FastaParser();
        }
        return m_instance;
    }

    //opens the file
    public boolean open(String file) throws Exception {
        if (br == null) {
            br = new BufferedReader(new FileReader(file));
            m_line = "";
        }
        return br.ready();
    }

    //closes the filestream
    public final void close() {
        m_line = "";
        try {
            br.close();
        } catch (IOException e) {
            //e.printStackTrace();
        }
        br = null;
    }

    //parses and returns the next FASTA entry.
    public final FastaEntry next_entry() throws IOException {
        if (m_line.isEmpty()) {
            m_line = br.readLine();
        }
        if (m_line == null) {
            return new FastaEntry("", "");
        }
        String header = m_line;
        StringBuilder sequenceBuilder = new StringBuilder();
        while ((m_line = br.readLine()) != null && !m_line.startsWith(">")) {
            sequenceBuilder.append(Utils.make_iso_sequence(m_line));
        }
        if (m_line == null || !m_line.startsWith(">")) {
            m_line = "";
        }
        return new FastaEntry(header, sequenceBuilder.toString());
    }
}