package org.bigbio.pgatk.pepgenome.io;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

public abstract class GenomeAnnotationParser {

    // Fields and methods common to GTFParser and GFF3Parser

    //input stream
    BufferedReader reader;

    private FileInputStream ifs;

    //current line
    String line;

    private static Logger log = LoggerFactory.getLogger(GenomeAnnotationParser.class);

    //returns true if in the GTF at position 6 there is a + (plus strand)
    static boolean is_first_strand(List<String> tokens) {
        return tokens.get(6).equals("+");
    }

    //returns true if position 2 in the GTF says "CDS"
    static boolean is_cds(List<String> tokens) {
        return tokens.get(2).equals("CDS");
    }

    //returns true if position 2 in the GTF says "exon"
    static boolean is_exon(List<String> tokens) {
        return tokens.get(2).equals("exon");
    }

    //returns true if position 2 in the GTF/GFF3 says "transcript" or "mRNA"
    static boolean is_next_transcript(List<String> tokens) {
        return (tokens.get(2).equals("transcript") || tokens.get(2).equals("mRNA"));
    }

    // Returns true if line feature type is gene
    static boolean is_next_gene(List<String> tokens) {
        return tokens.get(2).equals("gene");
    }

    //opens filestream, returns true if successful
    boolean open(String file) throws Exception {
        if (reader == null) {
            line = "";
            ifs = new FileInputStream(file);
            reader = new BufferedReader(new InputStreamReader(ifs));
        }
        boolean status = true;
        try{
            status = reader.ready();
        }catch (IOException ex){
            log.debug("The genome annotation file stream is closed -- " + file);
            ifs = new FileInputStream(file);
            reader = new BufferedReader(new InputStreamReader(ifs));
        }
        return status;
    }

    //closes the filestream
    void close() throws Exception {
        line = "";
        reader.close();
        ifs.close();
    }


}
