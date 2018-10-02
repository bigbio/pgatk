package org.bigbio.pgatk.pepgenome.kmer;

import org.bigbio.pgatk.pepgenome.common.ProteinEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;

import java.util.ArrayList;
import java.util.Map;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 *
 * @author ypriverol on 01/10/2018.
 */
public interface IKmerMap {


    /**
     * Digests and adds a protein to the map.
     * @param protein Protein Entry
     */
    void add_protein(ProteinEntry protein);

    /**
     * Searches for all matches (imperfect matching, set via PEPTIDE_MAPPER) and returns a map that contains all finds.
     * @param peptide_string Peptide String to search
     * @return Transcript map
     */
    Map<String, TranscriptsT> find_peptide(String peptide_string);

    /**
     * Inserts a found peptide into the current gene id map.
     * @param entry Entry
     * @param mismatches Number of missmatches
     */
    void insert_into_gene_id_map(IKmerEntry entry, ArrayList<Integer> mismatches);

    /**
     * Inserts a found peptide into the current gene id map.
     * @param entry KmerEntry
     * @param mismatches Number of mismatches accepted
     * @param offset offset
     */
    void insert_into_gene_id_map(IKmerEntry entry, ArrayList<Integer> mismatches, int offset);

    /**
     * Returns true if a kmer (key) is in the digested proteins
     * @param key find a kmer in the Map
     * @return True if the key is in the map
     */
    boolean contains(String key);

    /**
     * Returns the number of fragments that were created during digestion
     * @return Number of fragments (peptides)
     */
    int size();

}
