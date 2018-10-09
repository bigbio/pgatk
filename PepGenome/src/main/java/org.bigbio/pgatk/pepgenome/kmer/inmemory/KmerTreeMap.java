package org.bigbio.pgatk.pepgenome.kmer.inmemory;

import lombok.extern.slf4j.Slf4j;
import org.apache.commons.lang3.ArrayUtils;
import org.bigbio.pgatk.pepgenome.PossibleKeyGenerator;
import org.bigbio.pgatk.pepgenome.common.ProteinEntry;
import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;
import org.bigbio.pgatk.pepgenome.common.PositionMismatchT;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.kmer.IKmerEntry;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import org.bigbio.pgatk.pepgenome.kmer.KmerEntry;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

@Slf4j
public class KmerTreeMap implements IKmerMap, Serializable {

    private static final long serialVersionUID = 4276302531650610935L;
    //kmerMap
    //in the gencode fasta file, data analysis showed that the average size of the list is 18 elements (for approx. 93 000 proteins)
    //and the map size is at approx 1.8m elements.
    private Map<String, IKmerEntry[]> m_kmers = new TreeMap<>();

    //gene_id map
    private Map<String, TranscriptsT> m_gene_id_map = new TreeMap<>();

    //the keygenerator will generate keys for one and two mismatch matching.
    private PossibleKeyGenerator m_key_gen;

    private ProteinMatcher proteinMatcher;

    public KmerTreeMap() {

        this.m_key_gen = new PossibleKeyGenerator(this);
        if ((GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES > 1) && GenomeMapper.PEPTIDE_MAPPER.ONE_IN_FIVE_MODE) {
            proteinMatcher = new MatcherOneInFiveMode();
        } else {
            proteinMatcher = new MatcherNormal();
        }



    }

    //digests and adds a protein to the map.
    public void add_protein(ProteinEntry protein) {
        //this function digests a protein sequence and builds a KmerTreeMap with the bits.
        String protein_sequence = protein.get_sequence();

        if (protein_sequence.length() >= GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH) {
            //<= n!
            int n = protein_sequence.length() - GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH;
            String key;
            IKmerEntry[] kmerEntries;
            KmerEntry curr_kmer;

            for (int i = 0; i <= n; i++) {
                key = Utils.getCppStyleSubString(protein_sequence, i, GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH);
                //if the kmer doesnt exist yet, create kmerEntries
                kmerEntries = m_kmers.computeIfAbsent(key, k -> new KmerEntry[0]);
                curr_kmer = new KmerEntry(Utils.getCppStyleSubStringByShift(protein_sequence, i), protein, i);
                //add the current KmerEntry
                m_kmers.put(key, ArrayUtils.add(kmerEntries, curr_kmer));
            }
        }
    }

    //searches for all matches (imperfect matching, set via PEPTIDE_MAPPER) and returns a map that contains all finds.
    public final Map<String, TranscriptsT> find_peptide(String peptide_string) {
        //this function generates a gene_id_map.
        //this map will map a gene_id to all related transcript ids and all the peptides and their
        //position in the proteinsequence

        if (!m_gene_id_map.isEmpty()) {
            m_gene_id_map.clear();
        }

        int set_key_returned = m_key_gen.set_original_key(peptide_string);
        int backwards_multiplier = 0;

        String curr_key;
        ArrayList<Integer> mismatches = new ArrayList<>();

        int peptide_length = peptide_string.length();

        if (set_key_returned >= 0) {
            while ((curr_key = m_key_gen.get_next_key()) != null) {
                IKmerEntry[] kmerEntries = m_kmers.get(curr_key);
                if (kmerEntries != null) {
                    for (IKmerEntry entry : kmerEntries) {
                        if (set_key_returned == 0) {
                            if (proteinMatcher.match(peptide_string, entry, mismatches, peptide_length)) {
                                insert_into_gene_id_map(entry, mismatches);
                            }
                        } else if (set_key_returned == 1) {
                            //this mode is used when only allowed_mismatches + 1 keys are generated.
                            //see PossibleKeyGenerator::set_original_key
                            int offset = backwards_multiplier * GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH;
                            if (proteinMatcher.match_backwards(peptide_string, entry, mismatches, peptide_length, offset)) {
                                insert_into_gene_id_map(entry, mismatches, offset);
                            }
                        }
                        mismatches.clear();
                    }
                }
                backwards_multiplier++;
            }
        }
        return m_gene_id_map;
    }

    //inserts a found peptide into the current gene id map.
    public void insert_into_gene_id_map(IKmerEntry entry, ArrayList<Integer> mismatches) {
        insert_into_gene_id_map(entry, mismatches, 0);
    }

    //inserts a found peptide into the current gene id map.
    public void insert_into_gene_id_map(IKmerEntry entry, ArrayList<Integer> mismatches, int offset) {
        //inserts a found position into the gene id map
        int pos_in_protein = entry.m_pos_in_protein() - offset;
        String gene_id = entry.m_p_protein().get_gene_id();
        String transcript_id = entry.m_p_protein().get_transcript_id();

        ArrayList<PositionMismatchT> pMismatchTS = m_gene_id_map.computeIfAbsent(gene_id, k -> new TranscriptsT())
                .getM_entries().computeIfAbsent(transcript_id, j -> new ArrayList<>());

        //care: if the PositionMismatchT struct is changed to accomodate more than 2 mismatches this has to be updated as well.
        pMismatchTS.add(new PositionMismatchT(
                pos_in_protein,
                (mismatches.size() > 0) ? pos_in_protein + mismatches.get(0) : -1,
                (mismatches.size() > 1) ? pos_in_protein + mismatches.get(1) : -1));
    }

    // returns true if a kmer (key) is in the digested proteins
    public final boolean contains(String key) {
        return m_kmers.containsKey(key);
    }

    //returns the number of fragments that were created during digestion
    public final int size() {
        return m_kmers.size();
    }


    interface ProteinMatcher {
        boolean match(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength);

        boolean match_backwards(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength, int offset);
    }

    class MatcherNormal implements ProteinMatcher {

        @Override
        //the basic forward matching algorithm
        public boolean match(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength) {
            int protein_length = kmerEntry.m_p_protein().get_sequence().length() - kmerEntry.m_pos_in_protein();
            if (peptideLength <= protein_length) {
                String p_protein = Utils.getCppStyleSubStringByShift(kmerEntry.m_p_protein().get_sequence(), kmerEntry.m_pos_in_protein());
                for (int i = 0; i < peptideLength; i++) {
                    if (peptideString.charAt(i) != p_protein.charAt(i)) {
                        mismatches.add(i);
                        if (mismatches.size() > GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                            return false;
                        }
                    }
                }
                return true;
            }
            return false;
        }

        @Override
        //backwards matching functionality. this is possible by using cstrings.
        public boolean match_backwards(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength, int offset) {
            if (kmerEntry.m_pos_in_protein() >= offset) {
                int protein_length = kmerEntry.m_p_protein().get_sequence().length() - kmerEntry.m_pos_in_protein() + offset;
                if (peptideLength <= protein_length) {
                    String p_protein = Utils.getCppStyleSubStringByShift(kmerEntry.m_p_protein().get_sequence(), (kmerEntry.m_pos_in_protein() - offset));
                    for (int i = 0; i < peptideLength; i++) {
                        if (peptideString.charAt(i) != p_protein.charAt(i)) {
                            mismatches.add(i);
                            if (mismatches.size() > GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                                return false;
                            }
                        }
                    }
                    return true;
                }
            }
            return false;
        }
    }

    class MatcherOneInFiveMode implements ProteinMatcher {

        @Override
        //forward matching with one in five stop criterion
        public boolean match(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength) {
            int protein_length = kmerEntry.m_p_protein().get_sequence().length() - kmerEntry.m_pos_in_protein();
            if (peptideLength <= protein_length) {
                String p_protein = Utils.getCppStyleSubStringByShift(kmerEntry.m_p_protein().get_sequence(), kmerEntry.m_pos_in_protein());
                for (int i = 0; i < peptideLength; i++) {
                    if (peptideString.charAt(i) != p_protein.charAt(i)) {
                        if (mismatches.size() != 0) {
                            if ((mismatches.get(mismatches.size() - 1) - i) <= GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                                return false;
                            }
                        }
                        mismatches.add(i);
                        if (mismatches.size() > GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                            return false;
                        }
                    }
                }
                return true;
            }
            return false;
        }

        @Override
        //backwards matching functionality with one in five stop criterion.
        public boolean match_backwards(String peptideString, IKmerEntry kmerEntry, ArrayList<Integer> mismatches, int peptideLength, int offset) {
            if (kmerEntry.m_pos_in_protein() >= offset) {
                int protein_length = kmerEntry.m_p_protein().get_sequence().length() - kmerEntry.m_pos_in_protein() + offset;
                if (peptideLength <= protein_length) {
                    String p_protein = Utils.getCppStyleSubStringByShift(kmerEntry.m_p_protein().get_sequence(), (kmerEntry.m_pos_in_protein() - offset));
                    for (int i = 0; i < peptideLength; i++) {
                        if (peptideString.charAt(i) != p_protein.charAt(i)) {
                            if (mismatches.size() != 0) {
                                if ((mismatches.get(mismatches.size() - 1) - i) <= GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                                    return false;
                                }
                            }
                            mismatches.add(i);
                            if (mismatches.size() > GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES) {
                                return false;
                            }
                        }
                    }
                    return true;
                }
            }
            return false;
        }
    }
}
