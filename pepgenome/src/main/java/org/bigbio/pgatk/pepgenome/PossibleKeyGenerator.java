package org.bigbio.pgatk.pepgenome;

import org.bigbio.pgatk.pepgenome.common.Tuple;
import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class PossibleKeyGenerator implements Serializable {
    private static final long serialVersionUID = 809902192340291246L;
    //string used to generate keys.
    private String m_key = "";

    //are there any generated keys left?
    private boolean m_keys_generated = false;

    // kmermap, to check if a certain key exists in there,
    //otherwise the key will not be added to the generated keys.
    private IKmerMap m_kmers;

    //current set of generated keys.
    private List<String> m_keys = new ArrayList<>();

    //iterator pointing to the current element in m_keys.
    private int m_curr_index;

    public PossibleKeyGenerator(IKmerMap k) {
        this.m_kmers = k;
    }

    //returns
    //-1: key is shorter than KMER_LENGTH
    // 0: key is shorter than (ALLOWED MISMATCHES + 1) * KMER_LENGTH
    // 1: key is >= KMER_LENGTH
    public final int set_original_key(String key) {
        if (key.length() >= ((GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES + 1) * GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH)) {
            //if key.length >= (allowed_mismatches+1)*kmer_lenght
            //in this case there is no way there are more than allowed_mismatches in every kmer.
            //one of those has to be matching perfectly:
            //protein:  TESTTESTTEST allowed mismatches = 2
            //key:	    TEXTTEXTTEST  length = 15, 2 mismatches on the protein
            //mismatches +1 * kmerlength = 15 therefore this if will be used
            //splitting the key into TEXT, TEXT, and TEST. the TEST does match the protein perfectly.
            //therefore the peptide will be found. this uses the match_backwards function.
            //if the key is longer this approach will also work. in the 2 mismatches case it generates 3 keys,
            //as opposed to the brute force key generation which will produce 3331 unique keys.
            for (int i = 0; i <= GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES; i++) {
                int start = i * GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH;
                int end = ((start + GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH) > key.length()) ? key.length() : (start + GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH);
                m_keys.add(key.substring(start, end));
            }
            m_keys_generated = true;
            m_curr_index = 0;
            return 1;
        }
        if (key.length() >= GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH) {
            set_short_original_key(key.substring(0, GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH));
            return 0;
        }
        return -1;
    }

    //used to set the kmerlength long key.
    private void set_short_original_key(String key) {
        if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES > 0) {
            m_key = key;
            m_keys_generated = false;
        } else {
            m_keys.add(key);
            m_keys_generated = true;
            m_curr_index = 0;
        }
    }

    //assigns the next key to the passed string.
    //returns true if there is another key otherwise returns false.
    //the first call might be a bit slower as the keys are generated in a lazy way.
    public String get_next_key() {
        if (!m_keys_generated) {
            //if the keys are not yet generated
            if (m_key.length() > 0 && GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS.length > 0) {
                //generate them and return the first one (lazy generation)
                generate_keys();
                m_keys_generated = true;
                m_curr_index = 0;
                return get_next_key();
            }
        } else {
            //otherwise just return the next key
            if (m_curr_index < m_keys.size()) {
                String key = m_keys.get(m_curr_index);
                ++m_curr_index;
                return key;
            }
            //and if the end is reached, clear out the old keys and wait for the next peptide.
            m_keys_generated = false;
            m_key = "";
            m_keys.clear();
        }
        return null;
    }

    //generates keys for one mismatch matching.
    private void generate_keys_one_mismatch() {
        //generates all keys with one mismatch.
        for (int kmer_it = 0; kmer_it < GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH; kmer_it++) {
            if (kmer_it >= m_key.length()) {
                break;
            }
            for (int aa_it = 0; aa_it < GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS.length; aa_it++) {
                StringBuilder tmp = new StringBuilder(m_key);
                tmp.setCharAt(kmer_it, GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS[aa_it]);
                String tmpStr = tmp.toString();
                if (m_kmers.contains(tmpStr)) {
                    m_keys.add(tmpStr);
                }
            }
        }
    }

    private List<Tuple<Character, Character>> generateKeysTwoMismatchesCombinations = new ArrayList<>();
    private boolean generateKeysTwoMismatchesCombinationsGenerated = false;

    //generates keys for two mismatch matching.
    private void generate_keys_two_mismatches() {
        //generates all keys with 2 mismatches. this function is costly and should be avoided if possible.
        if (!generateKeysTwoMismatchesCombinationsGenerated) {
            for (int i = 0; i < GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS.length; ++i) {
                for (int j = 0; j < GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS.length; j++) {
                    generateKeysTwoMismatchesCombinations.add(new Tuple<>(GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS[i], GenomeMapper.PEPTIDE_MAPPER.ALLOWED_AMINO_ACIDS[j]));
                }
            }
            generateKeysTwoMismatchesCombinationsGenerated = true;
        }

        for (int i = 0; i < GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH - 1; i++) {
            if (i >= m_key.length()) {
                break;
            }
            for (int j = i; j < GenomeMapper.PEPTIDE_MAPPER.KMER_LENGTH; j++) {
                if (j >= m_key.length()) {
                    break;
                }
                StringBuilder tmp = new StringBuilder(m_key);
                for (Tuple<Character, Character> pair : generateKeysTwoMismatchesCombinations) {
                    char first = pair.getKey();
                    char second = pair.getValue();
                    tmp.setCharAt(i, first);
                    tmp.setCharAt(j, second);
                    String tmpStr = tmp.toString();
                    if (m_kmers.contains(tmpStr)) {
                        m_keys.add(tmpStr);
                    }
                }
            }
        }
    }

    //calls the other generator functions.
    private void generate_keys() {
        if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES == 1) {
            generate_keys_one_mismatch();
        } else if (GenomeMapper.PEPTIDE_MAPPER.ALLOWED_MISMATCHES == 2) {
            if (GenomeMapper.PEPTIDE_MAPPER.ONE_IN_FIVE_MODE) {
                generate_keys_one_mismatch();
            } else {
                generate_keys_two_mismatches();
            }
        }
        //add to here if three mismatches are required at some point.
    }

}
