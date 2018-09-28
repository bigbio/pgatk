package bigbio.pgatk.jpogo;

//a kmer Entry holds information on a kmer.
//it has pointers to the proteinsequence and associated protein to which it belongs, 
//and knows its position therein
public class KmerEntry {
    //the pointer to the first letter of the key in the protein sequence
    public String m_p_key;
    //the pointer to the protein.
    public ProteinEntry m_p_protein;
    //the (0 based) index of the first letter of the kmer in the protein string.
    public int m_pos_in_protein;

    public KmerEntry(String key, ProteinEntry protein, int pos_in_protein) {
        this.m_p_key = key;
        this.m_p_protein = protein;
        this.m_pos_in_protein = pos_in_protein;
    }

    public KmerEntry() {
        this.m_p_key = null;
        this.m_p_protein = null;
        this.m_pos_in_protein = 0;
    }
}
///#endif
