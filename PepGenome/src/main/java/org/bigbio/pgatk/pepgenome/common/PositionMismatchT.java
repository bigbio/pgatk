package org.bigbio.pgatk.pepgenome.common;

//this is made for a maximum of 2 mismatches but can easily be extended 
//(for example just replace the ints with a vaector of int.)
//care: this change has to be applied to KmerTreeMap::insert_into_gene_id_map

//this struct holds information about the position of a peptide in a protein.
//it also knows how many and where mismatches occured in the matching process.
public class PositionMismatchT {
    //position of the peptide in the proteinsequence (zerobased)
    private int m_position_in_protein;
    //first mismatch position in the proteinsequence (zerobased)
    private int m_first_mismatch_positon;
    //second mismatch position in the proteinsequence (zerobased)
    private int m_second_mismatch_positon;

    //ctr
    public PositionMismatchT(int posInProtein, int firstMismatch, int secondMismatch) {
        this.m_position_in_protein = posInProtein;
        this.m_first_mismatch_positon = firstMismatch;
        this.m_second_mismatch_positon = secondMismatch;
    }

    //returns the positon of the peptide in the protein or -1 if the peptide does not match.
    public int position_in_protein() {
        return m_position_in_protein;
    }

    //returns the position of the first mismatch or -1 if there was no mismatch.
    public int first() {
        return m_first_mismatch_positon;
    }

    //returns the position of the second mismatch or -1 if there was no mismatch.
    public int second() {
        return m_second_mismatch_positon;
    }
}