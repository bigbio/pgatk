package org.bigbio.pgatk.pepgenome.common;


import java.util.Comparator;

//allows for the comparision of two peptide coordinate pointers.
//it dereferences them and calls the operator<
public class PeptidecoordsPCompare implements Comparator<PeptideCoordinates> {

    @Override
    public int compare(PeptideCoordinates lhs, PeptideCoordinates rhs) {
        if (lhs.lessThan(rhs)) {
            return -1;
        }
        if (rhs.lessThan(lhs)) {
            return 1;
        }
        return 0;
    }
}