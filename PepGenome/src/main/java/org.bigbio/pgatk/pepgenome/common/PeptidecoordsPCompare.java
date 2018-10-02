package org.bigbio.pgatk.pepgenome.common;


import java.io.Serializable;
import java.util.Comparator;

//allows for the comparision of two peptide coordinate pointers.
//it dereferences them and calls the operator<
public class PeptidecoordsPCompare implements Comparator<PeptideCoordinates>, Serializable {

    private static final long serialVersionUID = -3542945510049682487L;

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