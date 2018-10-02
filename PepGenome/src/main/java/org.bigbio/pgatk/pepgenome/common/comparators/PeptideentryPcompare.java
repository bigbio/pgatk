package org.bigbio.pgatk.pepgenome.common.comparators;

import org.bigbio.pgatk.pepgenome.common.PeptideEntry;

import java.io.Serializable;
import java.util.Comparator;

//this comparator lets you compare PeptideEntry pointers.
//it just calls the operator < for the dereferenced pointers.
public class PeptideentryPcompare implements Comparator<PeptideEntry>, Serializable {

    private static final long serialVersionUID = -8833743639265606264L;

    @Override
    public int compare(PeptideEntry lhs, PeptideEntry rhs) {
        if (lhs.lessThan(rhs)) {
            return -1;
        }
        if (rhs.lessThan(lhs)) {
            return 1;
        }
        return 0;
    }
}
