package org.bigbio.pgatk.pepgenome;

import java.util.Comparator;

//this comparator lets you compare PeptideEntry pointers.
//it just calls the operator < for the dereferenced pointers.
public class PeptideentryPcompare implements Comparator<PeptideEntry> {

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
