package org.bigbio.pgatk.pepgenome.common.comparators;

import org.bigbio.pgatk.pepgenome.common.PeptideEntry;

import java.io.Serializable;
import java.util.Comparator;

/**
 * This comparator lets you compare PeptideEntry pointers.
 * it just calls the operator < for the dereferenced pointers.
 *
 * @author ypriverol
 */

public class PeptideentryPcompare implements Comparator<PeptideEntry>, Serializable {

    private static final long serialVersionUID = -8833743639265606264L;

    @Override
    public int compare(PeptideEntry o1, PeptideEntry o2) {
        if (o1.lessThan(o2)) {
            return -1;
        }
        if (o2.lessThan(o1)) {
            return 1;
        }
        return 0;
    }
}
