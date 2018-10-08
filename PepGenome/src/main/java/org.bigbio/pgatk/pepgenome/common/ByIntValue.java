package org.bigbio.pgatk.pepgenome.common;

import javafx.util.Pair;

import java.io.Serializable;
import java.util.Comparator;

public class ByIntValue implements Comparator<Pair<String, Integer>>, Serializable {

    private static final long serialVersionUID = -6930072683573085416L;

    @Override
    public int compare(Pair<String, Integer> o1, Pair<String, Integer> o2) {
        if (o1.getValue() < o2.getValue()) {
            return -1;
        }
        if (o2.getValue() < o1.getValue()) {
            return 1;
        }
        return 0;
    }
}
