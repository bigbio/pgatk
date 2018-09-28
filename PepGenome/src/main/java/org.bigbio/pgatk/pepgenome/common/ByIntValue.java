package org.bigbio.pgatk.pepgenome.common;

import javafx.util.Pair;

import java.util.Comparator;

public class ByIntValue implements Comparator<Pair<String, Integer>> {

    @Override
    public int compare(Pair<String, Integer> lhs, Pair<String, Integer> rhs) {
        if (lhs.getValue() < rhs.getValue()) {
            return -1;
        }
        if (rhs.getValue() < lhs.getValue()) {
            return 1;
        }
        return 0;
    }
}
