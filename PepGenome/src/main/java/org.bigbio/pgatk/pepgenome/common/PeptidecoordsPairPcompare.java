package org.bigbio.pgatk.pepgenome.common;

import javafx.util.Pair;

import java.io.Serializable;
import java.util.Comparator;

//allows for the comparision of std::pair<common.PeptideCoordinates*, common.GenomeCoordinates> pointers.
//dereferences and compares the std::pair::first (common.PeptideCoordinates::operator<)
public class PeptidecoordsPairPcompare implements Comparator<Pair<PeptideCoordinates, GenomeCoordinates>>, Serializable {

    private static final long serialVersionUID = -7691292743574950228L;

    @Override
    public int compare(Pair<PeptideCoordinates, GenomeCoordinates> lhs, Pair<PeptideCoordinates, GenomeCoordinates> rhs) {
        if (lhs.getKey().lessThan(rhs.getKey())) {
            return -1;
        }
        if (rhs.getKey().lessThan(lhs.getKey())) {
            return 1;
        }
        return 0;
    }
}