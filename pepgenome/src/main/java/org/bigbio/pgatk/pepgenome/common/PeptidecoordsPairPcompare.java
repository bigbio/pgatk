package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;
import java.util.Comparator;

//allows for the comparision of std::pair<common.PeptideCoordinates*, common.GenomeCoordinates> pointers.
//dereferences and compares the std::pair::first (common.PeptideCoordinates::operator<)
public class PeptidecoordsPairPcompare implements Comparator<Tuple<PeptideCoordinates, GenomeCoordinates>>, Serializable {

    private static final long serialVersionUID = -7691292743574950228L;

    @Override
    public int compare(Tuple<PeptideCoordinates, GenomeCoordinates> o1, Tuple<PeptideCoordinates, GenomeCoordinates> o2) {
        if (o1.getKey().lessThan(o2.getKey())) {
            return -1;
        }
        if (o2.getKey().lessThan(o1.getKey())) {
            return 1;
        }
        return 0;
    }
}