package bigbio.pgatk.jpogo.common;

import javafx.util.Pair;

import java.util.Comparator;

//allows for the comparision of std::pair<common.PeptideCoordinates*, common.GenomeCoordinates> pointers.
//dereferences and compares the std::pair::first (common.PeptideCoordinates::operator<)
public class PeptidecoordsPairPcompare implements Comparator<Pair<PeptideCoordinates, GenomeCoordinates>> {

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