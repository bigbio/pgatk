package org.bigbio.pgatk.pepgenome;

import java.util.Comparator;

//comparator for MapEntry pointers. dereferences the pointers and calls the operator< method.
public class MapentryPCompare implements Comparator<MapEntry> {

    @Override
    public int compare(MapEntry lhs, MapEntry rhs) {
        if (lhs.lessThan(rhs)) {
            return -1;
        }
        if (rhs.lessThan(lhs)) {
            return 1;
        }
        return 0;
    }
}