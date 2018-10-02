package org.bigbio.pgatk.pepgenome.common.comparators;

import org.bigbio.pgatk.pepgenome.common.MapEntry;

import java.io.Serializable;
import java.util.Comparator;

//comparator for MapEntry pointers. dereferences the pointers and calls the operator< method.
public class MapentryPCompare implements Comparator<MapEntry>, Serializable {

    private static final long serialVersionUID = 2616763124059915920L;

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