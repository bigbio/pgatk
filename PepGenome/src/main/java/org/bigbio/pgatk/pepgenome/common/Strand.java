package org.bigbio.pgatk.pepgenome.common;

import java.util.HashMap;
import java.util.Map;

//possible strands.
public enum Strand {
    fwd(1),
    rev(-1),
    unk(0);

    private int intValue;
    private static Map<Integer, Strand> mappings;

    private static Map<Integer, Strand> getMappings() {
        if (mappings == null) {
            synchronized (Strand.class) {
                if (mappings == null) {
                    mappings = new HashMap<>();
                }
            }
        }
        return mappings;
    }

    Strand(int value) {
        intValue = value;
        getMappings().put(value, this);
    }

    public int getValue() {
        return intValue;
    }

    public static Strand forValue(int value) {
        return getMappings().get(value);
    }
}