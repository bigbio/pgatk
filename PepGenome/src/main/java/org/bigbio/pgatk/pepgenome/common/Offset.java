package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

//possible offsets.
public enum Offset implements Serializable {
    off1(1),
    off2(2),
    off3(3);

    private int intValue;
    private static Map<Integer, Offset> mappings;

    private static Map<Integer, Offset> getMappings() {
        if (mappings == null) {
            synchronized (Offset.class) {
                if (mappings == null) {
                    mappings = new HashMap<>();
                }
            }
        }
        return mappings;
    }

    Offset(int value) {
        intValue = value;
        getMappings().put(value, this);
    }

    public int getValue() {
        return intValue;
    }

    public static Offset forValue(int value) {
        return getMappings().get(value);
    }
}