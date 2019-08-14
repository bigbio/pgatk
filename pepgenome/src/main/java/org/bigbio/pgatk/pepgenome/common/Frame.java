package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

//possible frames
public enum Frame implements Serializable {
    frame1(1),
    frame2(2),
    frame3(3),
    unknown(0);

    private int intValue;
    private static Map<Integer, Frame> mappings;

    private static Map<Integer, Frame> getMappings() {
        if (mappings == null) {
            synchronized (Frame.class) {
                if (mappings == null) {
                    mappings = new HashMap<>();
                }
            }
        }
        return mappings;
    }

    Frame(int value) {
        intValue = value;
        getMappings().put(value, this);
    }

    public int getValue() {
        return intValue;
    }

    public static Frame forValue(int value) {
        return getMappings().get(value);
    }
}