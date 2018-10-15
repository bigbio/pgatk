package org.bigbio.pgatk.pepgenome.common;

import java.util.HashMap;
import java.util.Map;

public enum Assembly {
    none(0),
    primary(1),
    patchhaploscaff(2);

    private int intValue;
    private static Map<Integer, Assembly> mappings;

    private static Map<Integer, Assembly> getMappings() {
        if (mappings == null) {
            synchronized (Assembly.class) {
                if (mappings == null) {
                    mappings = new HashMap<>();
                }
            }
        }
        return mappings;
    }

    Assembly(int value) {
        intValue = value;
        getMappings().put(value, this);
    }

    public int getValue() {
        return intValue;
    }

    public static Assembly forValue(int value) {
        return getMappings().get(value);
    }
}