package org.bigbio.pgatk.pepgenome.common;

import java.util.*;

public class CustomMultimap<K, V> {

    private TreeMap<K, ArrayList<V>> map;

    private CustomMultimap() {
        map = new TreeMap<>();
    }

    public static <K, V> CustomMultimap<K, V> create() {
        return new CustomMultimap();
    }

    public void put(K key, V value) {
        ArrayList<V> vs = map.computeIfAbsent(key, j -> new ArrayList<>());
        vs.add(value);
    }

    public ArrayList<V> get(K key) {
        ArrayList<V> vs = map.get(key);
        if (vs == null) {
            vs = new ArrayList<>();
        }
        return vs;
    }

    public boolean isEmpty() {
        return map.isEmpty();
    }

    public Set<K> keySet() {
        return map.keySet();
    }

    public Collection<ArrayList<V>> values() {
        return map.values();
    }

    public Set<Map.Entry<K, ArrayList<V>>> entrySet() {
        return map.entrySet();
    }
}
