package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

//possible chromosomes
public class Chromosome implements Serializable {

	private static final long serialVersionUID = -3608154103569444727L;
	
	private static Map<String, Integer> chrToInt;
    private static Map<Integer, String> intToChr;
    private static Set<String> scaffoldnames;

    private String name;
    
    public Chromosome(String name) {
    	if(getChrToInt().containsKey(name) || getScaffoldNames().contains(name)) {
    		this.name = name;
    	} else {
    		this.name = "NA";
    	}
    }

    private static Map<String, Integer> getChrToInt() {
        if (chrToInt == null) {
            synchronized (Chromosome.class) {
                if (chrToInt == null) {
                	chrToInt = new HashMap<>();
                }
            }
        }
        return chrToInt;
    }
    
    private static Map<Integer, String> getIntToChr() {
        if (intToChr == null) {
            synchronized (Chromosome.class) {
                if (intToChr == null) {
                	intToChr = new HashMap<>();
                }
            }
        }
        return intToChr;
    }
    
    private static Set<String> getScaffoldNames() {
    	if (scaffoldnames == null) {
    		synchronized (Chromosome.class) {
    			if (scaffoldnames == null) {
    				scaffoldnames = new HashSet<>();
    			}
    		}
    	}
    	return scaffoldnames;
    }


    public int getValue() {
    	if(getChrToInt().containsKey(name)) {
    		return getChrToInt().get(name);
    	} else {
    		return -1;
    	}
    }
    
    public String getName() {
    	return name;
    }
    
    public boolean isScaffold() {
    	return getScaffoldNames().contains(name);
    }
    
    public boolean isNA() {
    	return name.equals("NA");
    }

    public static String forValue(int value) {
    	if(getIntToChr().containsKey(value)) {
    		return getIntToChr().get(value);
    	} else {
    		return "NA";
    	}
    }
    
    public static int forName(String name) {
    	if(getChrToInt().containsKey(name)) {
    		return getChrToInt().get(name);
    	} else {
    		return -1;
    	}
    }
    
    public static void addChr(String name) {
    	String tmpname = name;
    	if(tmpname.startsWith("chr") || tmpname.startsWith("Chr")) {
    		tmpname = tmpname.substring(3);
    	}
    	if(!getChrToInt().containsKey(tmpname)) {
    		getChrToInt().put(tmpname, getChrToInt().size()+1);
    		getIntToChr().put(getChrToInt().get(tmpname), tmpname);
    		if(tmpname.equals("M")) {
    			getChrToInt().put("MT", getChrToInt().size()+1);
        		getIntToChr().put(getChrToInt().get("MT"), "MT");
    		} else if(tmpname.equals("MT")) {
    			getChrToInt().put("M", getChrToInt().size()+1);
        		getIntToChr().put(getChrToInt().get("M"), "M");
    		}
    	}
    }
    
    public static void addScaffold(String name) {
    	if(!getScaffoldNames().contains(name)) {
    		getScaffoldNames().add(name);
    	}
    }
}