package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;

import java.util.HashMap;
import java.util.Map;

import static org.bigbio.pgatk.pepgenome.common.Strand.*;

public class EnumStringMapper {
    private static Map<String, String> ptmToColours = new HashMap<String, String>() {

        private static final long serialVersionUID = 7898239362684593527L;

        {
        put("phospho", "255,51,51");
        put("acetyl", "204,102,0");
        put("amidated", "255,153,51");
        put("oxidation", "204,204,0");
        put("methyl", "0,204,0");
        put("glygly", "51,255,51");
        put("gg", "51,255,51");
        put("sulfo", "0,204,204");
        put("palmitoyl", "51,153,255");
        put("formyl", "0,0,204");
        put("deamidated", "51,51,255");
    }};

    public static String enumToString(Strand strand) {
        return enumToString(strand, true);
    }

    public static String enumToString(Strand strand, boolean numeric) {
        if (numeric) {
            switch (strand) {
                case fwd:
                    return "1";
                case rev:
                    return "-1";
                default:
                    return "0";
            }
        }
        switch (strand) {
            case fwd:
                return "+";
            case rev:
                return "-";
            default:
                return ".";
        }
    }

    public static String enumToString(Chromosome chr) {
    	String name = chr.getName();
    	if(name.startsWith("chr") || name.startsWith("Chr")) {
    		return name.substring(3);
    	} else {
    		return name;
    	}
    }

    public static String enumToChrString(Chromosome chr) {
        String name = enumToString(chr);
        if(!name.startsWith("chr") && !name.startsWith("Chr") && !chr.isScaffold() && !chr.isNA()) {
        	return "chr" + name;
        } else {
        	return name;
        }
        		
    }

    public static Strand string_to_strand(String str) {
        if (str.equals("-1") || str.equals("-")) {
            return rev;
        }
        if (str.equals("1") || str.equals("+")) {
            return fwd;
        }
        return unk;
    }

    public static Frame string_to_frame(String str) {
        if (str.equals("1")) {
            return Frame.frame1;
        }
        if (str.equals("2")) {
            return Frame.frame2;
        }
        if (str.equals("3") || str.equals("0")) {
            return Frame.frame3;
        }
        return Frame.unknown;
    }

    public static Chromosome string_to_chromosome(String str) {
        String substr;
        if (str.startsWith("chr") || str.startsWith("Chr")) {
            substr = str.substring(3);
        } else {
            substr = str;
        }
        if(GenomeMapper.PEPTIDE_MAPPER.CHR_FROM_GENOME_FASTA==false) {
        	Chromosome.addChr(substr);
        }
        return new Chromosome(substr);
    }

    public static String ptmToColour(String ptmPSIname) {
        String name = ptmPSIname.toLowerCase();
        String color = ptmToColours.get(name);
        if (color == null) {
            return "255,51,153";
        }
        return color;
    }


}