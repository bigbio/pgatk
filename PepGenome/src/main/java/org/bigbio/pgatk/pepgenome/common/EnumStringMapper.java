package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;

import java.util.HashMap;
import java.util.Map;

import static org.bigbio.pgatk.pepgenome.common.Chromosome.*;
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
        switch (chr) {
            case chr1:
                return "1";
            case chr1A:
                return "1A";
            case chr1B:
                return "1B";
            case chr2:
                return "2";
            case chr2A:
                return "2A";
            case chr2a:
                return "2a";
            case chr2B:
                return "2B";
            case chr2b:
                return "2b";
            case chr3:
                return "3";
            case chr4:
                return "4";
            case chr4A:
                return "4A";
            case chr5:
                return "5";
            case chr6:
                return "6";
            case chr7:
                return "7";
            case chr8:
                return "8";
            case chr9:
                return "9";
            case chr10:
                return "10";
            case chr11:
                return "11";
            case chr12:
                return "12";
            case chr13:
                return "13";
            case chr14:
                return "14";
            case chr15:
                return "15";
            case chr16:
                return "16";
            case chr17:
                return "17";
            case chr18:
                return "18";
            case chr19:
                return "19";
            case chr20:
                return "20";
            case chr21:
                return "21";
            case chr22:
                return "22";
            case chr23:
                return "23";
            case chr24:
                return "24";
            case chr25:
                return "25";
            case chr26:
                return "26";
            case chr27:
                return "27";
            case chr28:
                return "28";
            case chr29:
                return "29";
            case chr30:
                return "30";
            case chr31:
                return "31";
            case chr32:
                return "32";
            case chr33:
                return "33";
            case chr34:
                return "34";
            case chr35:
                return "35";
            case chr36:
                return "36";
            case chr37:
                return "37";
            case chr38:
                return "38";
            case chr39:
                return "39";
            case chr40:
                return "40";
            case chrI:
            	return "I";
            case chrII:
            	return "II";
            case chrIII:
            	return "III";
            case chrIV:
            	return "IV";
            case chrV:
            	return "V";
            case chrVI:
            	return "VI";
            case chrVII:
            	return "VII";
            case chrVIII:
            	return "VIII";
            case chrIX:
            	return "IX";
            case chrX:
                return "X";
            case chrXI:
            	return "XI";
            case chrXII:
            	return "XII";
            case chrXIII:
            	return "XIII";
            case chrXIV:
            	return "XIV";
            case chrXV:
            	return "XV";
            case chrXVI:
            	return "XVI";
            case chrY:
                return "Y";
            case chrXY:
                return "XY";
            case chrX1:
                return "X1";
            case chrX2:
                return "X2";
            case chrX3:
                return "X3";
            case chrX5:
                return "X5";
            case chrA1:
                return "A1";
            case chrA2:
                return "A2";
            case chrA3:
                return "A3";
            case chrB1:
                return "B1";
            case chrB2:
                return "B2";
            case chrB3:
                return "B3";
            case chrB4:
                return "B4";
            case chrC1:
                return "C1";
            case chrC2:
                return "C2";
            case chrD1:
                return "D1";
            case chrD2:
                return "D2";
            case chrD3:
                return "D3";
            case chrD4:
                return "D4";
            case chrE1:
                return "E1";
            case chrE2:
                return "E2";
            case chrE3:
                return "E3";
            case chrF1:
                return "F1";
            case chrF2:
                return "F2";
            case chrLG2:
                return "LG2";
            case chrLG5:
                return "LG5";
            case chrLGE22:
                return "LGE22";
            case chrW:
                return "W";
            case chrZ:
                return "Z";
            case chrM:
                return "MT";
            case chrMito:
            	return "Mito";
            case scaffold:
                return "";
            default:
                return "unknown";
        }
    }

    public static String enumToChrString(Chromosome chr) {
        switch (chr) {
            case chr1:
                return "chr1";
            case chr1A:
                return "chr1A";
            case chr1B:
                return "chr1B";
            case chr2:
                return "chr2";
            case chr2A:
                return "chr2A";
            case chr2a:
                return "chr2a";
            case chr2B:
                return "chr2B";
            case chr2b:
                return "chr2b";
            case chr3:
                return "chr3";
            case chr4:
                return "chr4";
            case chr4A:
                return "chr4A";
            case chr5:
                return "chr5";
            case chr6:
                return "chr6";
            case chr7:
                return "chr7";
            case chr8:
                return "chr8";
            case chr9:
                return "chr9";
            case chr10:
                return "chr10";
            case chr11:
                return "chr11";
            case chr12:
                return "chr12";
            case chr13:
                return "chr13";
            case chr14:
                return "chr14";
            case chr15:
                return "chr15";
            case chr16:
                return "chr16";
            case chr17:
                return "chr17";
            case chr18:
                return "chr18";
            case chr19:
                return "chr19";
            case chr20:
                return "chr20";
            case chr21:
                return "chr21";
            case chr22:
                return "chr22";
            case chr23:
                return "chr23";
            case chr24:
                return "chr24";
            case chr25:
                return "chr25";
            case chr26:
                return "chr26";
            case chr27:
                return "chr27";
            case chr28:
                return "chr28";
            case chr29:
                return "chr29";
            case chr30:
                return "chr30";
            case chr31:
                return "chr31";
            case chr32:
                return "chr32";
            case chr33:
                return "chr33";
            case chr34:
                return "chr34";
            case chr35:
                return "chr35";
            case chr36:
                return "chr36";
            case chr37:
                return "chr37";
            case chr38:
                return "chr38";
            case chr39:
                return "chr39";
            case chr40:
                return "chr40";
            case chrI:
            	return "chrI";
            case chrII:
            	return "chrII";
            case chrIII:
            	return "chrIII";
            case chrIV:
            	return "chrIV";
            case chrV:
            	return "chrV";
            case chrVI:
            	return "chrVI";
            case chrVII:
            	return "chrVII";
            case chrVIII:
            	return "chrVIII";
            case chrIX:
            	return "chrIX";
            case chrX:
                return "chrX";
            case chrXI:
            	return "chrXI";
            case chrXII:
            	return "chrXII";
            case chrXIII:
            	return "chrXIII";
            case chrXIV:
            	return "chrXIV";
            case chrXV:
            	return "chrXV";
            case chrXVI:
            	return "chrXVI";
            case chrY:
                return "chrY";
            case chrXY:
                return "chrXY";
            case chrX1:
                return "chrX1";
            case chrX2:
                return "chrX2";
            case chrX3:
                return "chrX3";
            case chrX5:
                return "chrX5";
            case chrA1:
                return "chrA1";
            case chrA2:
                return "chrA2";
            case chrA3:
                return "chrA3";
            case chrB1:
                return "chrB1";
            case chrB2:
                return "chrB2";
            case chrB3:
                return "chrB3";
            case chrB4:
                return "chrB4";
            case chrC1:
                return "chrC1";
            case chrC2:
                return "chrC2";
            case chrD1:
                return "chrD1";
            case chrD2:
                return "chrD2";
            case chrD3:
                return "chrD3";
            case chrD4:
                return "chrD4";
            case chrE1:
                return "chrE1";
            case chrE2:
                return "chrE2";
            case chrE3:
                return "chrE3";
            case chrF1:
                return "chrF1";
            case chrF2:
                return "chrF2";
            case chrLG2:
                return "chrLG2";
            case chrLG5:
                return "chrLG5";
            case chrLGE22:
                return "chrLGE22";
            case chrW:
                return "chrW";
            case chrZ:
                return "chrZ";
            case chrM:
                return "chrM";
            case chrMito:
            	return "chrMito";
            case scaffold:
                return "scaffold";
            default:
                return "unknown";
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

        if (substr.equals("I")) {
            return chrI;
        }
        if (substr.equals("II")) {
            return chrII;
        }
        if (substr.equals("III")) {
            return chrIII;
        }
        if (substr.equals("IV")) {
            return chrIV;
        }
        if (substr.equals("V")) {
            return chrV;
        }
        if (substr.equals("VI")) {
            return chrVI;
        }
        if (substr.equals("VII")) {
            return chrVII;
        }
        if (substr.equals("VIII")) {
            return chrVIII;
        }
        if (substr.equals("IX")) {
            return chrIX;
        }
        if (substr.equals("X")) {
            return chrX;
        }
        if (substr.equals("XI")) {
            return chrXI;
        }
        if (substr.equals("XII")) {
            return chrXII;
        }
        if (substr.equals("XIII")) {
            return chrXIII;
        }
        if (substr.equals("XIV")) {
            return chrXIV;
        }
        if (substr.equals("XV")) {
            return chrXV;
        }
        if (substr.equals("XVI")) {
            return chrXVI;
        }
        if (substr.equals("W")) {
            return chrW;
        }
        if (substr.equals("Z")) {
            return chrZ;
        }
        if (substr.equals("A1")) {
            return chrA1;
        }
        if (substr.equals("A2")) {
            return chrA2;
        }
        if (substr.equals("A3")) {
            return chrA3;
        }
        if (substr.equals("B1")) {
            return chrB1;
        }
        if (substr.equals("B2")) {
            return chrB2;
        }
        if (substr.equals("B3")) {
            return chrB3;
        }
        if (substr.equals("B4")) {
            return chrB4;
        }
        if (substr.equals("C1")) {
            return chrC1;
        }
        if (substr.equals("C2")) {
            return chrC2;
        }
        if (substr.equals("D1")) {
            return chrD1;
        }
        if (substr.equals("D2")) {
            return chrD2;
        }
        if (substr.equals("D3")) {
            return chrD3;
        }
        if (substr.equals("D4")) {
            return chrD4;
        }
        if (substr.equals("E1")) {
            return chrE1;
        }
        if (substr.equals("E2")) {
            return chrE2;
        }
        if (substr.equals("E3")) {
            return chrE3;
        }
        if (substr.equals("F1")) {
            return chrF1;
        }
        if (substr.equals("F2")) {
            return chrF2;
        }
        if (substr.equals("LGE64")) {
            return chrLGE64;
        }
        if (substr.equals("2A")) {
            return chr2A;
        }
        if (substr.equals("2B")) {
            return chr2B;
        }
        if (substr.equals("2a")) {
            return chr2a;
        }
        if (substr.equals("2b")) {
            return chr2b;
        }
        if (substr.equals("X1")) {
            return chrX1;
        }
        if (substr.equals("X2")) {
            return chrX2;
        }
        if (substr.equals("X3")) {
            return chrX3;
        }
        if (substr.equals("X5")) {
            return chrX5;
        }
        if (substr.equals("1A")) {
            return chr1A;
        }
        if (substr.equals("1B")) {
            return chr1B;
        }
        if (substr.equals("4A")) {
            return chr4A;
        }
        if (substr.equals("LG2")) {
            return chrLG2;
        }
        if (substr.equals("LG5")) {
            return chrLG5;
        }
        if (substr.equals("LGE22")) {
            return chrLGE22;
        }
        if (substr.equals("Y")) {
            return chrY;
        }
        if (substr.equals("XY")) {
            return chrXY;
        }
        if (substr.equals("M") || substr.equals("MT")) {
            return chrM;
        }
        if (substr.equals("Mito")) {
        	return chrMito;
        }
        if (substr.length() >= 2) {
            Chromosome scaffoldResult = GenomeMapper.ScaffoldIdentifier.findScaffoldTwoCode(substr);
            if(scaffoldResult != null)
                return scaffoldResult;
        }
        if (substr.length() >= 4) {
            Chromosome scaffoldResult = GenomeMapper.ScaffoldIdentifier.findScaffoldThreeCode(substr);
            if(scaffoldResult != null)
                return scaffoldResult;
        }
        int chromNum = Integer.parseInt(substr);
        if (chromNum == 1) {
            chromNum = 1;
        } else if (chromNum == 2) {
            chromNum = 4;
        } else if (chromNum == 3 || chromNum == 4) {
            chromNum = chromNum + 6;
        } else {
            chromNum = chromNum + 7;
        }
        return Chromosome.forValue(chromNum);
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