package bigbio.pgatk.jpogo.common;

import java.util.HashMap;
import java.util.Map;

//possible chromosomes
public enum Chromosome {
    chr1(1),
    chr1A(2),
    chr1B(3),
    chr2(4),
    chr2A(5),
    chr2a(6),
    chr2B(7),
    chr2b(8),
    chr3(9),
    chr4(10),
    chr4A(11),
    chr5(12),
    chr6(13),
    chr7(14),
    chr8(15),
    chr9(16),
    chr10(17),
    chr11(18),
    chr12(19),
    chr13(20),
    chr14(21),
    chr15(22),
    chr16(23),
    chr17(24),
    chr18(25),
    chr19(26),
    chr20(27),
    chr21(28),
    chr22(29),
    chr23(30),
    chr24(31),
    chr25(32),
    chr26(33),
    chr27(34),
    chr28(35),
    chr29(36),
    chr30(37),
    chr31(38),
    chr32(39),
    chr33(40),
    chr34(41),
    chr35(42),
    chr36(43),
    chr37(44),
    chr38(45),
    chr39(46),
    chr40(47),
    chrX(48),
    chrY(49),
    chrXY(50),
    chrX1(51),
    chrX2(52),
    chrX3(53),
    chrX5(54),
    chrA1(55),
    chrA2(56),
    chrA3(57),
    chrB1(58),
    chrB2(59),
    chrB3(60),
    chrB4(61),
    chrC1(62),
    chrC2(63),
    chrD1(64),
    chrD2(65),
    chrD3(66),
    chrD4(67),
    chrE1(68),
    chrE2(69),
    chrE3(70),
    chrF1(71),
    chrF2(72),
    chrLG2(73),
    chrLG5(74),
    chrLGE22(75),
    chrW(76),
    chrZ(77),
    chrM(78),
    chrLGE64(79),
    chrNA(-1),
    scaffold(0);

    private int intValue;
    private static Map<Integer, Chromosome> mappings;

    private static Map<Integer, Chromosome> getMappings() {
        if (mappings == null) {
            synchronized (Chromosome.class) {
                if (mappings == null) {
                    mappings = new HashMap<>();
                }
            }
        }
        return mappings;
    }

    Chromosome(int value) {
        intValue = value;
        getMappings().put(value, this);
    }

    public int getValue() {
        return intValue;
    }

    public static Chromosome forValue(int value) {
        return getMappings().get(value);
    }
}