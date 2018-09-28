package org.bigbio.pgatk.pepgenome.common;

public class TaxonomyIdentifiers {
    private String geneId;
    private String transcriptId;
    private String exonId;
    private int length;

    public TaxonomyIdentifiers(String geneId, String transcriptId, String exonId, int length) {
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.exonId = exonId;
        this.length = length;
    }

    public static final TaxonomyIdentifiers cow = new TaxonomyIdentifiers("ENSBTAG", "ENSBTAT", "ENSBTAE", 18);
    public static final TaxonomyIdentifiers marmoset = new TaxonomyIdentifiers("ENSCJAG", "ENSCJAT", "ENSCJAE", 18);
    public static final TaxonomyIdentifiers dog = new TaxonomyIdentifiers("ENSCAFG", "ENSCAFT", "ENSCFAE", 18);
    public static final TaxonomyIdentifiers vervet_AGM = new TaxonomyIdentifiers("ENSCSAG", "ENSCSAT", "ENSCSAE", 18);
    public static final TaxonomyIdentifiers cintestinalis = new TaxonomyIdentifiers("ENSCING", "ENSCINT", "ENSCINE", 18);
    public static final TaxonomyIdentifiers horse = new TaxonomyIdentifiers("ENSECAG", "ENSECAT", "ENSECAE", 18);
    public static final TaxonomyIdentifiers cat = new TaxonomyIdentifiers("ENSFCAG", "ENSFCAT", "ENSFCAE", 18);
    public static final TaxonomyIdentifiers chicken = new TaxonomyIdentifiers("ENSGALG", "ENSGALT", "ENSGALE", 18);
    public static final TaxonomyIdentifiers gorilla = new TaxonomyIdentifiers("ENSGGOG", "ENSGGOT", "ENSGGOE", 18);
    public static final TaxonomyIdentifiers human = new TaxonomyIdentifiers("ENSG", "ENST", "ENSE", 15);
    public static final TaxonomyIdentifiers macaque = new TaxonomyIdentifiers("ENSMMUG", "ENSMMUT", "ENSMMUE", 18);
    public static final TaxonomyIdentifiers turkey = new TaxonomyIdentifiers("ENSMGAG", "ENSMGAT", "ENSMGAE", 18);
    public static final TaxonomyIdentifiers opossum = new TaxonomyIdentifiers("ENSMODG", "ENSMODT", "ENSMODE", 18);
    public static final TaxonomyIdentifiers mouse = new TaxonomyIdentifiers("ENSMUSG", "ENSMUST", "ENSMUSE", 18);
    public static final TaxonomyIdentifiers platypus = new TaxonomyIdentifiers("ENSOANG", "ENSOANT", "ENSOANE", 18);
    public static final TaxonomyIdentifiers rabbit = new TaxonomyIdentifiers("ENSOCUG", "ENSOCUT", "ENSOCUE", 18);
    public static final TaxonomyIdentifiers medaka = new TaxonomyIdentifiers("ENSORLG", "ENSORLT", "ENSORLE", 18);
    public static final TaxonomyIdentifiers sheep = new TaxonomyIdentifiers("ENSOARG", "ENSOART", "ENSOARE", 18);
    public static final TaxonomyIdentifiers chimp = new TaxonomyIdentifiers("ENSPTRG", "ENSPTRT", "ENSPTRE", 18);
    public static final TaxonomyIdentifiers olivebaboon = new TaxonomyIdentifiers("ENSPANG", "ENSPANT", "ENSPANE", 18);
    public static final TaxonomyIdentifiers orangutan = new TaxonomyIdentifiers("ENSPPYG", "ENSPPYT", "ENSPPYE", 18);
    public static final TaxonomyIdentifiers rat = new TaxonomyIdentifiers("ENSRNOG", "ENSRNOT", "ENSRNOE", 18);
    public static final TaxonomyIdentifiers pig = new TaxonomyIdentifiers("ENSSSCG", "ENSSSCT", "ENSSSCE", 18);
    public static final TaxonomyIdentifiers zebrafinch = new TaxonomyIdentifiers("ENSTGUG", "ENSTGUT", "ENSTGUE", 18);
    public static final TaxonomyIdentifiers tetraodon = new TaxonomyIdentifiers("ENSTNIG", "ENSTNIT", "ENSTNIE", 18);

    public String getGeneId() {
        return geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public String getExonId() {
        return exonId;
    }

    public int getLength() {
        return length;
    }
}