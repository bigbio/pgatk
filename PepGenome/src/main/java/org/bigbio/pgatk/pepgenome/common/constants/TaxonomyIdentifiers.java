package org.bigbio.pgatk.pepgenome.common.constants;

/**
 * Taxonomy identifiers for Gene, Transcript, Exon, length of the Id
 *
 * @author ypriverol
 */
public enum TaxonomyIdentifiers {

    COW("ENSBTAG", "ENSBTAT", "ENSBTAE", 18),
    MARMOSET("ENSCJAG", "ENSCJAT", "ENSCJAE", 18),
    DOG("ENSCAFG", "ENSCAFT", "ENSCFAE", 18),
    VERVERT_AGM("ENSCSAG", "ENSCSAT", "ENSCSAE", 18),
    CITESTINALIS("ENSCING", "ENSCINT", "ENSCINE", 18),
    HORSE("ENSECAG", "ENSECAT", "ENSECAE", 18),
    CAT("ENSFCAG", "ENSFCAT", "ENSFCAE", 18),
    CHICKEN("ENSGALG", "ENSGALT", "ENSGALE", 18),
    GORILLA("ENSGGOG", "ENSGGOT", "ENSGGOE", 18),
    HUMAN("ENSG", "ENST", "ENSE", 15),
    MACAQUE("ENSMMUG", "ENSMMUT", "ENSMMUE", 18),
    TURKEY("ENSMGAG", "ENSMGAT", "ENSMGAE", 18),
    OPOSUM("ENSMODG", "ENSMODT", "ENSMODE", 18),
    MOUSE("ENSMUSG", "ENSMUST", "ENSMUSE", 18),
    PLATYPUS("ENSOANG", "ENSOANT", "ENSOANE", 18),
    RABBIT("ENSOCUG", "ENSOCUT", "ENSOCUE", 18),
    MEDAKA("ENSORLG", "ENSORLT", "ENSORLE", 18),
    SHEEP("ENSOARG", "ENSOART", "ENSOARE", 18),
    CHIMP("ENSPTRG", "ENSPTRT", "ENSPTRE", 18),
    OLIVEBABOOM("ENSPANG", "ENSPANT", "ENSPANE", 18),
    ORANGUTAN("ENSPPYG", "ENSPPYT", "ENSPPYE", 18),
    RAT("ENSRNOG", "ENSRNOT", "ENSRNOE", 18),
    PIG("ENSSSCG", "ENSSSCT", "ENSSSCE", 18),
    ZEBRAFINCH("ENSTGUG", "ENSTGUT", "ENSTGUE", 18),
    ZEBRAFISH("ENSDARG", "ENSDART", "ENSDARE", 18),
    TETRAODON("ENSTNIG", "ENSTNIT", "ENSTNIE", 18),
	YEAST("Y","Y","Y",8);

    private String geneId;
    private String transcriptId;
    private String exonId;
    private int length;

    TaxonomyIdentifiers(String geneId, String transcriptId, String exonId, int length) {
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.exonId = exonId;
        this.length = length;
    }

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