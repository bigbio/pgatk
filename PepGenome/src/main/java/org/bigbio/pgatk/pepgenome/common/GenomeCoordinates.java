package org.bigbio.pgatk.pepgenome.common;

import java.io.Serializable;

//extension of coordinates, holds genomic coordinates.
public class GenomeCoordinates extends Coordinates<GenomeCoordinates> implements Serializable {

    private static final long serialVersionUID = -857165864692238920L;
    private String transcriptid;
    private String exonid;

    //holds the chromosome.
    private Chromosome chr;
    //holds the scaffolding.
    private String chrscaf;
    //holds the strand.
    private Strand strand;
    //holds the frame.
    private Frame frame;

    public GenomeCoordinates(GenomeCoordinates obj) {
        super(obj);
        this.transcriptid = obj.transcriptid;
        this.exonid = obj.exonid;
        this.chr = obj.chr;
        this.chrscaf = obj.chrscaf;
        this.strand = obj.strand;
        this.frame = obj.frame;
    }

    public GenomeCoordinates() {
        super();
    }

    @Override
    public int compareTo(GenomeCoordinates o) {
        if (lessThan(o)) {
            return -1;
        } else if (o.lessThan(this)) {
            return 1;
        }
        return 0;
    }

    @Override
    public boolean equals(Object obj) {
        return compareTo((GenomeCoordinates) obj) == 0;
    }

    private boolean comapre2(GenomeCoordinates lhs, GenomeCoordinates rhs) {
        if (lhs.chr.isScaffold() && rhs.chr.isScaffold() && lhs.chrscaf.equals(rhs.chrscaf)) {
            return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end >= rhs.start;
        }
        if (lhs.chr.isScaffold() && rhs.chr.isScaffold() && !lhs.chrscaf.equals(rhs.chrscaf)) {
            return lhs.chrscaf.compareTo(rhs.chrscaf) < 0;
        }
        if (lhs.chr.isScaffold() && !rhs.chr.isScaffold()) {
            return false;
        }
        if (!lhs.chr.isScaffold() && rhs.chr.isScaffold()) {
            return true;
        }
        if (lhs.chr == rhs.chr) {
            return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end >= rhs.start;
        }
        return lhs.chr.getValue() < rhs.chr.getValue();
    }

    public boolean equals2(GenomeCoordinates rhs) {
        return ((!chr.isScaffold() && chr == rhs.chr) || (chr.isScaffold() && chrscaf.equals(rhs.chrscaf))) && start >= rhs.start && end <= rhs.end;
    }

    private boolean lessThan(GenomeCoordinates rhs) {
        if (chr.isScaffold() && rhs.chr.isScaffold() && chrscaf.equals(rhs.chrscaf)) {
            if (start == rhs.start) {
                return end < rhs.end;
            }
            return start < rhs.start;
        }
        if (chr.isScaffold() && rhs.chr.isScaffold() && !chrscaf.equals(rhs.chrscaf)) {
            return chrscaf.compareTo(rhs.chrscaf) < 0;
        }
        if (chr.isScaffold() && !rhs.chr.isScaffold()) {
            return false;
        }
        if (!chr.isScaffold() && rhs.chr.isScaffold()) {
            return true;
        }
        if (chr == rhs.chr) {
            if (start == rhs.start) {
                return end < rhs.end;
            }
            return start < rhs.start;
        }
        return chr.getValue() < rhs.chr.getValue();
    }

    public String getTranscriptid() {
        return transcriptid;
    }

    public void setTranscriptid(String transcriptid) {
        this.transcriptid = transcriptid;
    }

    public String getExonid() {
        return exonid;
    }

    public void setExonid(String exonid) {
        this.exonid = exonid;
    }

    public Chromosome getChr() {
        return chr;
    }

    public void setChr(Chromosome chr) {
        this.chr = chr;
    }

    public String getChrscaf() {
        return chrscaf;
    }

    public void setChrscaf(String chrscaf) {
        this.chrscaf = chrscaf;
    }

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    public Frame getFrame() {
        return frame;
    }

    public void setFrame(Frame frame) {
        this.frame = frame;
    }

    @Override
    public String toString() {
        return "GenomeCoordinates{" +
                "transcriptid='" + transcriptid + '\'' +
                ", exonid='" + exonid + '\'' +
                ", chr=" + chr +
                ", chrscaf='" + chrscaf + '\'' +
                ", strand=" + strand +
                ", frame=" + frame +
                ", start=" + start +
                ", end=" + end +
                ", Nterm=" + Nterm +
                ", Cterm=" + Cterm +
                '}';
    }
}