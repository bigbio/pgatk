package org.bigbio.pgatk.pepgenome.common;

import lombok.Data;

import java.io.Serializable;

//holds coordinates.
@Data
public class Coordinates<T> implements Comparable<T>, Serializable {

    private static final long serialVersionUID = -2348800999526373164L;
    //the start position
    protected int start;
    //the end position
    protected int end;
    //nterm offset (see enum common.Offset)
    protected Offset Nterm;
    //cterm offset (see enum common.Offset)
    protected Offset Cterm;

    public Coordinates() {
    }

    // lt opeator. returns true if lhs.start is lesser than rhs.start and lhs.end is smaller than rhs.end and lhs.end is smaller than rhs.start
    //otherwise returns false.
    private boolean lt(Coordinates lhs, Coordinates rhs) {
        return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end < rhs.start;
    }

    //equal oparator. returns true if lhs.start is bigger than or equal to rhs.start and lhs.end smaller than or equal to rhs.end
    //otherwise returns false.
    public boolean equalsTo(Coordinates rhs) {
        return start >= rhs.start && end <= rhs.end;
    }

    @Override
    public int compareTo(T o) {
        if (lt(this, (Coordinates) o)) {
            return -1;
        } else if (lt((Coordinates) o, this)) {
            return 1;
        }
        return 0;
    }

    @Override
    public boolean equals(Object obj) {
        return compareTo((T) obj) == 0;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public Offset getNterm() {
        return Nterm;
    }

    public void setNterm(Offset nterm) {
        Nterm = nterm;
    }

    public Offset getCterm() {
        return Cterm;
    }

    public void setCterm(Offset cterm) {
        Cterm = cterm;
    }

    public Coordinates(Coordinates obj) {
        this.start = obj.start;
        this.end = obj.end;
        Nterm = obj.Nterm;
        Cterm = obj.Cterm;
    }

    @Override
    public String toString() {
        return "Coordinates{" +
                "start=" + start +
                ", end=" + end +
                ", Nterm=" + Nterm +
                ", Cterm=" + Cterm +
                '}';
    }
}