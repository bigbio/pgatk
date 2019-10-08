package org.bigbio.pgatk.pepgenome.common;



import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * The PTMEntry contains information about post translational modifications
 * like methlation and at which position they appear.
 *
 * @author ypriverol
 */

public class PTMEntry implements Serializable {
    private static final long serialVersionUID = -7103530404088146538L;

    //holds the name of the modification
    private String psiName;

    //holds the lowest start coord
    private int peptideStartCoord;

    //holds the highest end coord.
    private int peptideEndCoord;

    //holds all added ptms and their genomic locations.
    private TreeSet<Tuple<PeptideCoordinates, GenomeCoordinates>> ptmCoord = new TreeSet<>(new PeptidecoordsPairPcompare());

    public PTMEntry() {
        this.peptideStartCoord = -1;
        this.peptideEndCoord = -1;
    }

    public PTMEntry(String name, int start, int end) {
        this.psiName = name;
        this.peptideStartCoord = start;
        this.peptideEndCoord = end;
    }

    public PTMEntry(PTMEntry rhs) {
        this.psiName = rhs.psiName;
        this.peptideStartCoord = rhs.peptideStartCoord;
        this.peptideEndCoord = rhs.peptideEndCoord;
        this.ptmCoord = rhs.ptmCoord;
    }

    @Override
    public boolean equals(Object obj) {
        return equalsTo((PTMEntry) obj);
    }

    //returns true if the ptm_psi_name, start and end coords are the same,
    //otherwise returns false.
    public boolean equalsTo(PTMEntry rhs) {
        return psiName.equals(rhs.psiName) && peptideStartCoord == rhs.peptideStartCoord && peptideEndCoord == rhs.peptideEndCoord;
    }

    //returns a pair of start and end coord for the current PTM.
    public final Tuple<Integer, Integer> get_range() {
        return new Tuple<>(peptideStartCoord, peptideEndCoord);
    }

    //decreases m_peptide_start coord if coord is smaller than that
    //or increases m_peptide start coor if coord is larger than that.
    //sets both to coord if they havent been set before.
    public final void add_coord(int coord) {
        if (peptideStartCoord == -1 && peptideEndCoord == -1) {
            peptideStartCoord = coord;
            peptideEndCoord = coord;
        } else if (coord <= peptideStartCoord) {
            peptideStartCoord = coord;
        } else if (coord >= peptideEndCoord) {
            peptideEndCoord = coord;
        }
    }

    //returns a list of all found genome coordinates.
    public final List<Tuple<PeptideCoordinates, GenomeCoordinates>> get_genome_coordinates() {
        return new ArrayList<>(ptmCoord);
    }

    //adds genome coordinates.
    public final void add_genome_coordinates(PeptideCoordinates coords, GenomeCoordinates ptmcoords) {
        ptmCoord.add(new Tuple<>(coords, ptmcoords));
    }
}