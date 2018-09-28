package org.bigbio.pgatk.pepgenome;

import org.bigbio.pgatk.pepgenome.common.GenomeCoordinates;
import org.bigbio.pgatk.pepgenome.common.PeptideCoordinates;
import org.bigbio.pgatk.pepgenome.common.PeptidecoordsPairPcompare;
import javafx.util.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

//the PTMEntry contains information about post translational modifications like methlation.
//and at which position they appear.
public class PTMEntry {
    //holds the name of the modification
    private String m_ptm_psi_name;
    //holds the lowest start coord
    private int m_peptide_start_coord;
    //holds the highest end coord.
    private int m_peptide_end_coord;
    //holds all added ptms and their genomic locations.
    private TreeSet<Pair<PeptideCoordinates, GenomeCoordinates>> m_ptm_coord = new TreeSet<>(new PeptidecoordsPairPcompare());

    public PTMEntry() {
        this.m_peptide_start_coord = -1;
        this.m_peptide_end_coord = -1;
    }

    public PTMEntry(String name, int start, int end) {
        this.m_ptm_psi_name = name;
        this.m_peptide_start_coord = start;
        this.m_peptide_end_coord = end;
    }

    public PTMEntry(PTMEntry rhs) {
        this.m_ptm_psi_name = rhs.m_ptm_psi_name;
        this.m_peptide_start_coord = rhs.m_peptide_start_coord;
        this.m_peptide_end_coord = rhs.m_peptide_end_coord;
        this.m_ptm_coord = rhs.m_ptm_coord;
    }

    @Override
    public boolean equals(Object o) {
        return equalsTo((PTMEntry) o);
    }

    //returns true if the ptm_psi_name, start and end coords are the same,
    //otherwise returns false.
    public boolean equalsTo(PTMEntry rhs) {
        return m_ptm_psi_name.equals(rhs.m_ptm_psi_name) && m_peptide_start_coord == rhs.m_peptide_start_coord && m_peptide_end_coord == rhs.m_peptide_end_coord;
    }

    //returns a pair of start and end coord for the current PTM.
    public final Pair<Integer, Integer> get_range() {
        return new Pair<>(m_peptide_start_coord, m_peptide_end_coord);
    }

    //decreases m_peptide_start coord if coord is smaller than that
    //or increases m_peptide start coor if coord is larger than that.
    //sets both to coord if they havent been set before.
    public final void add_coord(int coord) {
        if (m_peptide_start_coord == -1 && m_peptide_end_coord == -1) {
            m_peptide_start_coord = coord;
            m_peptide_end_coord = coord;
        } else if (coord <= m_peptide_start_coord) {
            m_peptide_start_coord = coord;
        } else if (coord >= m_peptide_end_coord) {
            m_peptide_end_coord = coord;
        }
    }

    //returns a list of all found genome coordinates.
    public final List<Pair<PeptideCoordinates, GenomeCoordinates>> get_genome_coordinates() {
        return new ArrayList<>(m_ptm_coord);
    }

    //adds genome coordinates.
    public final void add_genome_coordinates(PeptideCoordinates coords, GenomeCoordinates ptmcoords) {
        m_ptm_coord.add(new Pair<>(coords, ptmcoords));
    }
}