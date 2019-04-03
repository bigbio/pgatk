package org.bigbio.pgatk.pepgenome.kmer;

import org.bigbio.pgatk.pepgenome.common.ProteinEntry;

import java.io.Serializable;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 *
 * @author ypriverol on 01/10/2018.
 */
public interface IKmerEntry extends Serializable {

    /**
     * Retrieve the Protein entry.
     * @return ProteinEntry
     */
    ProteinEntry m_p_protein();

    /**
     * Retrieve the starting position of the kmer in the Protein sequence
     * @return //the (0 based) index of the first letter of the kmer in the protein string.
     */
    int m_pos_in_protein();
}
