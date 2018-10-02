package org.bigbio.pgatk.pepgenome.kmer.indb;

import org.apache.commons.lang3.SerializationUtils;
import org.bigbio.pgatk.pepgenome.kmer.IKmerEntry;
import org.bigbio.pgatk.pepgenome.kmer.KmerEntry;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DataInput2;
import org.mapdb.DataOutput2;
import org.mapdb.Serializer;

import java.io.IOException;
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
public class KmerSerializer implements Serializer<IKmerEntry>, Serializable {

    private static final long serialVersionUID = 7188928270499043150L;

    @Override
    public void serialize(@NotNull DataOutput2 dataOutput2, @NotNull IKmerEntry iKmerEntry) throws IOException {
        byte[] byteProtein = SerializationUtils.serialize(iKmerEntry);
        dataOutput2.write(byteProtein);
    }

    @Override
    public IKmerEntry deserialize(@NotNull DataInput2 dataInput2, int i) throws IOException {
        return (KmerEntry) SerializationUtils.deserialize(dataInput2.internalByteArray());
    }
}
