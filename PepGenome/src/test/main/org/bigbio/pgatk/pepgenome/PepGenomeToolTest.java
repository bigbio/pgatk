package org.bigbio.pgatk.pepgenome;


import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.zip.GZIPInputStream;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 *
 * @author ypriverol on 28/09/2018.
 */

public class PepGenomeToolTest {

    String fileIn = null;
    String fileFasta = null;
    String fileGTF = null;

    @Before
    public void setUp() throws Exception {

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/Testfile_small.txt")).toURI()).getAbsolutePath();
        fileFasta = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/minimal_gencode.v25.pc_translations.fa")).toURI()).getAbsolutePath();
        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());

        fileGTF = unGzip(inputGZfile).getAbsolutePath();

    }

    @Test
    public void main() {
        System.out.println("main");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);

        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

    }

    public static File unGzip(File infile) throws IOException {
        GZIPInputStream gin = new GZIPInputStream(new FileInputStream(infile));
        FileOutputStream fos = null;
        try {
            File outFile = new File(infile.getParent(), infile.getName().replaceAll("\\.gz$", ""));
            fos = new FileOutputStream(outFile);
            byte[] buf = new byte[100000];
            int len;
            while ((len = gin.read(buf)) > 0) {
                fos.write(buf, 0, len);
            }

            fos.close();
            return outFile;
        } finally {
            if (gin != null) {
                gin.close();
            }
            if (fos != null) {
                fos.close();
            }
        }
    }
}
