package org.bigbio.pgatk.pepgenome;

import org.apache.log4j.Logger;
import org.junit.Assert;
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
 * @author ypriverol on 03/10/2018.
 */
public class PepGeonomeMzTabTest  {


    private static final Logger log = Logger.getLogger(PepGenomeToolTest.class);

    String fileIn = null;
    String fileFasta = null;
    String fileGTF = null;


    @Before
    public void setUp() throws Exception {

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("mztab/control_exo_rep1_high_mol_weight.dat-pride.mztab")).toURI()).getAbsolutePath();
        fileFasta = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/minimal_gencode.v25.pc_translations.fa")).toURI()).getAbsolutePath();
        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());
        fileGTF = unGzip(inputGZfile).getAbsolutePath();

    }

    @Test
    public void mzTabTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-inf");
        argList.add("mztab");


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".txt", ".bed"));

        List<List<String>> bedLines = getBedLines(outputBed);
        Assert.assertEquals(4577, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    private void deleteOnExits() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.deleteOnExit();
    }

    public static File unGzip(File infile) throws IOException {
        File outFile = new File(infile.getParent(), infile.getName().replaceAll("\\.gz$", ""));
        try (GZIPInputStream gin = new GZIPInputStream(new FileInputStream(infile)); FileOutputStream fos = new FileOutputStream(outFile)) {
            byte[] buf = new byte[100000];
            int len;
            while ((len = gin.read(buf)) > 0) {
                fos.write(buf, 0, len);
            }

            fos.close();
            return outFile;
        }
    }

    private List<List<String>> getBedLines(File inFile) throws IOException {
        BufferedReader buf = new BufferedReader(new FileReader(inFile));
        List<List<String>> bedLines = new ArrayList<>();
        String line;
        while((line = buf.readLine()) != null){
            String[] lines = line.split("\t");
            bedLines.add(Arrays.asList(lines));
        }
        return bedLines;
    }
}
