package org.bigbio.pgatk.pepgenome;

import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

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

        fileFasta = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("mztab/gencode.v25.pc_translations.fa.gz")).toURI())).getAbsolutePath();

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("mztab/sample.mztab")).toURI()).getAbsolutePath();

        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());
        fileGTF = TestUtils.unGzip(inputGZfile).getAbsolutePath();

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

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(4577, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    @Test
    public void mzTabTestMissMatches() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("mm");
        argList.add("2");
        argList.add("mmmode");
        argList.add("true");
        argList.add("-inf");
        argList.add("mztab");


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".txt", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(4577, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }


    private void deleteOnExits() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.deleteOnExit();
    }
}
