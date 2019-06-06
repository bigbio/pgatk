package org.bigbio.pgatk.pepgenome;


import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.Before;
import org.junit.FixMethodOrder;
import org.junit.Test;
import org.junit.runners.MethodSorters;

import java.io.File;
import java.io.IOException;
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
 * @author ypriverol on 28/09/2018.
 */
@FixMethodOrder(MethodSorters.NAME_ASCENDING)
public class PepGenomeToolTest {


    private static final Logger log = Logger.getLogger(PepGenomeToolTest.class);

    String fileIn = null;
    String fileFasta = null;
    String fileGTF = null;
    String fileCPogo = null;


    @Before
    public void setUp() throws Exception {

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/Testfile_small.txt")).toURI()).getAbsolutePath();
        fileFasta = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/minimal_gencode.v25.pc_translations.fa")).toURI()).getAbsolutePath();
        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());

        fileGTF = TestUtils.unGzip(inputGZfile).getAbsolutePath();

        fileCPogo = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/cpogo/Testfile_small.bed")).toURI()).getAbsolutePath();

    }

    @Test
    public void mainInMemory() throws IOException {
        log.info("InMemoryTest");
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

        File outputBed = new File(fileIn.replace(".txt", ".bed"));
        File cPogoBed = new File(fileCPogo);

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = TestUtils.getBedLines(cPogoBed);
        Assert.assertEquals(bedLines.size(), cPogoLines.size());

        for (List<String> bedLine1 : bedLines) {
            boolean found = false;
            for (List<String> cbedLine : cPogoLines) {
                if (bedLine1.get(3).equalsIgnoreCase(cbedLine.get(3))) {
                    found = compareBedLines(bedLine1, cbedLine);
                    log.info(bedLine1.get(3) + " -- " + found);
                    if (found)
                        break;
                }
            }
            Assert.assertTrue(found);
        }

        deleteAfterTest();
        log.info(" ");

    }

    @Test
    public void mainInDB() throws IOException {
        log.info("InDBTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-inm");
        argList.add("1");

        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".txt", ".bed"));
        File cPogoBed = new File(fileCPogo);

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = TestUtils.getBedLines(cPogoBed);
        Assert.assertEquals(bedLines.size(), cPogoLines.size());

        for (List<String> bedLine1 : bedLines) {
            boolean found = false;
            for (List<String> cbedLine : cPogoLines) {
                if (bedLine1.get(3).equalsIgnoreCase(cbedLine.get(3))) {
                    found = compareBedLines(bedLine1, cbedLine);
                    log.info(bedLine1.get(3) + " -- " + found);
                    if (found)
                        break;
                }
            }
            Assert.assertTrue(found);
        }

        deleteAfterTest();
        log.info(" ");

    }

    @Test
    public void mainInMemorySparkMode() throws IOException {
        log.info("InMemoryTest-SparkMode");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-spark_master");
        argList.add("local[*]");

        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".txt", ".bed"));
        File cPogoBed = new File(fileCPogo);

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = TestUtils.getBedLines(cPogoBed);
        Assert.assertEquals(bedLines.size(), cPogoLines.size());

        for (List<String> bedLine1 : bedLines) {
            boolean found = false;
            for (List<String> cbedLine : cPogoLines) {
                if (bedLine1.get(3).equalsIgnoreCase(cbedLine.get(3))) {
                    found = compareBedLines(bedLine1, cbedLine);
                    log.info(bedLine1.get(3) + " -- " + found);
                    if (found)
                        break;
                }
            }
            Assert.assertTrue(found);
        }

        deleteAfterTest();
        log.info(" ");

    }

    @Test
    public void mainInDBSparkMode() throws IOException {
        log.info("InDBTest-SparkMode");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-inm");
        argList.add("1");
        argList.add("-spark_master");
        argList.add("local[*]");

        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".txt", ".bed"));
        File cPogoBed = new File(fileCPogo);

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = TestUtils.getBedLines(cPogoBed);
        Assert.assertEquals(bedLines.size(), cPogoLines.size());

        for (List<String> bedLine1 : bedLines) {
            boolean found = false;
            for (List<String> cbedLine : cPogoLines) {
                if (bedLine1.get(3).equalsIgnoreCase(cbedLine.get(3))) {
                    found = compareBedLines(bedLine1, cbedLine);
                    log.info(bedLine1.get(3) + " -- " + found);
                    if (found)
                        break;
                }
            }
            Assert.assertTrue(found);
        }

        deleteAfterTest();
        log.info(" ");

    }

    private boolean compareBedLines(List<String> bedLine, List<String> cbedLine) {
        return cbedLine.containsAll(bedLine);
    }

    private void deleteAfterTest() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.delete();
    }

}
