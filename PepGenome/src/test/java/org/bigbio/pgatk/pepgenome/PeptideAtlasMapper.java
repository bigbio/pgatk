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
 * @author ypriverol on 05/10/2018.
 */
public class PeptideAtlasMapper {

    private static final Logger log = Logger.getLogger(PepGenomeToolTest.class);

    String fileIn = null;
    String fileFasta = null;
    String fileGTF = null;
    private String fileUnmappedIn;
    private String fileInZebrafish;
    private String fileFastaZebrafish;
    private String fileGTFZebrafish;


    @Before
    public void setUp() throws Exception {

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("peptideatlas/peptideatlas-500.tsv")).toURI()).getAbsolutePath();
        fileUnmappedIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("peptideatlas/peptideatlas-unmapped-500.tsv")).toURI()).getAbsolutePath();
        fileFasta = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("mztab/gencode.v25.pc_translations.fa.gz")).toURI())).getAbsolutePath();
        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());
        fileGTF = TestUtils.unGzip(inputGZfile).getAbsolutePath();

        fileInZebrafish = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/taxon-different.tsv")).toURI()).getAbsolutePath();
        fileFastaZebrafish = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Danio_rerio.GRCz11.pep.all.fa.gz")).toURI())).getAbsolutePath();
        inputGZfile = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Danio_rerio.GRCz11.94.gtf.gz")).toURI()));
        fileGTFZebrafish = inputGZfile.getAbsolutePath();


    }

    @Test
    public void peptideAtlasTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-inf");
        argList.add("peptideatlas");


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileIn.replace(".tsv", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(109869, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    @Test
    public void peptideUnMappedTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileUnmappedIn);
        argList.add("-fasta");
        argList.add(fileFasta);
        argList.add("-gtf");
        argList.add(fileGTF);
        argList.add("-inf");
        argList.add("peptideatlas");


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileUnmappedIn.replace(".tsv", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(0, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    @Test
    public void zebrafishTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileInZebrafish);
        argList.add("-fasta");
        argList.add(fileFastaZebrafish);
        argList.add("-gtf");
        argList.add(fileGTFZebrafish);
        argList.add("-inf");
        argList.add("peptideatlas");
        argList.add("-species");
        argList.add("7955");


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileInZebrafish.replace(".tsv", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(0, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    private void deleteOnExits() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.deleteOnExit();
    }
}
