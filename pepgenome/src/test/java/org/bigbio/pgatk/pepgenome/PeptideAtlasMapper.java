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
//    String fileGenomeFasta = null;
    String fileGTF = null;
    private String fileUnmappedIn;
    private String fileInZebrafish;
    private String fileFastaZebrafish;
//    private String fileGenomeFastaZebrafish;
    private String fileGTFZebrafish;
    private String fileInYeast;
    private String fileFastaYeast;
//    private String fileGenomeFastaYeast;
    private String fileGTFYeast;
    private String fileInBonobo;
    private String fileFastaBonobo;
//    private String fileGenomeFastaBonobo;
    private String fileGTFBonobo;
    private String fileInAlpaca;
    private String fileFastaAlpaca;
//    private String fileGenomeFastaAlpaca;
    private String fileGTFAlpaca;
    


    @Before
    public void setUp() throws Exception {

        fileIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("peptideatlas/peptideatlas-500.tsv")).toURI()).getAbsolutePath();
        fileUnmappedIn = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("peptideatlas/peptideatlas-unmapped-500.tsv")).toURI()).getAbsolutePath();
        fileFasta = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("mztab/gencode.v25.pc_translations.fa.gz")).toURI())).getAbsolutePath();
//        fileGenomeFasta = TestUtils.unGzip(new File(Objects.requireNonNull(new File("D:/Data/Genomes/DNA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")).toURI())).getAbsolutePath();
        File inputGZfile = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("small/gencode.v25.annotation.gtf.gz")).toURI());
        fileGTF = TestUtils.unGzip(inputGZfile).getAbsolutePath();

        fileInZebrafish = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/taxon-different.tsv")).toURI()).getAbsolutePath();
        fileFastaZebrafish = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Danio_rerio.GRCz11.pep.all.fa.gz")).toURI())).getAbsolutePath();
//        fileGenomeFastaZebrafish = TestUtils.unGzip(new File(Objects.requireNonNull(new File("D:/Data/Genomes/DNA/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz")).toURI())).getAbsolutePath();
        inputGZfile = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Danio_rerio.GRCz11.94.gtf.gz")).toURI()));
        fileGTFZebrafish = inputGZfile.getAbsolutePath();
        
        fileInYeast = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/yeast.tsv")).toURI()).getAbsolutePath();
        fileFastaYeast = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz")).toURI())).getAbsolutePath();
//        fileGenomeFastaYeast = TestUtils.unGzip(new File(Objects.requireNonNull(new File("D:/Data/Genomes/DNA/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz")).toURI())).getAbsolutePath();
        inputGZfile = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz")).toURI()));
        fileGTFYeast = inputGZfile.getAbsolutePath();
        
        fileInBonobo = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/bonobo.pogo")).toURI()).getAbsolutePath();
        fileFastaBonobo = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Pan_paniscus.panpan1.1.pep.all.fa.gz")).toURI())).getAbsolutePath();
//        fileGenomeFastaBonobo = TestUtils.unGzip(new File(Objects.requireNonNull(new File("D:/Data/Genomes/DNA/Pan_paniscus.panpan1.1.dna.toplevel.fa.gz")).toURI())).getAbsolutePath();
        inputGZfile = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Pan_paniscus.panpan1.1.94.gtf.gz")).toURI()));
        fileGTFBonobo = inputGZfile.getAbsolutePath();
        
        fileInAlpaca = new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/alpaca.pogo")).toURI()).getAbsolutePath();
        fileFastaAlpaca = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Vicugna_pacos.vicPac1.pep.all.fa.gz")).toURI())).getAbsolutePath();
//        fileGenomeFastaAlpaca = TestUtils.unGzip(new File(Objects.requireNonNull(new File("D:/Data/Genomes/DNA/Vicugna_pacos.vicPac1.dna.toplevel.fa.gz")).toURI())).getAbsolutePath();
        inputGZfile = TestUtils.unGzip(new File(Objects.requireNonNull(PepGenomeToolTest.class.getClassLoader().getResource("taxonomies/Vicugna_pacos.vicPac1.94.gtf.gz")).toURI()));
        fileGTFAlpaca = inputGZfile.getAbsolutePath();


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
//        argList.add("-genome");
//        argList.add(fileGenomeFasta);


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
//        argList.add("-genome");
//        argList.add(fileGenomeFastaZebrafish);


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileInZebrafish.replace(".tsv", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(7, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }
    
    @Test
    public void yeastTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileInYeast);
        argList.add("-fasta");
        argList.add(fileFastaYeast);
        argList.add("-gtf");
        argList.add(fileGTFYeast);
        argList.add("-inf");
        argList.add("peptideatlas");
//        argList.add("-genome");
//        argList.add(fileGenomeFastaYeast);


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileInYeast.replace(".tsv", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(3, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }
    
    @Test
    public void bonoboTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileInBonobo);
        argList.add("-fasta");
        argList.add(fileFastaBonobo);
        argList.add("-gtf");
        argList.add(fileGTFBonobo);
        argList.add("-inf");
        argList.add("tab");
//        argList.add("-genome");
//        argList.add(fileGenomeFastaBonobo);


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileInBonobo.replace(".pogo", ".bed"));

        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        Assert.assertEquals(7, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }
    
    @Test
    public void alpacaTest() throws IOException {
        log.info("InMemoryTest");
        List<String> argList = new ArrayList<>();

        argList.add("-in");
        argList.add(fileInAlpaca);
        argList.add("-fasta");
        argList.add(fileFastaAlpaca);
        argList.add("-gtf");
        argList.add(fileGTFAlpaca);
        argList.add("-inf");
        argList.add("tab");
//        argList.add("-genome");
//        argList.add(fileGenomeFastaAlpaca);


        String[] args = new String[argList.size()];
        argList.toArray(args);
        PepGenomeTool.main(args);

        File outputBed = new File(fileInAlpaca.replace(".pogo", ".bed"));
        List<List<String>> bedLines = TestUtils.getBedLines(outputBed);
        
        if(argList.contains("-genome")) {
        	outputBed = new File(fileInAlpaca.replace(".pogo", "_patch_hapl_scaff.bed"));
        	bedLines.addAll(TestUtils.getBedLines(outputBed));
        }
        
        Assert.assertEquals(2, bedLines.size());

        deleteOnExits();
        log.info(" ");

    }

    private void deleteOnExits() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.deleteOnExit();
    }
}
