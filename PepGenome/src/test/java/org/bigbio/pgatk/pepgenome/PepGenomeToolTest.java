package org.bigbio.pgatk.pepgenome;


import org.apache.log4j.Logger;
import org.junit.*;
import org.junit.runners.MethodSorters;


import java.io.*;
import java.util.*;
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

        fileGTF = unGzip(inputGZfile).getAbsolutePath();

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

        List<List<String>> bedLines = getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = getBedLines(cPogoBed);
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

        deleteOnExits();

    }

    private void deleteOnExits() {
        String fileBed = fileIn.replaceAll(".txt", ".bed");
        File fileInput = new File(fileBed);
        fileInput.deleteOnExit();
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

        List<List<String>> bedLines = getBedLines(outputBed);
        Assert.assertEquals(29, bedLines.size());


        List<List<String>> cPogoLines = getBedLines(cPogoBed);
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

        deleteOnExits();

    }
    private boolean compareBedLines(List<String> bedLine, List<String> cbedLine) {
        return  bedLine.stream().allMatch(cbedLine::contains);
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
