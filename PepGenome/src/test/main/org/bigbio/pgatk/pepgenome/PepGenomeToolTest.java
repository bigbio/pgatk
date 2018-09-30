package org.bigbio.pgatk.pepgenome;


import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

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

public class PepGenomeToolTest {

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
    public void main() throws IOException {
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

        File outputBed = new File(fileIn.replace(".txt", ".bed"));
        File cPogoBed = new File(fileCPogo);

        List<List<String>> bedLines = getBedLines(outputBed);
        Assert.assertTrue(bedLines.size() == 29);


        List<List<String>> cPogoLines = getBedLines(cPogoBed);
        Assert.assertTrue(bedLines.size() == cPogoLines.size());

        for(int i = 0; i < bedLines.size(); i++){
            boolean found = false;
            List<String> bedLine = bedLines.get(i);
            for(int j = 0; j < cPogoLines.size(); j ++){
                List<String> cbedLine = cPogoLines.get(j);
                if(bedLine.get(3).equalsIgnoreCase(cbedLine.get(3))){
                    found = compareBedLines( bedLine, cbedLine);
                    System.out.println(bedLine.get(3) + " -- " + found);
                    if(found)
                        break;
                }
            }
            Assert.assertTrue(found);
        }

    }

    private boolean compareBedLines(List<String> bedLine, List<String> cbedLine) {
        return  bedLine.stream().allMatch(num -> cbedLine.contains(num));
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

    private List<List<String>> getBedLines(File inFile) throws IOException {
        BufferedReader buf = new BufferedReader(new FileReader(inFile));
        List<List<String>> bedLines = new ArrayList<>();
        String line = null;
        while((line = buf.readLine()) != null){
            String[] lines = line.split("\t");
            bedLines.add(Arrays.asList(lines));
        }
        return bedLines;
    }
}
