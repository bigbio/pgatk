package org.bigbio.pgatk.pepgenome;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
 * @author ypriverol on 05/10/2018.
 */
public class TestUtils {

    public static List<List<String>> getBedLines(File inFile) throws IOException {
        BufferedReader buf = new BufferedReader(new FileReader(inFile));
        List<List<String>> bedLines = new ArrayList<>();
        String line;
        while((line = buf.readLine()) != null){
            String[] lines = line.split("\t");
            bedLines.add(Arrays.asList(lines));
        }
        return bedLines;
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
}
