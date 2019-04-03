package org.bigbio.pgatk.pepgenome.io.custom;

import org.apache.log4j.Logger;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.io.PeptideInputReader;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

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
public class PeptideAtlasPeptideParser implements PeptideInputReader {

    private static final Logger log = Logger.getLogger(PeptideAtlasPeptideParser.class);


    /**
     * Read function. this reads the peptides input and sets the wheels in motion.
     * this function will set the wheels in motion to find the peptides in the proteins.
     *
     * @param file input file.
     */
    public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        File mzTabFile = new File(file);
        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        FileInputStream ifs = new FileInputStream(file);
        BufferedReader reader = new BufferedReader(new InputStreamReader(ifs));

        if (!mzTabFile.canRead()) {
            log.error("could not read '" + file + "'.");
            throw new IOException("The file doesn't not exists -- " + file);
        }

        try {

            String peptide_string;
            String iso_seq_without_ptms;
            String tissue;
            int sigPSMs;
            double quant;
            Map<String, TranscriptsT> gene_id_map;

            String line;
            while ((line = reader.readLine()) != null) {
                if ((line.startsWith("peptide_accession"))) {
                    continue;
                }
                ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t", false)));
                //using only the tokens needed.
                tissue = tokens.get(12);
                peptide_string = tokens.get(1);
                sigPSMs = Integer.parseInt(tokens.get(3));
                quant = Double.parseDouble(tokens.get(3));
                //the matching will only use the amino acids.
                iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));
                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it);
                    }
                    if (gene_id_map.isEmpty()){
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, file, sigPSMs, quant);
                    }
                }
            }
            ofs.flush();
            ofs.close();

        } catch (IOException e) {
            log.error("Could not create mzTab file reader", e);
            throw new IOException("The file doesn't not exists -- " + file);
        }

    }
}
