package org.bigbio.pgatk.pepgenome.io;

import org.apache.log4j.Logger;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import uk.ac.ebi.pride.jmztab.model.PSM;
import uk.ac.ebi.pride.jmztab.utils.MZTabFileParser;

import java.io.*;
import java.util.ArrayList;
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
 * @author ypriverol on 03/10/2018.
 */
public class MzTabInputPeptideFileParser implements PeptideInputReader{

    private static final org.apache.log4j.Logger log = Logger.getLogger(MzTabInputPeptideFileParser.class);


    /**
     * Read function. this reads the peptides input and sets the wheels in motion.
     * this function will set the wheels in motion to find the peptides in the proteins.
     *
     * @param file input file.
     */
    public  void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        File mzTabFile = new File(file);
        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        if (!mzTabFile.canRead()) {
            log.error("could not read '" + file + "'.");
            throw new IOException("The file doesn't not exists -- " + file);
        }

        try {
            MZTabFileParser parser = new MZTabFileParser(mzTabFile, new FileOutputStream(mzTabFile.getAbsolutePath() + "errors.out"));

            String peptideString;
            String isoSeqWithoutPtms;
            int sigPSMs;
            double quant;
            Map<String, TranscriptsT> gene_id_map;

            for(PSM psm: parser.getMZTabFile().getPSMs()){
                peptideString = psm.getSequence();
                quant = 1.0;
                sigPSMs = 1;
                //the matching will only use the amino acids.
                isoSeqWithoutPtms = Utils.make_iso_sequence(Utils.remove_ptms(peptideString));

                if (!coordwrapper.isPeptidePresent(isoSeqWithoutPtms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    gene_id_map = k.find_peptide(isoSeqWithoutPtms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptideString, file, sigPSMs, gene_id_map.size(), ofs, quant, it);
                    }
                    if (gene_id_map.isEmpty()){
                        ofs.write(("No-Gene" + "\t" + peptideString + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + file + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(isoSeqWithoutPtms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptideString, file, sigPSMs, quant);
                    }
                }
            }
        } catch (IOException e) {
            log.error("Could not create mzTab file reader", e);
            throw new IOException("The file doesn't not exists -- " + file);
        }

    }
}
