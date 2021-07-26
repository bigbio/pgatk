package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.PepGenomeTool;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * This class parses peptide input files from Tab delimited files. This classes only read the tab delimited files define by
 * Sample | Peptide | PSMs | Quant
 *
 * @author ypriverol
 */

public class TabInputPeptideFileParser implements PeptideInputReader, Serializable {

    public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {
      normalRead(file, coordwrapper, mapping, unmappedoutput, k);
    }

    //read function. this reads the peptides input and sets the wheels in motion.
    //this function will set the wheels in motion to find the peptides in the proteins.
    private void normalRead(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        FileInputStream ifs = new FileInputStream(file);
        BufferedReader reader = new BufferedReader(new InputStreamReader(ifs));
        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        String peptide_string;
        String tissue;
        String iso_seq_without_ptms;
        int sigPSMs;
        double quant;
        Map<String, TranscriptsT> gene_id_map;

        // ||For use with peptide filter mode||
        int allowedMismatches = 0;
        String targetTranscriptIDs = "";


        String line;
        while ((line = reader.readLine()) != null) {

            if ((line.toLowerCase().startsWith("experiment")) || (line.toLowerCase().startsWith("sample"))) {
                continue;
            }
            // Tokenize peptide input
            ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t", false)));
            //using only the tokens needed.
            // Extracts tissue
            tissue = tokens.get(0).trim();
            // Extracts peptide string
            peptide_string = tokens.get(1).trim();
            // Extracts PSMs
            String sigPsmStr = tokens.get(2).trim();
            if (sigPsmStr.length() == 0) {
                sigPsmStr = "0";
            }
            sigPSMs = Integer.parseInt(sigPsmStr);
            // Extracts Quant
            String quantStr = tokens.get(3).trim();
            if (quantStr.length() == 0) {
                quantStr = "0";
            }
            quant = Double.parseDouble(quantStr);

            // If peptide filter mode is on, and the peptide line has six fields filled.
            if (PepGenomeTool.usePeptideFilter && tokens.size()==6) {

                // Extract Allowed Mismatches
                allowedMismatches = Integer.parseInt(tokens.get(4).trim());

                // Extract Target Trans ID - If left blank, map to all
                targetTranscriptIDs = tokens.get(5).trim();
                if (targetTranscriptIDs.equals("")) {
                    targetTranscriptIDs = "all";
                }

            }

            // Clearing the tokens list for the next iteration.
            tokens.clear();

            if (sigPSMs > 0) {
                // Remove PTMs from the iso sequence, only using amino acids
                iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    // Match peptide, produce gene id map using KmerTreeMap/KmerSortedMap
                    if (PepGenomeTool.usePeptideFilter) {
                        // Using peptide filter mode
                        gene_id_map = k.find_peptide(iso_seq_without_ptms, targetTranscriptIDs, allowedMismatches);
                    } else {
                        // Default
                        gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    }

                    // Note: k.getIsVariant: Extract variant status (check whether peptide contains mismatches)||
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {

                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it, k.getIsVariant());
                    }
                    if (gene_id_map.isEmpty()) {
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    // if the peptide already exists its genomic coordinates dont have to be recalculated.
                    // only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sigPSMs, quant, k.getIsVariant());
                    }
                }
            }
        }
        ofs.close();
        reader.close();
        ifs.close();
    }

}
