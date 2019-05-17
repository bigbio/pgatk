package org.bigbio.pgatk.pepgenome.io;

import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SparkSession;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.SparkConfig;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
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

public class TabInputPeptideFileParser implements PeptideInputReader {

    public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {
        if (SparkConfig.getInstance().getMaster() != null) {
            sparkRead(file, coordwrapper, mapping, unmappedoutput, k);
        } else {
            normalRead(file, coordwrapper, mapping, unmappedoutput, k);
        }

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

        String line;
        while ((line = reader.readLine()) != null) {
            if ((line.toLowerCase().startsWith("experiment")) || (line.toLowerCase().startsWith("sample"))) {
                continue;
            }
            ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t", false)));
            //using only the tokens needed.
            tissue = tokens.get(0);
            peptide_string = tokens.get(1);
            //std::cout << "Mapping following peptide: " << peptide_string << std::endl;
            sigPSMs = Integer.parseInt(tokens.get(2));
            quant = Double.parseDouble(tokens.get(3));

            //clearing the tokens list for the next iteration.
            tokens.clear();

            if (sigPSMs > 0) {
                //the matching will only use the amino acids.
                iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it);
                    }
                    if (gene_id_map.isEmpty()) {
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sigPSMs, quant);
                    }
                }
            }
        }
        ofs.close();
        reader.close();
        ifs.close();
    }

    private void sparkRead(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        SparkSession sparkSession = SparkSession.builder().
                master(SparkConfig.getInstance().getMaster())
                .appName("pgatk tab input file parser")
                .getOrCreate();

        Dataset<Row> tsv = sparkSession.read()
                .option("sep", "\t")
//                .option("header", "true")
                .csv(file);

        List<Row> rows = tsv.collectAsList();
        for (Row r : rows) {
            String tissue = r.getString(0);
            if ((tissue.toLowerCase().startsWith("experiment")) || (tissue.toLowerCase().startsWith("sample"))) {
                continue;
            }
            String peptide_string = r.getString(1);
            int sigPSMs = Integer.parseInt(r.getString(2));
            double quant = Double.parseDouble(r.getString(3));

            if (sigPSMs > 0) {
                //the matching will only use the amino acids.
                String iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    Map<String, TranscriptsT> gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it);
                    }
                    if (gene_id_map.isEmpty()) {
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sigPSMs, quant);
                    }
                }
            }
        }

        ofs.close();
    }
}