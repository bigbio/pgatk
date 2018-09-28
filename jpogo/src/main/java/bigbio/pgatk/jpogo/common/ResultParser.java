package bigbio.pgatk.jpogo.common;

import bigbio.pgatk.jpogo.CoordinateWrapper;
import bigbio.pgatk.jpogo.KmerMap;
import bigbio.pgatk.jpogo.MappedPeptides;
import bigbio.pgatk.jpogo.PeptideEntry;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

//this class parses peptide input files and maps them to the genome.
public class ResultParser {

    //read function. this reads the peptides input and sets the wheels in motion.
    //this function will set the wheels in motion to find the peptides in the proteins.
    public static void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, KmerMap k) throws Exception {

        FileInputStream ifs = new FileInputStream(file);
        BufferedReader reader = new BufferedReader(new InputStreamReader(ifs));
        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        String peptide_string;
        String tissue;
        String iso_seq_without_ptms;
        int sig_PSMs;
        double quant;
        Map<String, TranscriptsT> gene_id_map;

        String line;
        while ((line = reader.readLine()) != null) {
            if ((line.startsWith("Experiment")) || (line.startsWith("Sample"))) {
                continue;
            }
            ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t", false)));
            //using only the tokens needed.
            tissue = tokens.get(0);
            peptide_string = tokens.get(1);
            //std::cout << "Mapping following peptide: " << peptide_string << std::endl;
            sig_PSMs = Integer.parseInt(tokens.get(2));
            quant = Double.parseDouble(tokens.get(3));

            //clearing the tokens list for the next iteration.
            tokens.clear();

            if (sig_PSMs > 0) {
                //the matching will only use the amino acids.
                iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                if (!coordwrapper.peptide_already_exists(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sig_PSMs, gene_id_map.size(), ofs, quant, it);
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sig_PSMs, quant);
                    }
                }
            }
        }
        ofs.close();
        reader.close();
        ifs.close();
    }
}