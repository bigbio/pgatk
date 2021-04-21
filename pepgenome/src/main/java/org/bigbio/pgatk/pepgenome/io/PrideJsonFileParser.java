package org.bigbio.pgatk.pepgenome.io;

import io.github.bigbio.pgatk.io.pride.ArchiveSpectrum;
import io.github.bigbio.pgatk.io.pride.PrideJsonIterableReader;
import lombok.extern.slf4j.Slf4j;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;
import uk.ac.ebi.pride.jmztab.model.PSM;
import uk.ac.ebi.pride.jmztab.utils.MZTabFileParser;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;

@Slf4j
public class PrideJsonFileParser implements PeptideInputReader, Serializable {


  @Override
  public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

    File prideJsonFile = new File(file);
    FileOutputStream ofs = new FileOutputStream(unmappedoutput);

    if (!prideJsonFile.canRead()) {
      log.error("could not read '" + file + "'.");
      throw new IOException("The file doesn't not exists -- " + file);
    }

    try {

      PrideJsonIterableReader reader = new PrideJsonIterableReader(prideJsonFile);

      String peptideString;
      String isoSeqWithoutPtms;
      int sigPSMs;
      double quant;
      Map<String, TranscriptsT> gene_id_map;

      while(reader.hasNext()){

        ArchiveSpectrum spectrum = (ArchiveSpectrum) reader.next();

        peptideString = spectrum.getPeptideSequence();
        quant = spectrum.getPeptideIntensity() != null? spectrum.getPeptideIntensity():1.0;
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
  }catch (IOException e) {
      log.error("Could not create mzTab file reader", e);
      throw new IOException("The file doesn't not exists -- " + file);
    }
  }
}
