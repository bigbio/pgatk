package org.bigbio.pgatk.pepgenome.io;

import io.github.bigbio.pgatk.io.common.PgatkIOException;
import io.github.bigbio.pgatk.io.pride.*;
import lombok.extern.slf4j.Slf4j;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

@Slf4j
public class PeptideAvroFileParser implements PeptideInputReader, Serializable {

  PrideAvroWriter writer;


  public PeptideAvroFileParser(String fileOutput) {
    try {
      writer = new PrideAvroWriter(new File(fileOutput));
    } catch (PgatkIOException e) {
      e.printStackTrace();
    }
  }

  @Override
  public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

    File prideJsonFile = new File(file);
    FileOutputStream ofs = new FileOutputStream(unmappedoutput);

    if (!prideJsonFile.canRead()) {
      log.error("could not read '" + file + "'.");
      throw new IOException("The file doesn't not exists -- " + file);
    }

    try {

      PrideAvroReader reader = new PrideAvroReader(prideJsonFile);

      String peptideString;
      String isoSeqWithoutPtms;
      int sigPSMs;
      double quant;
      Map<String, TranscriptsT> gene_id_map;



      while(reader.hasNext()){

        AnnotatedSpectrum spectrum = reader.next();

        String tissue = file;

        if(spectrum.getSample() != null && spectrum.getSample().size() > 0){
          for(AvroTuple sampleProperty: spectrum.getSample()){
            if(sampleProperty.getKey().equalsIgnoreCase("organism part")){
              tissue = sampleProperty.getValue();
            }
          }
        }

        peptideString = spectrum.getPepSequence();
        quant = spectrum.getPeptideIntensity() != null? spectrum.getPeptideIntensity():1.0;
        sigPSMs = 1;

        //the matching will only use the amino acids.
        isoSeqWithoutPtms = Utils.make_iso_sequence(Utils.remove_ptms(peptideString));

        if (!coordwrapper.isPeptidePresent(isoSeqWithoutPtms)) {
          //the gene_id_map.find_peptide function will match the peptide.
          gene_id_map = k.find_peptide(isoSeqWithoutPtms);
          for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
            mapping.add_peptide(coordwrapper, peptideString, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it, k.getIsVariant());
          }

          if (gene_id_map.isEmpty()){
            ofs.write(("No-Gene" + "\t" + peptideString + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant +
              "\n").getBytes());
          }

          Optional<PeptideEntry> peptide = coordwrapper.get_existing_peptides_at(peptideString) != null && coordwrapper.get_existing_peptides_at(peptideString).size() >0 ?
            coordwrapper.get_existing_peptides_at(peptideString).stream().findFirst(): Optional.empty();

          if(peptide.isPresent()){
            List<GeneCoordinates> summaryCoords = Utils.coordinates_info(peptide.get(), peptideString, true, true);
            log.info(summaryCoords.toString());
            spectrum.setGeneLocalizations(summaryCoords);
            if(summaryCoords.size() > 0){
              List<String> geneAccessions = spectrum.getGeneAccessions();
              if (geneAccessions == null)
                geneAccessions = new ArrayList<>();
              List<String> finalGeneAccessions = geneAccessions;
              summaryCoords.forEach(x-> {
                finalGeneAccessions.add(x.getGeneAccession());
                finalGeneAccessions.add(x.getTranscriptAccession());
                x.getExonInfoList().forEach( y -> finalGeneAccessions.add(y.getExonAccession()));
              });
              if(geneAccessions.size() > 0){
                HashSet<String> geneSet = new HashSet<>(finalGeneAccessions);
                spectrum.setGeneAccessions(new ArrayList<>(geneSet));
              }
            }
          }

        } else {
          //if the peptide already exists its genomic coordinates dont have to be recalculated.
          //only the tags and PTMs have to be added
          ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(isoSeqWithoutPtms);
          for (PeptideEntry aRefVec : refVec) {
            aRefVec.add_peptide(peptideString, file, sigPSMs, quant, k.getIsVariant());
          }
        }
        writer.write(spectrum);
      }
      writer.close();
      reader.close();
  }catch (IOException e) {
      log.error("Could not create Archive Spectrum file reader", e);
      throw new IOException("The file doesn't not exists -- " + file);
    }
  }
}
