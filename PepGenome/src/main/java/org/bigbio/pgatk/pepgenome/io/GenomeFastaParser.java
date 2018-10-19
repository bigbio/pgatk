package org.bigbio.pgatk.pepgenome.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;

import org.bigbio.pgatk.pepgenome.common.Chromosome;
import org.bigbio.pgatk.pepgenome.common.constants.GenomeMapper;

public class GenomeFastaParser  implements Serializable {

	private static final long serialVersionUID = 350507011419545734L;

	public static void readGenomeFASTA(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">")) {
					String[] split = line.split(" ");
					if (line.contains("dna:chromosome") || line.contains("dna:genescaffold")) {
						Chromosome.addChr(split[0].substring(1));
					} else if (line.contains("dna:scaffold")) {
						Chromosome.addChr("scaffold");
					}
				}
			}
			Chromosome.addChr("scaffold");
			GenomeMapper.PEPTIDE_MAPPER.CHR_FROM_GENOME_FASTA = true;
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
