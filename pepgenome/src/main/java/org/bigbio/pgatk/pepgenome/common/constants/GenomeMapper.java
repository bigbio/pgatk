package org.bigbio.pgatk.pepgenome.common.constants;

/**
 *  These are the default parameters to map the Peptides to a Genome. Some of these parameters
 *  can be change on request by the main tool.
 *
 * @author ypriverol
 *
 */
public class GenomeMapper {

    public static class PEPTIDE_MAPPER {

        /**
         * KMER_LENGTH holds the size of the kmers in the KmerTreeMap.
         * The default value is 5. tests showed that any other kmersize slows down.
         * Significantly if the number of mappings is high.
         */

        public static int KMER_LENGTH = 5;

        /**
         * Allowed mismatches holds the number of allowed mismatches and has to be between
         * 0 and 2. This can be modified with the -mm input parameter.
         */

        public static int ALLOWED_MISMATCHES = 0;

        /** Allowed amino acids holds the aminoa acids that the PossibleKeyGenerator
         * will use to generate Keys.
         **/
        public static char[] ALLOWED_AMINO_ACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

        /**
         * toggles whether 1 in 5 mode is on. this mode only works with
         * two mismatches. if 1 in 5 mode is on, only one mismatch is
         * allowed in every 5 amino acids. (this only works if KmerLength == 5)
         * can be toggled with the -mmmode switch.
         * */
        public static boolean ONE_IN_FIVE_MODE = false;
        
        /**
         * toggles whether chromosomes and scaffolds are extracted from genome FASTA file and
         * used in output for order of chromosomes and separation of assembly and scaffolds.
         * if genome fasta file is provided through parameter -genome CHR_FROM_GENOME_FASTA==true
         * otherwise chromosome order is extracted from GTF and no separation of assembly and scaffold enabled.
         * */
        public static boolean CHR_FROM_GENOME_FASTA = false;
    }
}