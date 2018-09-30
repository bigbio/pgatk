package org.bigbio.pgatk.pepgenome.common;

import java.util.HashMap;
import java.util.Map;

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
         * KMER_LENGTH holds the size of the kmers in the KmerMap.
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
         * toggles wheter 1 in 5 mode is on. this mode only works with
         * two mismatches. if 1 in 5 mode is on, only one mismatch is
         * allowed in every 5 amino acids. (this only works if KmerLength == 5)
         * can be toggled with the -mmmode switch.
         * */
        public static boolean ONE_IN_FIVE_MODE = false;
    }

    //this holds information on the used gene and transcript ids.
    public static class ID {
        //the gene id is the prefix (e. g. ENSG or ENSMUSG) used for the gene identifier.
        public static String GENE_ID = "ENSG";

        //the transcript id is the prefix (e. g. ENST or ENSMUST) used for the transcript identifier.
        public static String TRANSCRIPT_ID = "ENST";

        //the exon id is the prefix (e. g. ENSE or ENSMUSE) used for the exon identifier.
        public static String EXON_ID = "ENSE";

        //length holds the combined length of the prefix (see above) and the number.
        //(a default length of 11 is assumed for the number making for example ensembl numbers 15 characters long. (ENSG+11))
        public static int LENGTH = GENE_ID.length() + 11;
    }

    public static Map<String, TaxonomyIdentifiers> TAX = new HashMap<>();

    static {
        TAX.put("bos taurus", TaxonomyIdentifiers.cow);
        TAX.put("cow", TaxonomyIdentifiers.cow);
        TAX.put("9913", TaxonomyIdentifiers.cow);
        TAX.put("callithrix jacchus", TaxonomyIdentifiers.marmoset);
        TAX.put("marmoset", TaxonomyIdentifiers.marmoset);
        TAX.put("9483", TaxonomyIdentifiers.marmoset);
        TAX.put("canis lupus familiaris", TaxonomyIdentifiers.dog);
        TAX.put("dog", TaxonomyIdentifiers.dog);
        TAX.put("9615", TaxonomyIdentifiers.dog);
        TAX.put("chlorocebus sabaeus", TaxonomyIdentifiers.vervet_AGM);
        TAX.put("vervet_agm", TaxonomyIdentifiers.vervet_AGM);
        TAX.put("60711", TaxonomyIdentifiers.vervet_AGM);
        TAX.put("ciona intestinalis", TaxonomyIdentifiers.cintestinalis);
        TAX.put("cintestinalis", TaxonomyIdentifiers.cintestinalis);
        TAX.put("7719", TaxonomyIdentifiers.cintestinalis);
        TAX.put("equus caballus", TaxonomyIdentifiers.horse);
        TAX.put("horse", TaxonomyIdentifiers.horse);
        TAX.put("9796", TaxonomyIdentifiers.horse);
        TAX.put("felis catus", TaxonomyIdentifiers.cat);
        TAX.put("cat", TaxonomyIdentifiers.cat);
        TAX.put("9685", TaxonomyIdentifiers.cat);
        TAX.put("gallus gallus", TaxonomyIdentifiers.chicken);
        TAX.put("chicken", TaxonomyIdentifiers.chicken);
        TAX.put("9031", TaxonomyIdentifiers.chicken);
        TAX.put("gorilla gorilla gorilla", TaxonomyIdentifiers.gorilla);
        TAX.put("gorilla", TaxonomyIdentifiers.gorilla);
        TAX.put("9595", TaxonomyIdentifiers.gorilla);
        TAX.put("homo sapiens", TaxonomyIdentifiers.human);
        TAX.put("human", TaxonomyIdentifiers.human);
        TAX.put("9606", TaxonomyIdentifiers.human);
        TAX.put("macaca mulatta", TaxonomyIdentifiers.macaque);
        TAX.put("macaque", TaxonomyIdentifiers.macaque);
        TAX.put("9544", TaxonomyIdentifiers.macaque);
        TAX.put("meleagris gallopavo", TaxonomyIdentifiers.turkey);
        TAX.put("turkey", TaxonomyIdentifiers.turkey);
        TAX.put("9103", TaxonomyIdentifiers.turkey);
        TAX.put("monodelphis domestica", TaxonomyIdentifiers.opossum);
        TAX.put("opossum", TaxonomyIdentifiers.opossum);
        TAX.put("13616", TaxonomyIdentifiers.opossum);
        TAX.put("mus musculus", TaxonomyIdentifiers.mouse);
        TAX.put("mouse", TaxonomyIdentifiers.mouse);
        TAX.put("10090", TaxonomyIdentifiers.mouse);
        TAX.put("ornithorhynchus anatinus", TaxonomyIdentifiers.platypus);
        TAX.put("platypus", TaxonomyIdentifiers.platypus);
        TAX.put("9258", TaxonomyIdentifiers.platypus);
        TAX.put("oryctolagus cuniculus", TaxonomyIdentifiers.rabbit);
        TAX.put("rabbit", TaxonomyIdentifiers.rabbit);
        TAX.put("9986", TaxonomyIdentifiers.rabbit);
        TAX.put("oryzias latipes", TaxonomyIdentifiers.medaka);
        TAX.put("medaka", TaxonomyIdentifiers.medaka);
        TAX.put("8090", TaxonomyIdentifiers.medaka);
        TAX.put("ovis aries", TaxonomyIdentifiers.sheep);
        TAX.put("sheep", TaxonomyIdentifiers.sheep);
        TAX.put("9940", TaxonomyIdentifiers.sheep);
        TAX.put("pan troglodytes", TaxonomyIdentifiers.chimp);
        TAX.put("chimpanzee", TaxonomyIdentifiers.chimp);
        TAX.put("9598", TaxonomyIdentifiers.chimp);
        TAX.put("papio anubis", TaxonomyIdentifiers.olivebaboon);
        TAX.put("olivebaboon", TaxonomyIdentifiers.olivebaboon);
        TAX.put("9555", TaxonomyIdentifiers.olivebaboon);
        TAX.put("pongo abelii", TaxonomyIdentifiers.orangutan);
        TAX.put("orangutan", TaxonomyIdentifiers.orangutan);
        TAX.put("9601", TaxonomyIdentifiers.orangutan);
        TAX.put("rattus norvegicus", TaxonomyIdentifiers.rat);
        TAX.put("rat", TaxonomyIdentifiers.rat);
        TAX.put("10116", TaxonomyIdentifiers.rat);
        TAX.put("sus scrofa", TaxonomyIdentifiers.pig);
        TAX.put("pig", TaxonomyIdentifiers.pig);
        TAX.put("9823", TaxonomyIdentifiers.pig);
        TAX.put("taeniopygia guttata", TaxonomyIdentifiers.zebrafinch);
        TAX.put("zebra finch", TaxonomyIdentifiers.zebrafinch);
        TAX.put("59729", TaxonomyIdentifiers.zebrafinch);
        TAX.put("tetraodon nigroviridid", TaxonomyIdentifiers.tetraodon);
        TAX.put("tetraodon", TaxonomyIdentifiers.tetraodon);
        TAX.put("99883", TaxonomyIdentifiers.tetraodon);
    }
}