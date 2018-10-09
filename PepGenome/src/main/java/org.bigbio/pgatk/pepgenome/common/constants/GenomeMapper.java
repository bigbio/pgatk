package org.bigbio.pgatk.pepgenome.common.constants;

import org.bigbio.pgatk.pepgenome.common.Chromosome;
import org.bigbio.pgatk.pepgenome.common.constants.TaxonomyIdentifiers;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
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
        TAX.put("bos taurus", TaxonomyIdentifiers.COW);
        TAX.put("cow", TaxonomyIdentifiers.COW);
        TAX.put("9913", TaxonomyIdentifiers.COW);
        TAX.put("callithrix jacchus", TaxonomyIdentifiers.MARMOSET);
        TAX.put("marmoset", TaxonomyIdentifiers.MARMOSET);
        TAX.put("9483", TaxonomyIdentifiers.MARMOSET);
        TAX.put("canis lupus familiaris", TaxonomyIdentifiers.DOG);
        TAX.put("dog", TaxonomyIdentifiers.DOG);
        TAX.put("9615", TaxonomyIdentifiers.DOG);
        TAX.put("chlorocebus sabaeus", TaxonomyIdentifiers.VERVERT_AGM);
        TAX.put("vervet_agm", TaxonomyIdentifiers.VERVERT_AGM);
        TAX.put("60711", TaxonomyIdentifiers.VERVERT_AGM);
        TAX.put("ciona intestinalis", TaxonomyIdentifiers.CITESTINALIS);
        TAX.put("cintestinalis", TaxonomyIdentifiers.CITESTINALIS);
        TAX.put("7719", TaxonomyIdentifiers.CITESTINALIS);
        TAX.put("equus caballus", TaxonomyIdentifiers.HORSE);
        TAX.put("horse", TaxonomyIdentifiers.HORSE);
        TAX.put("9796", TaxonomyIdentifiers.HORSE);
        TAX.put("felis catus", TaxonomyIdentifiers.CAT);
        TAX.put("cat", TaxonomyIdentifiers.CAT);
        TAX.put("9685", TaxonomyIdentifiers.CAT);
        TAX.put("gallus gallus", TaxonomyIdentifiers.CHICKEN);
        TAX.put("chicken", TaxonomyIdentifiers.CHICKEN);
        TAX.put("9031", TaxonomyIdentifiers.CHICKEN);
        TAX.put("gorilla gorilla gorilla", TaxonomyIdentifiers.GORILLA);
        TAX.put("gorilla", TaxonomyIdentifiers.GORILLA);
        TAX.put("9595", TaxonomyIdentifiers.GORILLA);
        TAX.put("homo sapiens", TaxonomyIdentifiers.HUMAN);
        TAX.put("human", TaxonomyIdentifiers.HUMAN);
        TAX.put("9606", TaxonomyIdentifiers.HUMAN);
        TAX.put("macaca mulatta", TaxonomyIdentifiers.MACAQUE);
        TAX.put("macaque", TaxonomyIdentifiers.MACAQUE);
        TAX.put("9544", TaxonomyIdentifiers.MACAQUE);
        TAX.put("meleagris gallopavo", TaxonomyIdentifiers.TURKEY);
        TAX.put("turkey", TaxonomyIdentifiers.TURKEY);
        TAX.put("9103", TaxonomyIdentifiers.TURKEY);
        TAX.put("monodelphis domestica", TaxonomyIdentifiers.OPOSUM);
        TAX.put("opossum", TaxonomyIdentifiers.OPOSUM);
        TAX.put("13616", TaxonomyIdentifiers.OPOSUM);
        TAX.put("mus musculus", TaxonomyIdentifiers.MOUSE);
        TAX.put("mouse", TaxonomyIdentifiers.MOUSE);
        TAX.put("10090", TaxonomyIdentifiers.MOUSE);
        TAX.put("ornithorhynchus anatinus", TaxonomyIdentifiers.PLATYPUS);
        TAX.put("platypus", TaxonomyIdentifiers.PLATYPUS);
        TAX.put("9258", TaxonomyIdentifiers.PLATYPUS);
        TAX.put("oryctolagus cuniculus", TaxonomyIdentifiers.RABBIT);
        TAX.put("rabbit", TaxonomyIdentifiers.RABBIT);
        TAX.put("9986", TaxonomyIdentifiers.RABBIT);
        TAX.put("oryzias latipes", TaxonomyIdentifiers.MEDAKA);
        TAX.put("medaka", TaxonomyIdentifiers.MEDAKA);
        TAX.put("8090", TaxonomyIdentifiers.MEDAKA);
        TAX.put("ovis aries", TaxonomyIdentifiers.SHEEP);
        TAX.put("sheep", TaxonomyIdentifiers.SHEEP);
        TAX.put("9940", TaxonomyIdentifiers.SHEEP);
        TAX.put("pan troglodytes", TaxonomyIdentifiers.CHIMP);
        TAX.put("chimpanzee", TaxonomyIdentifiers.CHIMP);
        TAX.put("9598", TaxonomyIdentifiers.CHIMP);
        TAX.put("papio anubis", TaxonomyIdentifiers.OLIVEBABOOM);
        TAX.put("olivebaboon", TaxonomyIdentifiers.OLIVEBABOOM);
        TAX.put("9555", TaxonomyIdentifiers.OLIVEBABOOM);
        TAX.put("pongo abelii", TaxonomyIdentifiers.ORANGUTAN);
        TAX.put("orangutan", TaxonomyIdentifiers.ORANGUTAN);
        TAX.put("9601", TaxonomyIdentifiers.ORANGUTAN);
        TAX.put("rattus norvegicus", TaxonomyIdentifiers.RAT);
        TAX.put("rat", TaxonomyIdentifiers.RAT);
        TAX.put("10116", TaxonomyIdentifiers.RAT);
        TAX.put("sus scrofa", TaxonomyIdentifiers.PIG);
        TAX.put("pig", TaxonomyIdentifiers.PIG);
        TAX.put("9823", TaxonomyIdentifiers.PIG);
        TAX.put("taeniopygia guttata", TaxonomyIdentifiers.ZEBRAFISH);
        TAX.put("zebra finch", TaxonomyIdentifiers.ZEBRAFISH);
        TAX.put("59729", TaxonomyIdentifiers.ZEBRAFISH);
        TAX.put("tetraodon nigroviridid", TaxonomyIdentifiers.TETRAODON);
        TAX.put("tetraodon", TaxonomyIdentifiers.TETRAODON);
        TAX.put("99883", TaxonomyIdentifiers.TETRAODON);
    }

    /**
     * ScaffoldIdentifier
     *
     */
    public static class ScaffoldIdentifier {

        final static List<String> SCAFFOLD_IDENTIFIER_START = Arrays.asList(new String[]{"ACFV","AAEX", "AQIB", "JSUE", "AAGW", "AMGL", "AACZ", "AHZZ","AABR", "AAND"});
        final static List<String> SCAFFOLD_IDENTIFIER_END = Arrays.asList(new String[]{"hap1", "hap2", "ndom"});

        /**
         * Find scaffold for an specific Scafold notation.
         * @param substr
         * @return Scaffold
         */
        public static Chromosome findScaffold(String substr){
            Chromosome scaffold = null;
            for(String value: SCAFFOLD_IDENTIFIER_START){
                if(substr.startsWith(value))
                    scaffold = Chromosome.scaffold;
            }
            for(String value: SCAFFOLD_IDENTIFIER_END){
                if(substr.endsWith(value))
                    scaffold = Chromosome.scaffold;
            }
            return scaffold;

        }
    }
}