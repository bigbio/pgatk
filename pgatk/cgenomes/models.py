
class SNP(object):

    def __init__(self, gene=None, mrna=None, dna_mut=None, aa_mut=None, mutation_type=None):
        """
        Default constructor of SNP
        :param gene:
        :param dna_mut:
        :param aa_mut:
        :param mutation_type:
        """
        self.gene = gene
        self.mrna = mrna
        self.aa_mut = aa_mut
        self.type = mutation_type
        self.dna_mut = dna_mut
