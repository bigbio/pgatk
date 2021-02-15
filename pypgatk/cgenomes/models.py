
class SNP(object):

    def __init__(self, gene=None, mrna=None, dna_mut=None, aa_mut=None, type=None):
        """
        Default constructor of SNP
        :param gene:
        :param dna_mut:
        :param aa_mut:
        :param type:
        """
        self.gene = gene
        self.mrna = mrna
        self.aa_mut = aa_mut
        self.type = type
        self.dna_mut = dna_mut
