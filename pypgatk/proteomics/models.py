import re
from collections import deque
from itertools import chain


class Enzymes:

  def __init__(self):

    self.enzymes = {

      'Trypsin':{
        'name': 'Trypsin',
        'cleavage rule': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
        'accession':'MS:1001251'
      },
      'Arg-C': {
        'name': 'Arg-C',
        'cleavage rule': 'R',
        'accession': 'MS:1001303'
      },
      'Asp-N': {
        'name': 'Asp-N',
        'cleavage rule': '\w(?=D)',
        'accession': 'MS:1001303'
      },
      'Chymotrypsin':{
          'name':'Chymotrypsin',
          'cleavage rule': '([FY](?=[^P]))|(W(?=[^MP]))',
          'accession':'MS:1001306'
        },
    }

  @staticmethod
  def cleave(sequence: str, rule: str, max_missed_cleavages=0, min_acids: int = None, max_acids: int =None):

    """Cleaves a polypeptide sequence using a given rule.
    Taken from pyteomics.parser .

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str
        A string with a regular expression describing the C-terminal site of
        cleavage.
    missed_cleavages : int, optional
        The maximal number of allowed missed cleavages. Defaults to 0.
    Returns
    -------
    out : list
        A list of peptides.
    Examples
    --------
    >>> cleave('AKAKBK', 'KR', 0)
    ['AK', 'AK', 'BK']
    >>> cleave('AKAKBKCK', expasy_rules['trypsin'], 2)
    ['AK', 'AKAK', 'AK', 'AKAKBK', 'AKBK', 'BK', 'AKBKCK', 'BKCK', 'CK']
    """
    peptides = []
    cleavage_sites = deque([0], maxlen=max_missed_cleavages + 2)
    for i in chain([x.end() for x in re.finditer(rule, sequence)],[None]):
      cleavage_sites.append(i)
      for j in range(0, len(cleavage_sites) - 1):
        peptide = sequence[cleavage_sites[j]:cleavage_sites[-1]]
        if peptide == '':
          continue
        if peptide.find('X') != -1:
          continue
        if min_acids is not None and len(peptide) < min_acids:
          continue
        if max_acids is not None and len(peptide) > max_acids:
          continue
        peptides.append(peptide)
    return peptides

PYGPATK_ENZYMES = Enzymes()



#         'asp-n': '\w(?=D)',
#         'bnps-skatole': 'W',
#         'caspase 1': '(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
#         'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
#         'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
#         'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
#         'caspase 5': '(?<=[LW]EH)D',
#         'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
#         'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
#         'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
#         'caspase 9': '(?<=LEH)D',
#         'caspase 10': '(?<=IEA)D',
#         'chymotrypsin low specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
#         'chymotrypsin high specificity':
#           '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
#         'clostripain': 'R',
#         'cnbr': 'M',
#         'enterokinase': '(?<=[DN][DN][DN])K',
#         'factor xa': '(?<=[AFGILTVM][DE]G)R',
#         'formic acid': 'D',
#         'glutamyl endopeptidase': 'E',
#         'granzyme b': '(?<=IEP)D',
#         'hydroxylamine': 'N(?=G)',
#         'iodosobezoic acid': 'W',
#         'lysc': 'K',
#         'ntcb': '\w(?=C)',
#         'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
#                         '((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
#         'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
#                         '((?<=[^HKR][^P])[FL](?=\w[^P]))',
#         'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
#         'proteinase k': '[AEFILTVWY]',
#         'staphylococcal peptidase i': '(?<=[^E])E',
#         'thermolysin': '[^DE](?=[AFILMV])',
#         'thrombin': '((?<=G)R(?=G))|'
#                     '((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
