class Enzymes:

  def __init__(self):
    """Creates the Enzymes model."""
    self.enzymes = {
      'Trypsin': {'cleavage rule': '(?<=[KRX])(?!P)', 'PSIID': 'MS:1001251', 'cleavage sites': 'KR'},
      'Arg-C': {'cleavage rule': '(?<=[RX])(?!P)', 'PSIID': 'MS:1001303', 'cleavage sites': 'R'},
      'Arg-C/P': {'cleavage rule': '(?<=[RX])', 'PSIID': '', 'cleavage sites': 'R'},
      'Asp-N': {'cleavage rule': '(?=[DBX])', 'PSIID': 'MS:1001304', 'cleavage sites': 'DB'},
      'Asp-N/B': {'cleavage rule': '(?=[DX])', 'PSIID': '', 'cleavage sites': 'D'},
      'Asp-N_ambic': {'cleavage rule': '(?=[DBEZX])', 'PSIID': 'MS:1001305', 'cleavage sites': 'DBEZ'},
      'Chymotrypsin': {'cleavage rule': '(?<=[FYWLJX])(?!P)', 'PSIID': 'MS:1001306', 'cleavage sites': 'FYWLJ'},
      'Chymotrypsin/P': {'cleavage rule': '(?<=[FYWLJX])', 'PSIID': '', 'cleavage sites': 'FYWLJ'},
      'CNBr': {'cleavage rule': '(?<=[MX])', 'PSIID': 'MS:1001307', 'cleavage sites': 'M'},
      'Formic_acid': {'cleavage rule': '((?<=[DBX]))|((?=[DBX]))', 'PSIID': 'MS:1001308', 'cleavage sites': 'DB'},
      'Lys-C': {'cleavage rule': '(?<=[KX])(?!P)', 'PSIID': 'MS:1001309', 'cleavage sites': 'K'},
      'Lys-N': {'cleavage rule': '(?=[KX])', 'PSIID': '', 'cleavage sites': 'K'},
      'Lys-C/P': {'cleavage rule': '(?<=[KX])', 'PSIID': 'MS:1001310', 'cleavage sites': 'K'},
      'PepsinA': {'cleavage rule': '(?<=[FLJX])', 'PSIID': 'MS:1001311', 'cleavage sites': 'FLJ'},
      'TrypChymo': {'cleavage rule': '(?<=[FYWLJKRX])(?!P)', 'PSIID': 'MS:1001312', 'cleavage sites': 'FYWLJKR'},
      'Trypsin/P': {'cleavage rule': '(?<=[KRX])', 'PSIID': 'MS:1001313', 'cleavege sites': 'KR'},
      'V8-DE': {'cleavage rule': '(?<=[DBEZX])(?!P)', 'PSIID': 'MS:1001314', 'cleavege sites': 'DBEZX'},
      'V8-E': {'cleavage rule': '(?<=[EZX])(?!P)', 'PSIID': 'MS:1001315', 'cleavage sites': 'EZX'},
      'leukocyte elastase': {'cleavage rule': '(?<=[ALIJVX])(?!P)', 'PSIID': 'MS:1001915', 'cleavage sites': 'ALIJVX'},
      'proline endopeptidase': {'cleavage rule': '(?<=[HKRX][PX])(?!P)', 'PSIID': 'MS:1001916'
        , 'cleavage sites': 'HKRX'},
      'glutamyl endopeptidase': {'cleavage rule': '(?<=[DBEZX])', 'PSIID': 'MS:1001917', 'cleavage sites': 'DBEZX'},
      'alphaLP': {'cleavage rule': '(?<=[TASVX])', 'PSIID': '', 'cleavage sites': 'TASVX'},
      '2-iodobenzoate': {'cleavage rule': '(?<=[WX])', 'PSIID': 'MS:1001918', 'cleavage sites': 'WX'},
      'iodosobenzoate': {'cleavage rule': '(?<=W)', 'PSIID': '', 'cleavage sites': 'W'},
      'staphylococcal protease/D': {'cleavage rule': '(?<=[EZX])', 'PSIID': '', 'cleavage sites': 'EZX'},
      'proline-endopeptidase/HKR': {'cleavage rule': '(?<=[PX])', 'PSIID': '', 'cleavage sites': 'PX'},
      'Glu-C+P': {'cleavage rule': '(?<=[DBEZX])(?!P)', 'PSIID': '', 'cleavage sites': 'DBEZX'},
      'PepsinA + P': {'cleavage rule': '(?<=[FLJX])(?!P)', 'PSIID': '', 'cleavage sites': 'FLJX'},
      'cyanogen-bromide': {'cleavage rule': '(?<=[MX])', 'PSIID': '', 'cleavage sites': 'MX'},
      'Clostripain/P': {'cleavage rule': '(?<=[RX])', 'PSIID': '', 'cleavage sites': 'RX'},
      'elastase-trypsin-chymotrypsin': {'cleavage rule': '(?<=[ALIVKRWFYX])(?!P)', 'PSIID': '',
                                        'cleavage sites': 'ALIVKRWFYX'},
      'unspecific cleavage': {'cleavage rule': '(?<=[A-Z])', 'PSIID': 'MS:1001956', 'cleavage sites': ''}
    }


# Constants

PYGPATK_ENZYMES = Enzymes()
PYGPATK_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBJOUXZ"
