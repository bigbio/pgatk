from dataclasses import dataclass
from typing import Optional


@dataclass
class SNP:
    """Represents a single nucleotide polymorphism or mutation."""
    gene: Optional[str] = None
    mrna: Optional[str] = None
    dna_mut: Optional[str] = None
    aa_mut: Optional[str] = None
    mutation_type: Optional[str] = None
