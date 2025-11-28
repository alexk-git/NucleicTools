'''
    Specialized functions for processing DNA/RNA sequences
    Each function performs one processing of one sequence (str).

    Functions:
        is_nucleic_acid: Validate if sequence contains only valid nucleotides
        is_dna: Validate if sequence contains only DNA-valid nucleotides
        is_rna: Validate if sequence contains only RNA-valid nucleotides
        transcribe: Transcribe DNA to RNA (replace T with U)
        reverse: Reverse the sequence order
        complement: Return complementary sequence (A<->T, A<->U, C<->G)
        reverse_complement: Return reverse complementary sequence

    Example:
        is_nucleic_acid('TTUU') # False
        is_dna('AUG') # False
        is_rna('UGCAU') # True
        transcribe('ATG') # 'AUG'
        reverse('ATG') # 'GTA'
        complement('AtG') # 'TaC'
        reverse_complement('ATg') # 'cAT'
        reverse('ATG') # 'GTA'

    Raises:
        ValueError: if wrong sequence
'''

compliments_dna = {"A": "T", "a": "t", "T": "A", "t": "a", "G": "C", "g": "c", "C": "G", "c": "g"}
compliments_rna = {"A": "U", "a": "u", "U": "A", "u": "a", "G": "C", "g": "c", "C": "G", "c": "g"}
trinscribets = {"A": "A", "a": "a", "T": "U", "t": "u", "G": "G", "g": "g", "C": "C", "c": "c"}

def is_nucleic_acid(seq: str) -> bool:
    '''
    Check if sequence (seq) is valid DNA or RNA sequence.
    '''
    dna = set("atcg")
    rna = set("aucg")

    return set(seq.lower()) <= dna or set(seq.lower()) <= rna


def is_dna(seq: str) -> bool:
    '''
    Check if sequence (seq) is valid DNA sequence.
    '''
    dna = set("atcg")
    return set(seq.lower()) <= dna


def is_rna(seq: str) -> bool:
    '''
    Check if sequence (seq) is valid RNA sequence.
    '''
    rna = set("aucg")
    return set(seq.lower()) <= rna


def transcribe(seq: str) -> str:
    '''
    Return transcribed RNA from DNA sequence.
    '''
    if not is_dna(seq): raise ValueError(f"{seq} is not DNA, can transcribe only DNA.")
    else: return "".join([trinscribets[n] for n in seq])


def reverse(seq: str) -> str:
    '''
    Return reversed sequence.
    '''
    return seq[::-1]


def complement(seq: str) -> str:
    '''
    Return complimented sequence.
    '''
    if not is_nucleic_acid(seq):
        raise ValueError(f"{seq} is not DNA or RNA sqeuence!")

    if is_dna(seq): return "".join([compliments_dna[n] for n in seq])

    if is_rna(seq): return "".join([compliments_rna[n] for n in seq])


def reverse_complement(seq: str) -> str:
    '''
    Return reversed complimented sequence.
    '''
    return reverse(complement(seq))


if __name__ == "__main__":
    pass
