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


def is_nucleic_acid(seq: str) -> bool:
    '''
    Check if sequence (seq) is valid DNA or RNA sequence.
    '''
    dna = set("atcg")
    rna = set("aucg")

    return set(seq.lower()) <= dna or set(seq.lower()) <= rna:


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
    if not is_dna(seq):
        raise ValueError(f"{seq} is not DNA, can transcribe only DNA.")
    rez = []
    for n in seq:
        match n:
            case 'a':
                rez.append('a')
            case 'A':
                rez.append('A')
            case 't':
                rez.append('u')
            case 'T':
                rez.append('U')
            case 'c':
                rez.append('c')
            case 'C':
                rez.append('C')
            case 'g':
                rez.append('g')
            case 'G':
                rez.append('G')
            case 'u':
                rez.append('t')
            case 'U':
                rez.append('T')
            case _:
                return "not dna or rna sequence"
    return ''.join(rez)


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
    rez = []

    if is_dna(seq):
        for n in seq:
            match n:
                case 'a':
                    rez.append('t')
                case 'A':
                    rez.append('T')
                case 't':
                    rez.append('a')
                case 'T':
                    rez.append('A')
                case 'c':
                    rez.append('g')
                case 'C':
                    rez.append('G')
                case 'g':
                    rez.append('c')
                case 'G':
                    rez.append('C')
        if is_rna(seq):
            for n in seq:
                match n:
                    case 'a':
                        rez.append('u')
                    case 'A':
                        rez.append('U')
                    case 'u':
                        rez.append('a')
                    case 'U':
                        rez.append('A')
                    case 'c':
                        rez.append('g')
                    case 'C':
                        rez.append('G')
                    case 'g':
                        rez.append('c')
                    case 'G':
                        rez.append('C')

    return ''.join(rez)


def reverse_complement(seq: str) -> str:
    '''
    Return reversed complimented sequence.
    '''
    return reverse(complement(seq))


if __name__ == "__main__":
    pass
