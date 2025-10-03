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


def is_nucleic_acid(posl: str) -> bool:
    '''
    Check if posl is valid DNA or RNA sequence.
    '''
    dna = set("atcg")
    rna = set("aucg")

    if set(posl.lower()) <= dna or set(posl.lower()) <= rna:
        return True
    else:
        return False


def is_dna(posl: str) -> bool:
    '''
    Check if posl is valid DNA sequence.
    '''
    dna = set("atcg")
    if set(posl.lower()) <= dna:
        return True
    else:
        return False


def is_rna(posl: str) -> bool:
    '''
    Check if posl is valid RNA sequence.
    '''
    rna = set("aucg")
    if set(posl.lower()) <= rna:
        return True
    else:
        return False


def transcribe(posl: str) -> str:
    '''
    Return transcribed RNA from DNA sequence.
    '''
    if not is_dna(posl):
        raise ValueError(f"{posl} is not DNA, can transcribe only DNA.")
    rez = []
    for i in posl:
        match i:
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


def reverse(posl: str) -> str:
    '''
    Return reversed sequence.
    '''
    return posl[::-1]


def complement(posl: str) -> str:
    '''
    Return complimented sequence.
    '''
    if not is_nucleic_acid(posl):
        raise ValueError(f"{posl} is not DNA or RNA sqeuence!")
    rez = []

    if is_dna(posl):
        for i in posl:
            match i:
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
        if is_rna(posl):
            for i in posl:
                match i:
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


def reverse_complement(posl: str) -> str:
    '''
    Return reversed complimented sequence.
    '''
    return reverse(complement(posl))


if __name__ == "__main__":
    pass
