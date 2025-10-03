'''
    Specialized functions for processing DNA/RNA sequences
    Each function performs one processing of one sequence, the sequence is a string.

    Functions:
        nucl_count: count nucleotides in sequence
        gc_count: count GC-content in sequenct
        average_quality: count average quality of sequence from the sequence quality string

    Raises:
        ValueError: if wrong sequence

'''

from . import dna_rna_tools


def nucl_count(posl: str, nucl: str) -> int:
    '''
    Returns the number (amount) of nucleotide (nucl) in a sequence (posl)
    '''

    if dna_rna_tools.is_nucleic_acid(posl):
        posl = posl.lower()
        nucl = nucl.lower()
        nc = 0
        for i in range(len(posl)):
            if posl[i] == nucl:
                nc += 1

    return nc


def gc_count(read: str) -> float:
    '''
    Return GC-content (in %) of the read

    Arguments:
        read: string of nucleotides

    Returns:
        float number: percentage of guanine and cetosine to read length

    Raises:
        exceptions if something went wrong
    '''

    if dna_rna_tools.is_nucleic_acid(read):
        g_count = nucl_count(read, 'G')
        c_count = nucl_count(read, 'C')

        return (g_count + c_count)*100/len(read)


def average_quality(read: tuple) -> int:
    '''
    Return average quality (in phred33) of the read

    Arguments:
        read: string of nucleotides

    Returns:
        int number: quality number in phred33 score

        Symbol ! " # $ % & ' ( ) * +  ,  -  .  /  0  1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ?  @  A  B  C  D  E  F  G  H  I
         Score 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40

    Raises:
        exceptions if something went wrong
    '''

    scores = {'!': '0', '"': '1', '#': '2', '$': '3', '%': '4', '&': '5', "'": '6', '(': '7', ')': '8', '*': '9', '+': '10', ',': '11', '-': '12', '.': '13', '/': '14', '0': '15', '1': '16', '2': '17', '3': '18', '4': '19', '5': '20', '6': '21', '7': '22', '8': '23', '9': '24', ':': '25', ';': '26', '<': '27', '=': '28', '>': '29', '?': '30', '@': '31', 'A': '32', 'B': '33', 'C': '34', 'D': '35', 'E': '36', 'F': '37', 'G': '38', 'H': '39', 'I': '40'}

    seq_score = 0

    for i in range(len(read[0])):
        seq_score += int(scores[read[1][i]])

    return round(seq_score/len(read[0]))


if __name__ == "__main__":
    pass
