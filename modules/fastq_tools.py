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
from sys import exit

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

def read_seq_from_file(file_name: str) -> dict:
    '''
    Get open file object, read from in a sequence result and make a dictionary of it.
    Update the status of file context if needed.

    Arguments:
        file_name: opened file object for reading

    Global variable:
        status

    Returns:
        seq: dict of a form {id: (seq_read, seq_quality, seq_plus)}
        
    Raises:
        exceptions if something went wrong
    '''
    global status
    seq_id = file_name.readline()
    if not seq_id:
        status = False
        return {}

    else:
        seq_id = seq_id.strip()
        seq_seq = file_name.readline().strip()
        seq_plus = file_name.readline().strip()
        seq_qual = file_name.readline().strip()
        return {seq_id: (seq_seq, seq_qual, seq_plus)}
    

def write_seq_to_fle(file_name: str, seq: dict, ) -> None:
    '''
    Writes given seq (dict) of special shape/form to a given file

    Arguments:
        seq: dict of a form {id: (seq_read, seq_quality)}
        file_name: str of file name to write the dict content

    Returns:
        None
        in the file with file_name four rows to be writed:
        id
        seq_read
        +
        seq_quality
        

    Raises:
        exceptions if something went wrong
    '''
    with file_name.open("a", encoding="utf-8") as file_w:
        for key in seq.keys():
            file_w.write(key+'\n')
            file_w.write(seq[key][0]+'\n')
            file_w.write(seq[key][2]+'\n')
            file_w.write(seq[key][1]+'\n')

    return None

if __name__ == "__main__":
    pass
