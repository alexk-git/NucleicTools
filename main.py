from typing import Union
from os import path
from sys import exit
from pathlib import Path
import modules.fastq_tools


from abc import ABC, abstractmethod

'''
    Specialized classes for processing DNA/RNA sequences
    Each function performs one processing of one sequence, the sequence is a string.

    Functions:
        nucl_count: count nucleotides in sequence
        gc_count: count GC-content in sequenct
        average_quality: count average quality of sequence from the sequence quality string

    Raises:
        ValueError: if wrong sequence

'''

class BiologicalSequence(ABC):

     def __init__(self, sequence):
        self._sequence = str(sequence).upper()
        self._validate_alphabet()
     
     @abstractmethod
     def _validate_alphabet(self):
         pass

     def __len__(self):
        return len(self._sequence)

     def __getitem__(self, index):
        return self._sequence[index]

     def __str__(self):
        return self._sequence


     def __repr__(self):
         return f"{self.__class__.__name__}('{self._sequence}')"


     def __iter__(self):
         return iter(self._sequence)


class NucleicAcidSequence(BiologicalSequence):

    _compliments_dna = {"A": "T", "a": "t", "T": "A", "t": "a", "G": "C", "g": "c", "C": "G", "c": "g"}
    _compliments_rna = {"A": "U", "a": "u", "U": "A", "u": "a", "G": "C", "g": "c", "C": "G", "c": "g"}

    def __init__(self, seq):
        super().__init__(seq)


    def _validate_alphabet(self):
        '''
        Check if sequence (seq) is valid DNA or RNA sequence.
        '''
        dna = set("atcg")
        rna = set("aucg")

        if not (set(self._sequence.lower()) <= dna or set(self._sequence.lower()) <= rna):
            raise ValueError(f"{self._sequence} is not nucleic sequence.")


    def complement(self) -> NucleicAcidSequence:
        '''
        Return complimented sequence.
        '''
        dna = set("atcg")
        rna = set("aucg")

        if set(self._sequence.lower()) <= dna:
            return self.__class__("".join([self._compliments_dna[n] for n in self._sequence]))
        else:
            return self.__class__("".join([self._compliments_rna[n] for n in self._sequence]))


    def reverse(self) -> NucleicAcidSequence:
        '''
        Return reversed sequence.
        '''
        return self.__class__(self._sequence[::-1])


    def reverse_complement(self) -> NucleicAcidSequence:
        '''
        Return reversed complimented sequence.
        '''
        return self.reverse().complement()


    def nucl_count(self, nucl: str) -> int:
        '''
        Returns the number (amount) of nucleotide (nucl) in a sequence (posl)
        '''

        posl = self._sequence.lower()
        nucl = nucl.lower()
        nc = 0
        for i in range(len(posl)):
            if posl[i] == nucl:
                nc += 1

        return nc


    def gc_count(self) -> float:
        '''
        Return GC-content (in %) of the read

        Arguments:
            read: string of nucleotides

        Returns:
            float number: percentage of guanine and cetosine to read length

        Raises:
            exceptions if something went wrong
        '''

        g_count = self.nucl_count('G')
        c_count = self.nucl_count('C')

        return (g_count + c_count)*100/len(self._sequence)


class RNASequence(NucleicAcidSequence):

    def __init__(self, seq):
        super().__init__(seq)


    def _validate_alphabet(self):
        '''
        Check if sequence is valid RNA sequence.
        '''
        rna = set("aucg")
        if not set(self._sequence.lower()) <= rna:
            raise ValueError(f"{self._sequence} is not RNA.")


class DNASequence(NucleicAcidSequence):

    _transcribe_map = {"A": "A", "a": "a", "T": "U", "t": "u", "G": "G", "g": "g", "C": "C", "c": "c"}

    def __init__(self, seq):
        super().__init__(seq)


    def _validate_alphabet(self):
        '''
        Check if sequence is valid DNA sequence.
        '''
        dna = set("atcg")
        if not set(self._sequence.lower()) <= dna:
            raise ValueError(f"{self._sequence} is not DNA.")


    def transcribe(self) -> RNASequence:
        '''
        Return transcribed RNA from DNA sequence.
        '''
        return RNASequence("".join([self._transcribe_map[n] for n in self._sequence]))


class AminoAcidSequence(BiologicalSequence):

    def __init__(self, seq):
        super().__init__(seq)


    def _validate_alphabet(self):
        '''
        Check if sequence is valid amino-acid sequence.
        '''
        aa = set("ACDEFGHIKLMNPQRSTVWY")
        if not set(self._sequence.upper()) <= aa:
            raise ValueError(f"{self._sequence} is not amino-acid.")


    def molecular_weight(self) -> float:
        '''
        Mean molecular weight.
        '''

        aa_weights = {
            'A': 89.09, 'C': 121.15, 'D': 133.10, 'E': 147.13,
            'F': 165.19, 'G': 75.07, 'H': 155.16, 'I': 131.17,
            'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
            'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
            'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
        }

        weight = sum(aa_weights.get(aa, 0) for aa in self._sequence.upper())
        weight -= (len(self) - 1) * 18.02

        return weight




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

def write_genes_seq_to_fasta(genes_data: dict, output_file: str):
    """
    Writes genes to a FASTA file.

    Arguments:
        genes_data: dictionary with gene data
        output_file: name of the output FASTA file

    Returns:
        None

    """

    with open(output_file, 'w', encoding='utf-8') as file:
        for gene_num, gene_info in genes_data.items():
            header = f">gene_{gene_num}"
            if gene_info['gene']:
                header += f"|name_{gene_info['gene']}"

            file.write(header + "\n")
            file.write(gene_info['translation'] + "\n\n")

    return None

def find_genes_with_neighbors(genes_all: dict, genes: Union[int, tuple, list], n_before: int = 1, n_after: int = 1) -> dict:
    """
    Finds genes of interest and their neighbors in a dictionary.

    Arguments:
        genes_dict: dictionary with genes {'gene_name': {'gene': gene_name, 'gene_count': number, 'translation': seq}}
        genes: gene numbers to search (int, tuple, or list)
        n_before: how many genes before the target gene to include
        n_after: how many genes after the target gene to include

    Returns:
        dict: dictionary with found genes and their neighbors

    """

    if isinstance(genes, str):
        genes = [genes]

    rez = {}
    genes_to_check = list(genes_all.keys())

    for gene in genes:
        idx = genes_to_check.index(gene)

        left_end = max(idx - n_before, 0)
        right_end = min(idx + n_after + 1, len(genes_to_check))

        for i in range(left_end, right_end):
            key = genes_to_check[i]
            rez[key] = genes_all[key]

    return rez


def filter_fastq(input_fastq: str, output_fastq: str, gc_bounds: Union[int, tuple] = (0, 100), length_bounds: Union[int, tuple] = (0, 2**32), quality_threshold: int = 0) -> None:
    '''
    A function working with fastq sequences.
    All bounds is included.
    Quality in Phred33.

    Input:
        input_fastq: file with the input sequences
        output_fastq: file to store filtered sequences a special directory "filtered" will be created for it

    Arguments:
        gc_bounds: tuple = (0, 100) # bound included
        length_bounds: tuple = (0, 2**32) # bound included
        quality_threshold: int = 0 # in phred33

    Intermediate:
        seqs: dict of a fastq sequences key: sequence_name string, value: tupple: (sequence: str, quality: str)

    Returns:
        dictionary consisting only of sequences that satisfy all conditions.

    Raises:
        exceptions if something went wrong.
    '''

    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    path_to_write = Path("filtered", f"{output_fastq}")
    if path_to_write.is_file():
        is_overwtire = input(f"The file {path_to_write} already exists, want to overwrite it? Y/N")
        if is_overwtire in {'Y', 'y'}: path_to_write.unlink()
        else: exit()
            
    status = True
    
    with open(input_fastq, "r") as file:
        while status:
            seqs = modules.fastq_tools.read_seq_from_file(file)
            if not status:
                print(f"processing of the {input_fastq} is complete, filtering results are saved in {output_fastq}")
                break
            
            if len(seqs.keys()) == 0:
                print(f"processing of the {input_fastq} is complete, filtering results are saved in {output_fastq}")
                break
            

            if not isinstance(seqs, dict): raise TypeError("seqs must be a dictionary")
            if isinstance(gc_bounds, int): gc_bounds = (0, gc_bounds)
            if isinstance(length_bounds, int): length_bounds = (0, length_bounds)
        
            rez_gc_bounds = {}
        
            for key in seqs.keys():
                if (gc_bounds[0] <= modules.fastq_tools.gc_count(seqs[key][0])) and (modules.fastq_tools.gc_count(seqs[key][0]) <= gc_bounds[1]):
                    rez_gc_bounds[key] = seqs[key]
        
            rez_length_bounds = {}
        
            for key in rez_gc_bounds.keys():
                if (length_bounds[0] <= len(seqs[key][0])) and (len(seqs[key][0]) <= length_bounds[1]):
                    rez_length_bounds[key] = rez_gc_bounds[key]
        
            rez_quality_threshold = {}
        
            for key in rez_length_bounds.keys():
                if (modules.fastq_tools.average_quality(seqs[key]) >= quality_threshold):
                    rez_quality_threshold[key] = rez_length_bounds[key]

            if len(rez_quality_threshold.keys())>0: modules.fastq_tools.write_seq_to_fle(path_to_write, rez_quality_threshold)
            else: continue
    
    return None


if __name__ == "__main__":
    pass
