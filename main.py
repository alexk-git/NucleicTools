from typing import Union
from os import path
from pathlib import Path
from typing import Union, Optional
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from abc import ABC, abstractmethod
from logger_config import setup_logger

logger = setup_logger(__name__)

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


    def complement(self) -> 'NucleicAcidSequence':
        '''
        Return complimented sequence.
        '''
        dna = set("atcg")
        rna = set("aucg")

        if set(self._sequence.lower()) <= dna:
            return self.__class__("".join([self._compliments_dna[n] for n in self._sequence]))
        else:
            return self.__class__("".join([self._compliments_rna[n] for n in self._sequence]))


    def reverse(self) -> 'NucleicAcidSequence':
        '''
        Return reversed sequence.
        '''
        return self.__class__(self._sequence[::-1])


    def reverse_complement(self) -> 'NucleicAcidSequence':
        '''
        Return reversed complimented sequence.
        '''
        return self.reverse().complement()


    def nucl_count(self, nucl: str) -> int:
        '''
        Returns the number (amount) of nucleotide (nucl) in a sequence (posl)
        '''

        return self._sequence.lower().count(nucl.lower())


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

        return gc_fraction(self._sequence)


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


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, tuple] = (0, 100),
    length_bounds: Union[int, tuple] = (0, 2**32),
    quality_threshold: int = 0
) -> None:
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

    # bound preprocessing
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    # converting percentages to fractions
    # because gc_fraction (returns 0-1)
    gc_min, gc_max = gc_bounds[0] / 100.0, gc_bounds[1] / 100.0

    # creating output directory
    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / output_fastq

    # checking if file exists
    if output_path.exists():
        response = input(f"File {output_path} exists. Overwrite? (Y/N): ")
        if response.lower() != 'y':
            print("Operation cancelled.")
            return
        output_path.unlink()

    # list of records that have passed the filter for answer
    passed_records = []

    try:
        # reading FASTQ-file with SeqIO
        for record in SeqIO.parse(input_fastq, "fastq"):
            sequence = str(record.seq)
            length = len(sequence)

            # length filter
            if not (length_bounds[0] <= length <= length_bounds[1]):
                continue

            # gc-count filter
            gc_content = gc_fraction(sequence)
            if not (gc_min <= gc_content <= gc_max):
                continue

            # quality filter
            if 'phred_quality' in record.letter_annotations:
                avg_quality = sum(record.letter_annotations['phred_quality']) / length
                if avg_quality < quality_threshold:
                    continue
            else:
                # check threshold > 0 if no quality info
                if quality_threshold > 0:
                    continue

            # add to answer
            passed_records.append(record)

        # write filtered to file
        if passed_records:
            count = SeqIO.write(passed_records, output_path, "fastq")
            print(f"Filtered {count} sequences saved to {output_path}")
        else:
            print("No sequences passed the filters.")
            output_path.touch()  # creating empty file

    except Exception as e:
        raise RuntimeError(f"Error processing FASTQ file: {e}")


if __name__ == "__main__":
    pass
