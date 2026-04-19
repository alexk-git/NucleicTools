from typing import Union
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC, abstractmethod
from logger_config import setup_logger

logger = setup_logger(__name__)

'''
    Specialized classes for processing DNA/RNA sequences
    Each function performs one processing of one sequence, the sequence is a string.

    Functions:
        nucl_count: count nucleotides in sequence
        gc_count: count GC-content in sequence
        average_quality: count average quality of sequence from the sequence quality string

    Raises:
        ValueError: if wrong sequence

'''

class BiologicalSequence(ABC):

     def __init__(self, sequence):
        self._sequence = str(sequence).upper()
        self.validate_alphabet()

     @abstractmethod
     def validate_alphabet(self):
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

    _complement_map = {}

    def __init__(self, seq):
        super().__init__(seq)

    def validate_alphabet(self):
        '''
        Check if sequence (seq) is valid DNA or RNA sequence.
        '''
        dna = set("ATCG")
        rna = set("AUCG")

        if not (set(self._sequence) <= dna or set(self._sequence) <= rna):
            raise ValueError(f"{self._sequence} is not nucleic sequence.")

    def complement(self) -> 'NucleicAcidSequence':
        '''
        Return complemented sequence.
        '''
        if not self._complement_map:
            raise NotImplementedError("complement is not defined for base NucleicAcidSequence")

        return self.__class__("".join([self._complement_map[n] for n in self._sequence]))

    def reverse(self) -> 'NucleicAcidSequence':
        '''
        Return reversed sequence.
        '''
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        '''
        Return reversed complemented sequence.
        '''
        return self.reverse().complement()

    def nucl_count(self, nucl: str) -> int:
        '''
        Returns the number (amount) of nucleotide (nucl) in a sequence.
        '''
        return self._sequence.lower().count(nucl.lower())

    def gc_count(self) -> float:
        '''
        Return GC-content of the sequence as a fraction (0 to 1).

        Returns:
            float: fraction of guanine and cytosine in the sequence

        Raises:
            exceptions if something went wrong
        '''
        return gc_fraction(self._sequence)


class RNASequence(NucleicAcidSequence):

    _complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def __init__(self, seq):
        super().__init__(seq)

    def validate_alphabet(self):
        '''
        Check if sequence is valid RNA sequence.
        '''
        rna = set("AUCG")
        if not set(self._sequence) <= rna:
            raise ValueError(f"{self._sequence} is not RNA.")


class DNASequence(NucleicAcidSequence):

    _complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    _transcribe_map = {"A": "A", "T": "U", "G": "G", "C": "C"}

    def __init__(self, seq):
        super().__init__(seq)

    def validate_alphabet(self):
        '''
        Check if sequence is valid DNA sequence.
        '''
        dna = set("ATCG")
        if not set(self._sequence) <= dna:
            raise ValueError(f"{self._sequence} is not DNA.")

    def transcribe(self) -> RNASequence:
        '''
        Return transcribed RNA from DNA sequence.
        '''
        return RNASequence("".join([self._transcribe_map[n] for n in self._sequence]))


class AminoAcidSequence(BiologicalSequence):

    def __init__(self, seq):
        super().__init__(seq)

    def validate_alphabet(self):
        '''
        Check if sequence is valid amino-acid sequence.
        '''
        aa = set("ACDEFGHIKLMNPQRSTVWY")
        if not set(self._sequence) <= aa:
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

        weight = sum(aa_weights.get(aa, 0) for aa in self._sequence)
        weight -= (len(self) - 1) * 18.02

        return weight


def find_genes_with_neighbors(genes_all: dict, genes: Union[int, tuple, list], n_before: int = 1, n_after: int = 1) -> dict:
    """
    Finds genes of interest and their neighbors in a dictionary.

    Arguments:
        genes_all: dictionary with genes {'gene_name': {'gene': gene_name, 'gene_count': number, 'translation': seq}}
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
    quality_threshold: int = 0,
    overwrite: bool = False
) -> None:
    '''
    A function working with fastq sequences.
    All bounds are included.
    Quality in Phred33.

    Arguments:
        input_fastq: file with the input sequences
        output_fastq: file to store filtered sequences, a special directory "filtered" will be created for it
        gc_bounds: tuple = (0, 100), bounds included
        length_bounds: tuple = (0, 2**32), bounds included
        quality_threshold: int = 0, in phred33
        overwrite: bool = False, if True, overwrite existing output file

    Returns:
        None. Writes filtered sequences to the output file.

    Raises:
        RuntimeError: if error occurs during FASTQ processing
        FileExistsError: if output file exists and overwrite is False
    '''

    logger.info(f"Starting filter: input={input_fastq}, output={output_fastq}")
    logger.info(f"Parameters: gc_bounds={gc_bounds}, length_bounds={length_bounds}, quality_threshold={quality_threshold}")

    # bound preprocessing
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    # converting percentages to fractions
    # because gc_fraction returns 0-1
    gc_min, gc_max = gc_bounds[0] / 100.0, gc_bounds[1] / 100.0

    # creating output directory
    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / output_fastq

    # checking if file exists
    if output_path.exists():
        if not overwrite:
            logger.error(f"File {output_path} already exists and overwrite is False")
            raise FileExistsError(f"File {output_path} already exists. Use overwrite=True to overwrite.")
        output_path.unlink()

    # list of records that have passed the filter
    passed_records = []

    try:
        total_records = 0
        # reading FASTQ-file with SeqIO
        for record in SeqIO.parse(input_fastq, "fastq"):
            total_records += 1
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

        logger.info(f"Processed {total_records} records, filtered {len(passed_records)}")

        # write filtered to file
        if passed_records:
            count = SeqIO.write(passed_records, output_path, "fastq")
            print(f"Filtered {count} sequences saved to {output_path}")
        else:
            print("No sequences passed the filters.")
            output_path.touch()  # creating empty file

    except Exception as e:
        logger.error(f"Error processing FASTQ file: {e}")
        raise RuntimeError(f"Error processing FASTQ file: {e}")


if __name__ == "__main__":
    pass
