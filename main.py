from typing import Union
from os import path
from sys import exit
from pathlib import Path
import modules.dna_rna_tools
import modules.fastq_tools


def run_dna_rna_tools(*seq_data):
    '''
    One function to rule them all!

    Receives one or more comma-separated sequences as input and
    desired prodedure, can do the following:
        is_nucleic_acid - returns a Boolean result of a sequence check
        is_rna - returns a Boolean result of a RNA-sequence check
        is_dna - returns a Boolean result of a DNA-sequence check
        transcribe - returns the transcribed sequence
        reverse - returns the reversed sequence
        complement - returns the complementary sequence
        reverse_complement - returns the reverse complementary sequence

    Arguments:
        sequence: one or several, comma-separated
        last argument as always a desired prodecure

    Returns:
        one or list of sequences

    Raises:
        exceptions if something went wrong.
    '''
    to_do = seq_data[-1]
    rez = []

    match to_do:
        case 'is_nucleic_acid':
            for posl in seq_data[:-1]:
                rez.append(modules.dna_rna_tools.is_nucleic_acid(posl))
        case 'transcribe':
            for posl in seq_data[:-1]:
                rez.append(modules.dna_rna_tools.transcribe(posl))
        case 'reverse':
            for posl in seq_data[:-1]:
                rez.append(modules.dna_rna_tools.reverse(posl))
        case 'complement':
            for posl in seq_data[:-1]:
                rez.append(modules.dna_rna_tools.complement(posl))
        case 'reverse_complement':
            for posl in seq_data[:-1]:
                rez.append(modules.dna_rna_tools.reverse_complement(posl))
        case _:
            return f"There are no procerure for the {to_do} method!"

    if len(rez) == 1:
        return rez[0]
    else:
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
