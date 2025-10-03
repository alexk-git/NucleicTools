from typing import Union
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
            return f"Для для метода {to_do} функций нет!"

    if len(rez) == 1:
        return rez[0]
    else:
        return rez


def filter_fastq(seqs: dict, gc_bounds: Union[int, tuple] = (0, 100), length_bounds: Union[int, tuple] = (0, 2**32), quality_threshold: int = 0) -> dict:
    '''
    A function working with fastq sequences.
    All bounds is included.
    Quality in Phred33.

    Input:
        seqs: dict of a fastq sequences key: sequence_name string, value: tupple: (sequence: str, quality: str)

    Arguments:
        gc_bounds: tuple = (0, 100) # bound included
        length_bounds: tuple = (0, 2**32) # bound included
        quality_threshold: int = 0 # in phred33

    Returns:
        dictionary consisting only of sequences that satisfy all conditions.

    Raises:
        exceptions if something went wrong.
    '''


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

    return rez_quality_threshold


if __name__ == "__main__":
    pass
