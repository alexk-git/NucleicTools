
def run_dna_rna_tools(*seq_data):
    '''
    One function to rule them all!
    '''
    to_do = seq_data[-1]
    rez = []

    match to_do:
        case 'is_nucleic_acid':
            for posl in seq_data[:-1]:
                rez.append(is_nucleic_acid(posl))
        case 'transcribe':
            for posl in seq_data[:-1]:
                rez.append(transcribe(posl))
        case 'reverse':
            for posl in seq_data[:-1]:
                rez.append(reverse(posl))
        case 'complement':
            for posl in seq_data[:-1]:
                rez.append(complement(posl))
        case 'reverse_complement':
            for posl in seq_data[:-1]:
                rez.append(reverse_complement(posl))
        case _:
            return f"Для для метода {to_do} функций нет!"

    if len(rez) == 1:
        return rez[0]
    else:
        return rez

def filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold):
    pass
