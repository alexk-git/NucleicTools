'''
    Functions for reading bioinformatics files.

    convert_multiline_fasta_to_oneline: reads a fasta file supplied as input,
    in which the sequence (DNA/RNA/protein/etc.) can be split into several lines,
    and then saves it into a new fasta file in which each sequence fits one line.

    parse_blast_output: function reads the specified TXT file,
    and for each QUERY request (paragraph "Sequences producing significant alignments:"),
    selects the first row from the Description column.
    The resulting set of proteins is saved to a new file in a single column, sorted alphabetically.
    
'''

from pathlib import Path

def convert_multiline_fasta_to_oneline(input_fasta: str) -> None:
    '''
    A function processes fasta file: reads a fasta file supplied as input_fasta (str),
    where one sequence (DNA/RNA/protein/etc.) can be split into several lines,
    and then saves it into a new fasta file (output_fasta) where each sequence fits one line.

    Input:
        input_fasta: file with the input sequences

    Arguments:
        input_fasta: str

    Returns:
        output_fasta: a file where each sequence fits one line

    Raises:
        exceptions if something went wrong.
    '''

    path_to_write = Path(f"output_{input_fasta}")
    
    with open(input_fasta, "r") as file:
        rez = {}
        line = file.readline()
        while line:
            if line and line[0] == '>':
                key_line = line.strip()
                rez[key_line] = []
                n_line = file.readline().strip()
                while n_line and n_line[0] != '>':
                    if n_line:
                        rez[key_line].append(n_line)
                    n_line = file.readline().strip()
                    line = n_line
            else:
                line = file.readline().strip()
        
        with open(path_to_write, "w") as file_w:
            for k in rez.keys():
                file_w.write(k+'\n')
                file_w.write(''.join(rez[k])+'\n')

def parse_blast_output():
    '''
    A function processes fasta file: reads a fasta file supplied as input_fasta (str),
    where one sequence (DNA/RNA/protein/...) can be split into several lines,
    and then saves it into a new fasta file (output_fasta) where each sequence fits one line.

    Input:
        input_fasta: file with the input sequences

    Arguments:
        input_fasta: str

    Returns:
        output_fasta: a file where each sequence fits one line

    Raises:
        exceptions if something went wrong.
    '''
