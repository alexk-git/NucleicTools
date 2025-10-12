'''
    Functions for reading bioinformatics files.

    convert_multiline_fasta_to_oneline: reads a fasta file supplied as input,
    in which the sequence (DNA/RNA/protein/etc.) can be split into several lines,
    and then saves it into a new fasta file in which each sequence fits one line.

    select_genes_from_gbk_to_fasta: function receives a GBK-file as input, 
    extracts the specified number of genes before and after each gene of interest (gene), 
    and saves their protein sequence (translation) to a fasta file.
    
'''

from pathlib import Path
from typing import Union
import json
import modules.fastq_tools

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

def parse_blast_output(input_gbk: str, genes: Union[int, tuple, list], output_fasta: str, n_before: int = 1, n_after: int = 1):
    '''
    Receives a GBK-file as input, extracts the specified number of genes before and after each gene of interest (gene), 
    and saves their protein sequence (translation) to a fasta file.

    Arguments:
        input_gbk: path to the input GBK file
        genes: genes of interest, near which neighbors are searched (string or collection of strings).
        n_before, n_after: number of genes before and after (>=1)
        output_fasta: output file name.

    Returns:
        output_fasta: a file where each sequence fits one line

    Raises:
        exceptions if something went wrong.
    '''

    genes_parsed = {}
    gene_count = 0
    current_gene = None
    current_translation = ""
    in_translation = False
    in_cds = False
    path_to_write = Path(f"output_{input_gbk.split('.')[0]}")
    
    with open(input_gbk, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.rstrip()
            
            if line.startswith('     CDS'):
                if in_cds and current_gene is not None and current_translation:
                    gene_count += 1
                    genes_parsed[gene_count] = {
                        'gene': current_gene,
                        'translation': current_translation
                    }
                
                in_cds = True
                current_gene = None
                current_translation = ""
                in_translation = False
                continue
            
            if in_cds:
                if line and not line.startswith('                     '):
                    if current_gene is not None and current_translation:
                        gene_count += 1
                        genes_parsed[gene_count] = {
                            'gene': current_gene,
                            'translation': current_translation
                        }
                    
                    in_cds = line.startswith('     CDS')
                    current_gene = None
                    current_translation = ""
                    in_translation = False
                    if in_cds:
                        continue
                
                if '/gene=' in line and current_gene is None:
                    gene_part = line.split('/gene=')[1].strip()
                    if gene_part.startswith('"'):
                        current_gene = gene_part.split('"')[1]
                    else:
                        current_gene = gene_part.split()[0].strip('"')
                
                elif '/translation=' in line:
                    in_translation = True
                    translation_part = line.split('/translation=')[1].strip()
                    if translation_part.startswith('"'):
                        translation_part = translation_part[1:]
                    if translation_part.endswith('"'):
                        translation_part = translation_part[:-1]
                        in_translation = False
                    
                    current_translation = translation_part
                
                elif in_translation:
                    clean_line = line.strip()
                    
                    if '"' in clean_line:
                        clean_line = clean_line.split('"')[0]
                        in_translation = False
                    
                    if clean_line.startswith('/'):
                        in_translation = False
                    else:
                        current_translation += clean_line

    if in_cds and current_gene is not None and current_translation:
        gene_count += 1
        genes_parsed[gene_count] = {
            'gene': current_gene,
            'translation': current_translation
        }

    #print(genes_parsed)

    # For future use: 
    # to avoid new-parsing the if the data is already collected in a convenient form (dictionary)
    with open(f'{input_gbk.split(".")[0]}.json', 'w', encoding='utf-8') as f:
        json.dump(genes_parsed, f, ensure_ascii=False, indent=4)

    genes_of_interests = modules.fastq_tools.find_genes_with_neighbors(genes_parsed, genes, n_before, n_after)

    modules.fastq_tools.write_genes_seq_to_fasta(genes_of_interests, output_fasta)

    return None

    








    