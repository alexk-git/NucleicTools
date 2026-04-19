'''
Functions for reading and writing bioinformatics files.

    convert_multiline_fasta_to_oneline: reads a multiline FASTA and saves it
    with each sequence on a single line.

    parse_blast_output: parses a GBK file, extracts genes of interest
    with their neighbors, and saves translations to a FASTA file.

    find_genes_with_neighbors: finds genes and their neighbors in a dictionary.

    write_genes_seq_to_fasta: writes gene translations to a FASTA file.
'''

from pathlib import Path
from typing import Union
import json


def convert_multiline_fasta_to_oneline(input_fasta: str) -> None:
    '''
    Reads a FASTA file where sequences may span multiple lines
    and saves a new FASTA file with each sequence on a single line.

    Arguments:
        input_fasta: path to the input FASTA file

    Returns:
        None. Writes output to output_{input_fasta}.
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
                file_w.write(k + '\n')
                file_w.write(''.join(rez[k]) + '\n')


def find_genes_with_neighbors(genes_dict: dict, genes: Union[int, tuple, list], n_before: int = 1, n_after: int = 1) -> dict:
    """
    Finds genes of interest and their neighbors in a dictionary.

    Arguments:
        genes_dict: dictionary with genes {n: {'gene': name, 'translation': seq}}
        genes: gene numbers to search (int, tuple, or list)
        n_before: how many genes before the target gene to include
        n_after: how many genes after the target gene to include

    Returns:
        dict: dictionary with found genes and their neighbors
    """
    rez = {}

    if isinstance(genes, int):
        target_genes = [genes]
    elif isinstance(genes, tuple):
        target_genes = list(genes)
    else:
        target_genes = genes

    genes_to_include = set()

    for gene_num in target_genes:
        if gene_num in genes_dict:
            start = max(1, gene_num - n_before)
            end = gene_num + n_after

            for i in range(start, end + 1):
                if i in genes_dict:
                    genes_to_include.add(i)

    for gene_num in sorted(genes_to_include):
        rez[gene_num] = genes_dict[gene_num]

    return rez


def write_genes_seq_to_fasta(genes_data: dict, output_file: str) -> None:
    """
    Writes gene translations to a FASTA file.

    Arguments:
        genes_data: dictionary with gene data
        output_file: path to the output FASTA file

    Returns:
        None. Writes output to the specified file.
    """

    with open(output_file, 'w', encoding='utf-8') as file:
        for gene_num, gene_info in genes_data.items():
            header = f">gene_{gene_num}"
            if gene_info['gene']:
                header += f"|name_{gene_info['gene']}"

            file.write(header + "\n")
            file.write(gene_info['translation'] + "\n\n")


def parse_blast_output(input_gbk: str, genes: Union[int, tuple, list], output_fasta: str, n_before: int = 1, n_after: int = 1) -> None:
    '''
    Parses a GBK file, extracts genes of interest with their neighbors,
    and saves their protein sequences (translations) to a FASTA file.

    Arguments:
        input_gbk: path to the input GBK file
        genes: genes of interest (int, tuple, or list)
        output_fasta: output FASTA file name
        n_before: how many genes before each target to include
        n_after: how many genes after each target to include

    Returns:
        None. Writes output FASTA and intermediate JSON files.
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
                    genes_parsed[current_gene] = {
                        'gene_count': gene_count,
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
                        genes_parsed[current_gene] = {
                            'gene_count': gene_count,
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
        genes_parsed[current_gene] = {
            'gene_count': gene_count,
            'gene': current_gene,
            'translation': current_translation
        }

    with open(f'{input_gbk.split(".")[0]}.json', 'w', encoding='utf-8') as f:
        json.dump(genes_parsed, f, ensure_ascii=False, indent=4)

    genes_of_interests = find_genes_with_neighbors(genes_parsed, genes, n_before, n_after)

    write_genes_seq_to_fasta(genes_of_interests, path_to_write)
