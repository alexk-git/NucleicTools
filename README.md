# NucleicTools
Tools for working with nucleotides and their sequences.

## Core Classes

### BiologicalSequence (ABC)
Abstract base class defining the interface for all biological sequences:

- `__len__`  length of sequence
- `__getitem__`  indexing and slicing
- `__str__`, `__repr__`  string representation
- `__iter__`  iteration over sequence
- `validate_alphabet()`  abstract method for alphabet validation

### NucleicAcidSequence (inherits from BiologicalSequence)
Base class for DNA and RNA sequences. The `_complement_map` attribute must be defined in subclasses.

- `complement()`  returns complementary sequence (raises `NotImplementedError` if called on base class directly)
- `reverse()`  returns reversed sequence
- `reverse_complement()`  returns reverse complement
- `gc_count()`  calculates GC-content as a fraction (0 to 1)
- `nucl_count(nucl)`  counts occurrences of a specific nucleotide

### DNASequence (inherits from NucleicAcidSequence)
- `transcribe()`  transcribes DNA to RNA (returns `RNASequence`)
- DNA-specific alphabet validation (A, T, G, C)

### RNASequence (inherits from NucleicAcidSequence)
- RNA-specific alphabet validation (A, U, G, C)

### AminoAcidSequence (inherits from BiologicalSequence)
- `molecular_weight()`  calculates protein molecular weight
- Amino acid alphabet validation (20 standard amino acids)

### Usage examples

```python
from main import DNASequence, RNASequence, AminoAcidSequence

# DNA
dna = DNASequence("ATGC")
print(dna.complement())          # TACG
print(dna.reverse_complement())  # GCAT
print(dna.gc_count())            # 0.5
rna = dna.transcribe()           # RNASequence("AUGC")

# RNA
rna = RNASequence("AUGC")
print(rna.complement())          # UACG

# Protein
protein = AminoAcidSequence("MA")
print(protein.molecular_weight())  # 220.28
```

## FASTQ Filtering

The function `filter_fastq` reads a FASTQ file, applies filters, and writes passing sequences to an output file in the `filtered/` directory.

### Arguments

- `input_fastq`  path to input FASTQ file
- `output_fastq`  output file name (saved in `filtered/` directory)
- `gc_bounds`  GC-content interval in percent, default `(0, 100)`. A single number is treated as upper bound.
- `length_bounds`  sequence length interval, default `(0, 2**32)`. A single number is treated as upper bound.
- `quality_threshold`  minimum average read quality (Phred33), default `0`
- `overwrite`  if `True`, overwrite existing output file, default `False`

All bounds are inclusive.

### Usage example

```python
from main import filter_fastq

filter_fastq(
    input_fastq="reads.fastq",
    output_fastq="filtered_reads.fastq",
    gc_bounds=(30, 70),
    length_bounds=(50, 500),
    quality_threshold=30,
    overwrite=True
)
```

## bio_files_processor

### convert_multiline_fasta_to_oneline
Reads a FASTA file where sequences may be split across multiple lines and saves a new FASTA file with each sequence on a single line.

### parse_blast_output
Receives a GBK file as input, extracts genes and their neighbors, and saves protein sequences (translations) to a FASTA file.

### find_genes_with_neighbors
Finds genes of interest and their neighbors in a gene dictionary.

```python
from main import find_genes_with_neighbors

result = find_genes_with_neighbors(genes_dict, ["geneA", "geneB"], n_before=2, n_after=2)
```

## Logging

The project uses Python's `logging` module. Logs are written to `logs/app.log` and to the console. The `filter_fastq` function logs informational messages (start, parameters, results) and errors.

## Error Handling

- `ValueError`  invalid sequence alphabet
- `NotImplementedError`  calling `complement()` on `NucleicAcidSequence` directly
- `FileExistsError`  output file already exists and `overwrite=False`
- `RuntimeError`  errors during FASTQ file processing

## Requirements

- Python 3.11+
- biopython

## License

MIT License  see [LICENSE](LICENSE) file for details.
