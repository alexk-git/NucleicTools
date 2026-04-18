# NucleicTools
Tools for working with nucleotides and their sequences.

# Modules

## Core Classes

### BiologicalSequence (ABC)
Abstract base class defining the interface for all biological sequences:

* `__len__` - length of sequence
* `__getitem__` - indexing and slicing
* `__str__`, `__repr__` - pretty printing
* `_validate_alphabet` - abstract method for alphabet validation

### NucleicAcidSequence (inherits from BiologicalSequence)
Base class for DNA and RNA sequences with:

* `complement()` - returns complementary sequence
* `reverse()` - returns reversed sequence
* `reverse_complement()` - returns reverse complement
* `gc_content()` - calculates GC-content in percent
* `nucl_count()` - counts specific nucleotides
* Automatic DNA/RNA type detection

### DNASequence (inherits from NucleicAcidSequence)
* `transcribe()` - transcribes DNA to RNA (returns RNASequence object)
* DNA-specific alphabet validation

### RNASequence (inherits from NucleicAcidSequence)
* RNA-specific alphabet validation

### AminoAcidSequence (inherits from BiologicalSequence)
* `molecular_weight()` - calculates protein molecular weight
* Amino acid alphabet validation

## DNA/RNA Tools
The aggregation function `run_dna_rna_tools` processes DNA or RNA sequences with various operations. It accepts an arbitrary number of arguments containing DNA or RNA sequences (str), as well as the name of the procedure to be executed (this is always the last argument), performs the specified operation on all the passed sequences and returns the result.

Available procedures:

* `is_nucleic_acid` - checks if sequence is DNA or RNA
* `is_rna` - checks if sequence is RNA
* `is_dna` - checks if sequence is DNA
* `transcribe` - transcribes DNA to RNA
* `reverse` - reverses sequence
* `complement` - returns complementary sequence
* `reverse_complement` - returns reverse complement

### Usage examples
```python
# Class-based approach
dna = DNASequence("ATG")
rna = dna.transcribe()  # Returns RNASequence("AUG")
print(dna.complement())  # TAC
print(dna.gc_content())  # 33.33

# Function-based approach
run_dna_rna_tools('ATG', 'transcribe')  # 'AUG'
run_dna_rna_tools('ATG', 'aT', 'reverse')  # ['GTA', 'Ta']
```

#### run_dna_rna_tools usage example
```
filter_fastq(
input_fastq="reads.fastq",
output_fastq="filtered_reads.fastq",
gc_bounds=(30, 70),        # 30-70% GC content
length_bounds=(50, 500),   # 50-500 bp length
quality_threshold=30       # Minimum average quality Phred30
)
```

### FASTQ Filtering

The filtration function `filter_fastq` accepts four arguments (seqs, gc_bounds, length_bounds, quality_threshold) and returns a dictionary similar to the input one, but consisting only of those sequences that satisfy all filtering conditions.

`seqs` is a dictionary of fastq sequences: key is a string name of the sequence, the value is a tuple of two strings: the sequence and the quality in phred33 scale.

`gc_bounds` - the GC-content interval (in percent), by default is (0, 100). If a single number is passed as an argument, it is considered to be the upper limit.

`length_bounds` - the length interval, by default is (0, 2**32).

`quality_threshold` - the average read quality threshold for filtering; by default, it is 0 (phred33 scale).

#### filter_fastq usage example

```
fastq_dict = {
    "read_1": ("ATCGATCGAT", "IIIIIIIIII"),
    "read_2": ("GGCCGGCCGG", "=;@B??@<>@"),
    "read_3": ("TTTTAAAAAA", "HHHHHHHHHH"),
    "read_4": ("ACGTACGTAC", "FFFFFFFFF#"),
    "read_5": ("GGGGGGGGGG", "GFFCFEEEFF"),
    "read_6": ("ATATATATAT", "EEEEEEEEEE"),
    "read_7": ("CGCGCGCGCG", "GGGGGGGGGG"),
    "read_8": ("AAAAAAAAAA", "DDDDDDDDDD"),
    "read_9": ("TAGCTAGCGA", "IIIIIIIIII"),
    "read_10": ("CCCCCCCCCC", "D@EDEFFB=D")
}

filter_fastq(fastq_dict, 100, 9, 37)
```

### bio_files_processor

#### convert_multiline_fasta_to_oneline: 
    reads a fasta file supplied as input,
    in which the sequence (DNA/RNA/protein/etc.) can be split into several lines,
    and then saves it into a new fasta file in which each sequence fits one line.

#### select_genes_from_gbk_to_fasta:
    function receives a GBK-file as input, 
    extracts the specified number of genes before and after each gene of interest (gene), 
    and saves their protein sequence (translation) to a fasta file.

### Error Handling
Functions raise appropriate exceptions for invalid inputs.

## Requirements:
Python 3.11+

biopython

pathlib

## License
MIT License - see [LICENSE](LICENSE) file for details.
