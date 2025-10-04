# NucleicTools
Tools for working with nucleotides and their sequences.

The aggregation function `run_dna_rna_tools` accepts an arbitrary number of arguments containing DNA or RNA sequences (str), as well as the name of the procedure to be executed (this is always the last argument; see the usage example). It then performs the specified operation on all the passed sequences and returns the result.

The following procedres are available:

`is_nucleic_acid` - returns a Boolean result of a sequence check

`is_rna` - returns a Boolean result of a RNA-sequence check

`is_dna` - returns a Boolean result of a DNA-sequence check

`transcribe` - returns the transcribed sequence

`reverse` - returns the reversed sequence

`complement` - returns the complementary sequence

`reverse_complement` - returns the reverse complementary sequence

## run_dna_rna_tools usage example
```
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False !!
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

The filtration function `filter_fastq` accepts four arguments: seqs, gc_bounds, length_bounds, quality_threshold:

`seqs` is a dictionary of fastq sequences: key is a string name of the sequence, the value is a tuple of two strings: the sequence and the quality in phred33 scale.

`gc_bounds` - the GC-content interval (in percent), by default is (0, 100). If a single number is passed as an argument, it is considered to be the upper limit.

`length_bounds` - the length interval, by default is (0, 2**32).

`quality_threshold` - the average read quality threshold for filtering; by default, it is 0 (phred33 scale).

returns a dictionary similar to the input one, but consisting only of those sequences that satisfy all filtering conditions.

## filter_fastq usage example

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

## Requirements:
tested on Python 3.11.13


