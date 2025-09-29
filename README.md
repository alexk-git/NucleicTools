# NocleicTools
Tools for working with nucleotides and their sequences.

The aggregatin function run_dna_rna_tools accepts an arbitrary number of arguments containing DNA or RNA sequences (str), as well as the name of the procedure to be executed (this is always the last argument; see the usage example). It then performs the specified operation on all the passed sequences and returns the result.

The following procedres are available:

`is_nucleic_acid` - returns a Boolean result of a sequence check

`is_rna` - returns a Boolean result of a RNA-sequence check

`is_dna` - returns a Boolean result of a DNA-sequence check

`transcribe` - returns the transcribed sequence

`reverse` - returns the reversed sequence

`complement` - returns the complementary sequence

`reverse_complement` - returns the reverse complementary sequence


## Usage example
```
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False !!
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```


