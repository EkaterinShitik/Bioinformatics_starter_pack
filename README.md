# Bioinformatics_starter_pack
## Program `bioinf_starter_pack.py`
The program `bioinf_starter_pack.py` provides basic tools to work with three main types of bioinformatic data:
- `dna_rna_tools.py` to process *nucleotide* *sequences*
- `protein_tools.py` to process *amino* *acid* *sequences*
- `fastq_tools.py` to filtrate the obtained sequences in *FASTQ format*


![Image alt](https://github.com/EkaterinShitik/the-second-repository/raw/main/Презентация-без-названия.jpg)


## Program `dna_rna_tools.py` 

### A tool to work with nucleotide sequences

Program `dna_rna_tools.py` automaticaly convert DNA or RNA sequences and calculate their *GC*-content. The conversion is performed according to [complementarity rule](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)). The calculation of *GC*-content is performed according to the [formula](https://en.wikipedia.org/wiki/GC-content#Determination).

### Usage

The programm is based on `run_dna_rna_tools` function that takes an arbitrary number of DNA or RNA sequences and a name of the function to be performed. The name of function must be inputed last. Use only one function at a time. Use only sequences that contain standard types of nucleotides: 

DNA - *A*, *G*, *C*, *T*

RNA - *A*, *G*, *C*, *U*

To start with the program run the following command:

`run_dna_rna_tools(sequences, func='function')`

Where:
- sequences - an arbitrary number of DNA or RNA sequences that must to be inputed in *string* type
- func - keyword argument, a type of functions to use that is inputed in *string* type
  
Before start, check the *Options* and *Examples*.


### Options
The program has five types of functions, for more information please see provided docstrings:

- `transcribe` — print transcribed sequence*
- `reverse` — print reversed sequence
- `complement` — print complementary sequence
- `reverse_complement` — print reversed complementary sequence
- `gc_count` — count *GC*-content in percentage
  
\* Reverse transcription is also taken into account (from RNA to DNA)
    
### Examples

```python
run_dna_rna_tools('ATG', func='transcribe') # 'AUG'
run_dna_rna_tools('ATG', func='reverse') # 'GTA'
run_dna_rna_tools('AtG', func='complement') # 'TaC'
run_dna_rna_tools('ATg', func='reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', func='reverse') # ['GTA', 'Ta']
run_dna_rna_tools('ATGgGCCtAA', func='gc_count') # '50.0%'
run_dna_rna_tools('ATTg', 'AuUgG', func='gc_count') # ['25.0%', '40.0%']
```

### Troubleshooting

|  Type of the problem                                             |  Probable cause
| ------------------------------------------------------------ |--------------------
| Output does not correspond the expected resultes             | The name of function is wrong. You see the results of another procedure
| run_dna_rna_tools() missing 1 required keyword-only argument: 'func'                          | The 'func' argument is not added
| ValueError: Invalid operation!                              | There is a mistake in the name of function
| ValueError: Sequence of number *n* is incorrect! | Sequence of this number *n* does not correspond to structure of DNA or RNA



## Program `protein_tools.py` 
### A tool to work with protein sequences

*Proteins* are under the constant focus of scientists. Currently, there are an enormous amount of tools to operate with nucleotide sequences, however, the same ones for proteins are extremely rare. 

`protein_tools.py` is an open-source program that facilitates working with protein sequences. 

### Usage
The programm is based on `run_protein_tools` function that takes an arbitrary number of **one-letter amino acid sequences**,  a name of procedure and a relevant argument. If you have three-letter amino acids sequences you could convert them by using `three_one_letter_code` procedure in advance. Three-lettter names of amino acids **must be separated with hyphen**.

To start with the program run the following command:

`run_protein_tools(sequences, procedure="procedure", ...)`

Where:
- sequences - an arbitrary number of amino acid sequences
- procedure - keyword argument, a type of procedure to use that is inputed in *string* type
- ... - an additional keyword arguments that are to be inputed in *string* type
- 
Before start, check the *Options* and *Examples*.
### Options

The program has five types of procedures, for more information please see provided docstrings:

 `three_one_letter_code`
 
- The main aim - to convert three-letter amino acid sequences to one-letter ones and vice-versa
- In case of three-to-one translation the names of amino acids **must be separated with hyphen**
- An additional argument: no

 `define_molecular_weight` 
 
- The main aim - to determine the exact molecular weight of protein sequences
- An additional argument: no

 `search_for_motifs` 

- The main aim - to search for the motif of interest in protein sequences
- An additional arguments: motif (*str*), overlapping (*bool*)

 `search_for_alt_frames` 

- The main aim - to look for alternative frames that start with methyonine or other non-canonical start amino acids
- Ignores the last three amino acids due to the insignicance of alternative frames of this length
- An additional argument: alt_start_aa (*str*)
- Use alt_start_aa **only for non-canonical start amino acids**
- Without alt_start_aa the procedure find alternative frames that start with methyonine

`convert_to_nucl_acids` 
 
- The main aim - to convert protein sequences to DNA, RNA or both nucleic acid sequences
- The program use the most frequent codons in human that could be found [here](https://www.genscript.com/tools/codon-frequency-table)
- An additional argument: nucl_acids (*str*)
- Use as nucl_acids only DNA, RNA or both (for more detailes, check *Examples*)


### Examples
```python
# three_one_letter_code
run_protein_tools('met-Asn-Tyr', 'Ile-Ala-Ala', procedure='three_one_letter_code')  # ['mNY', 'IAA']
run_protein_tools('mNY','IAA', procedure='three_one_letter_code')  # ['met-Asn-Tyr', 'Ile-Ala-Ala']


# define_molecular_weight
run_protein_tools('MNY','IAA', procedure='define_molecular_weight')  # {'MNY': 426.52, 'IAA': 273.35}


# check_for_motifs
run_protein_tools('mNY','IAA', procedure='search_for_motifs', motif='NY')
#Sequence: mNY
#Motif: NY
#Motif is present in protein sequence starting at positions: 1

#Sequence: IAA
#Motif: NY
#Motif is not present in protein sequence

{'mNY': [1], 'IAA': []}


# search_for_alt_frames
run_protein_tools('mNYQTMSPYYDMId', procedure='search_for_alt_frames')  # {'mNYQTMSPYYDMId': ['MSPYYDMId']}
run_protein_tools('mNYTQTSP', procedure='search_for_alt_frames', alt_start_aa='T')  # {'mNYTQTSP': ['TQTSP']}


# convert_to_nucl_acids
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'RNA')  # {'RNA': ['AUGAACUAU'], 'DNA': []}
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'DNA')  # {'RNA': [], 'DNA': ['ATGAACTAT']}
run_protein_tools('MNY', procedure='convert_to_nucl_acids', nucl_acids = 'both') # {'RNA': ['AUGAACUAU'], 'DNA': ['ATGAACTAT']}

```

### Troubleshooting

|  Type of the problem                                             |  Probable cause
| ------------------------------------------------------------ |--------------------
| Output does not correspond the expected resultes             | The name of procedure is wrong. You see the results of another procedure
| ValueError: No sequences provided                            | A list of sequences are not inputed
| ValueError: Wrong procedure                                  | The procedure does not exist in this program
| ValueError: Use three-letter aa only for "three_one_letter_code" procedure!  | Three-letter sequences were inputed. This type of sequences is accepted only by "three_one_letter_code" procedure
| ValueError: Invalid sequence given                           | The sequences do not correspond to standard amino acid code
| ValueError: Please provide desired motif                     | There are no an additional argument *motif* in `search_for_motifs`
| ValueError: Invalid alternative start AA                                 | There is more than one letter in an additional argument *alt_start_aa* in `search_for_alt_frames`
| ValueError: Please provide desired type of nucl_acids        | There are no an additional argument *nucl_acids* in `convert_to_nucl_acids`
| ValueError: Invalid nucl_acids argument                      | An additional argument in `convert_to_nucl_acids` is written incorrectly

## Contacts 

**Ekaterina Shitik** (shitik.ekaterina@gmail.com)

Doctor of medicine, molecular biologist with the main interests on gene engineering, AAV vectors and CRISPR/Cas9 technologies

**Vladimir Grigoriants** (vova.grig2002@gmail.com)

Bioinformatician, immunologist, MiLaborary inc. TCR-libraries QC developer 

**Vlada Tuliavko** (vladislavi2742@gmail.com)

MiLaboratory inc. manager&designer, immunologist

Приятного использования! ✨✨
