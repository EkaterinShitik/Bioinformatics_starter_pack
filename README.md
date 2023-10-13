# Bioinformatics_starter_pack
Bioinformatic_starter_pack contains two type of tools that could be useful for bioinformaticians.

**Program `bioinf_starter_pack.py`**

The program `bioinf_starter_pack.py` provides basic tools to work with three main types of bioinformatic data:
- `dna_rna_tools.py` to process *nucleotide* *sequences*
- `protein_tools.py` to process *amino* *acid* *sequences*
- `fastq_tools.py` to filtrate the obtained sequences in *FASTQ format*


![Image alt](https://github.com/EkaterinShitik/the-second-repository/raw/main/Презентация-без-названия.jpg)

**Program `bio_files_processor.py`**

The program `bio_files_processor.py.` provides basic functions to work with fasta and gbk files:
- `convert_multiline_fasta_to_oneline` function converts DNA/RNA/protein sequences from multiline format to one line
- `select_genes_from_gbk_to_fasta` function to select genes that flanking the genes of main interest

## Program `bioinf_starter_pack.py`
### Program `dna_rna_tools.py` 

#### A tool to work with nucleotide sequences

Program `dna_rna_tools.py` automaticaly convert DNA or RNA sequences and calculate their *GC*-content. The conversion is performed according to [complementarity rule](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)). The calculation of *GC*-content is performed according to the [formula](https://en.wikipedia.org/wiki/GC-content#Determination).

#### Usage

The programm is based on `run_dna_rna_tools` function that takes an arbitrary number of DNA or RNA sequences and a name of the function to be performed. The name of function must be inputed last. Use only one function at a time. Use only sequences that contain standard types of nucleotides: 

DNA - *A*, *G*, *C*, *T*

RNA - *A*, *G*, *C*, *U*

To start with the program run the following command:

`run_dna_rna_tools(sequences, func='function')`

Where:
- sequences - an arbitrary number of DNA or RNA sequences that must to be inputed in *string* type
- func - keyword argument, a type of functions to use that is inputed in *string* type
  
Before start, check the *Options* and *Examples*.


#### Options
The program has five types of functions, for more information please see provided docstrings:

- `transcribe` — print transcribed sequence*
- `reverse` — print reversed sequence
- `complement` — print complementary sequence
- `reverse_complement` — print reversed complementary sequence
- `gc_count` — count *GC*-content in percentage
  
\* Reverse transcription is also taken into account (from RNA to DNA)
    
#### Examples

```python
run_dna_rna_tools('ATG', func='transcribe') # 'AUG'
run_dna_rna_tools('ATG', func='reverse') # 'GTA'
run_dna_rna_tools('AtG', func='complement') # 'TaC'
run_dna_rna_tools('ATg', func='reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', func='reverse') # ['GTA', 'Ta']
run_dna_rna_tools('ATGgGCCtAA', func='gc_count') # '50.0%'
run_dna_rna_tools('ATTg', 'AuUgG', func='gc_count') # ['25.0%', '40.0%']
```

#### Troubleshooting

|  Type of the problem                                             |  Probable cause
| ------------------------------------------------------------ |--------------------
| Output does not correspond the expected resultes             | The name of function is wrong. You see the results of another procedure
| run_dna_rna_tools() missing 1 required keyword-only argument: 'func'                          | The 'func' argument is not added
| ValueError: Invalid operation!                              | There is a mistake in the name of function
| ValueError: Sequence of number *n* is incorrect! | Sequence of this number *n* does not correspond to structure of DNA or RNA



### Program `protein_tools.py` 
#### A tool to work with protein sequences

*Proteins* are under the constant focus of scientists. Currently, there are an enormous amount of tools to operate with nucleotide sequences, however, the same ones for proteins are extremely rare. 

`protein_tools.py` is an open-source program that facilitates working with protein sequences. 

#### Usage
The programm is based on `run_protein_tools` function that takes an arbitrary number of **one-letter amino acid sequences**,  a name of procedure and a relevant argument. If you have three-letter amino acids sequences you could convert them by using `three_one_letter_code` procedure in advance. Three-lettter names of amino acids **must be separated with hyphen**.

To start with the program run the following command:

`run_protein_tools(sequences, procedure="procedure", ...)`

Where:
- sequences - an arbitrary number of amino acid sequences
- procedure - keyword argument, a type of procedure to use that is inputed in *string* type
- ... - an additional keyword arguments that are to be inputed in *string* type
- 
Before start, check the *Options* and *Examples*.
#### Options

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


#### Examples
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

#### Troubleshooting

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

### Program `fastq_tools.py` 

#### A tool to work with nucleotide fragments in FASTQ format

Program `fastq_tools.py` selects nucleotide fragments in FASTQ format according to three requirements that could be determined by User:
- *GC*-content. The calculation of *GC*-content is performed according to the [formula](https://en.wikipedia.org/wiki/GC-content#Determination).
- *The* *length* of fragment sequences
- *The* *quality* *score* of sequencing. The program uses phred+33 score and converts the quality values according to this [rule](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)

#### Usage

The program is based on `select_fastq` function that takes the dictionary of nucleotide fragments in FASTQ format, the range of maximal and minimal bounds of GC-content, the range of maximal and minimal bounds of the length of fragments and the quality threshold. All of the keyword arguments have their default meaning(see below). 

To start with the program run the following command:

`select_fastq(seqs, gc_bounds=..., length_bounds=..., quality_threshold=...)`

Where:
- seqs - the dictionary of nucleotide fragments in FASTQ format
  
  Example: {'name' : ('sequence', 'quality values')}
  
- gc_bounds - keyword argument that determines maximal and minimal bounds of GC-content
  
  By default `gc_bounds=(0,100)`
  
  This argument could be inputed in *tuple*, *int* or *float* types. In case of using one number the function accepts it as 
  a maximum bound
  
  Example: `gc_bounds=40  # (0,40)`
  
- length_bounds - keyword argument that determines maximal and minimal bounds of the length of fragments
  
  By default `length_bounds=(0,4294967296)`
  
  This argument could be inputed in *tuple*, *int* or *float* types. In case of using one number the function accepts it as 
  a maximum bound
  
  Example: `length_bounds=90.8  # (0,90.8)`

- quality_threshold(int or float): the quality threshold of the main interest. It could be inputed in *int* and *float* types
  
  By default `quality_threshold=0`
  
  Example: `quality_threshold=10`

Without full names use arguments in a certain order

Example: 

`select_fastq(seqs, (0,100), (0,200), 0)`  # *the* *same* *as* 

`select_fastq(seqs, gc_bounds=(0,100), length_bounds=(0,200), quality_threshold=0)`
             
             
In case of using only one argument, provide its full name!

Example: `select_fastq(seqs, length_bounds=(50, 100))`
  
Before start, check the *Options* and *Examples*.


#### Options

There are three functions that are used in the program:
    
- is_in_gc_bounds(seq, gc_bounds) - takes a sequence in *string* type and check if the sequence falls in the range of GC-content bounds
- is_in_length_bounds(seq, length_bounds) - takes a sequence in *string* type and check if the sequence falls in the range of length bounds
- is_above_quality_threshold(quality_scores, quality_threshold) - takes a quality values in *string* type and check if the mean of quality values exceeds the quality threshold
    
#### Examples

```python
example = {'@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
 '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
  'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')}

select_fastq(example, gc_bounds=(40,60))
# Output:
# {'@SRX079804:1:SRR292678:1:1101:24563:24563':
# ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
#  'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')}

select_fastq(example, quality_threshold=90)
# Output:
# 'There are no sequences suited to requirements'
```

#### Troubleshooting

|  Type of the problem                                             |  Probable cause
| ------------------------------------------------------------ |--------------------
| Output does not correspond the expected resultes             | The arguments are inputed without full names in the incorrect order
| AttributeError: 'some' object has no attribute 'keys'                         | The 'seq' argument has incorrect 'some' type. It must be in dictionary type
| TypeError: select_fastq() got an unexpected keyword argument 'n'   | The name of argument 'n' is written incorrectly
| TypeError: '>=' not supported between instances of 'float' and 'str'| The arguments are inputed in incorrect type


## Program `bio_files_processor.py`
### Function `convert_multiline_fasta_to_oneline`
### Function `select_genes_from_gbk_to_fasta`

## Contacts 

**Ekaterina Shitik** (shitik.ekaterina@gmail.com)

Doctor of medicine, molecular biologist with the main interests on gene engineering, AAV vectors and CRISPR/Cas9 technologies

**Vladimir Grigoriants** (vova.grig2002@gmail.com)

Bioinformatician, immunologist, MiLaborary inc. TCR-libraries QC developer 

**Vlada Tuliavko** (vladislavi2742@gmail.com)

MiLaboratory inc. manager&designer, immunologist

Enjoy your use! ✨✨
