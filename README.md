# Bioinformatics_starter_pack

![Image alt](https://github.com/EkaterinShitik/the-second-repository/raw/main/Презентация-без-названия.jpg)


## Программа `dna_rna_tools.py` 

### Общая информация

Программа `dna_rna_tools.py` выполяет автоматизированный перевод последовательностей ДНК или РНК и подсчёт их *GC* состава. Перевод осуществляется в соответствии с принципом [комплементарности](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)). Подсчёт *GC* пар выполняется по формуле, которую можно посмотреть [здесь](https://en.wikipedia.org/wiki/GC-content#Determination).

#### Подробное описание

Основу программы `dna_rna_tools.py` составляет функция `run_dna_rna_tools`, которая принимает на вход произвольное количество последовательностей нуклеотидов ДНК или РНК и название процедуры, которую необходимо выполнить. 

**Список процедур:**

- `transcribe` — напечатать транскрибированную последовательность*
- `reverse` — напечатать перевёрнутую последовательность
- `complement` — напечатать комплементарную последовательность
- `reverse_complement` — напечатать обратную комплементарную последовательность
- `gc_count` — посчитать содержание нуклеотидов *G* и *C* в процентах
  
\* Обратная транскрипция в рамках данной процедуры также учитывается (РНК в ДНК)
    
**Пример использования**

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
run_dna_rna_tools('ATGgGCCtAA', 'gc_count') # 1 50.0%
run_dna_rna_tools('ATTg', 'AuUgG', 'gc_count') # 1 25.0% 2 40.0%
```

**Требования:**

- Все вводимые аргументы должны быть представлены в виде строковых переменных (*str*)
- Последовательности должны содержать только стандартные нуклеотиды: ДНК - *A*, *G*, *C*, *T*; РНК - *A*, *G*, *C*, *U*
- Название процедуры всегда указывается последним


## Protein_tools.py
### A tool to work with protein sequences

*Proteins* are under the constant focus of scientists. Currently, there are an enormous amount of tools to operate with nucleotide sequences, however, the same ones for proteins are extremely rare. 


`protein_tools.py` is an open-source program that facilitates working with protein sequences. 

### Usage
The programm is based on `run_protein_tools` function that takes the list of **one-letter amino acid sequences**,  a name of procedure and a relevant argument. If you have three-letter amino acids sequences you could convert them by using `three_one_letter_code` procedure in advance. Please convert your three-letter coded sequences with `three_one_letter_code` procedure before using any other procedures on them.

To start with the program run the following command:

`run_protein_tools(sequences, procedure="procedure", ...)`

Where:
- sequences - positional argument, a list of protein sequences
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
- 

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
```

`convert_to_nucl_acids` 
 


- The main aim - to convert protein sequences to DNA, RNA or both nucleic acid sequences
- The program use the most frequent codons in human that could be found [here](https://www.genscript.com/tools/codon-frequency-table)
- An additional argument: nucl_acids (*str*)
- Use as nucl_acids only DNA, RNA or both (for more detailes, check *Examples*)
```
"""
Convert protein sequences to RNA or DNA sequences.

Use the most frequent codons in human. The source - https://www.genscript.com/tools/codon-frequency-table\n
All nucleic acids (DNA and RNA) are showed in 5"-3" direction

Arguments:
- sequences (tuple[str] or list[str]): sequences to convert
- nucl_acids (str]: the nucleic acid that is prefered\n
Example: nucl_acids = "RNA" - convert to RNA\n
               nucl_acids = "DNA" - convert to DNA\n
               nucl_acids = "both" - convert to RNA and DNA
Return:
- dictionary: nucleic acids (str) as keys, collection of sequences (list) as values
"""
```

## Examples
```python
# three_one_letter_code
run_protein_tools(['met-Asn-Tyr', 'Ile-Ala-Ala'], procedure='three_one_letter_code')  # ['mNY', 'IAA']
run_protein_tools(['mNY','IAA'], procedure='three_one_letter_code')  # ['met-Asn-Tyr', 'Ile-Ala-Ala']


# define_molecular_weight
run_protein_tools(['MNY','IAA'], procedure='define_molecular_weight')  # {'MNY': 426.52, 'IAA': 273.35}


# check_for_motifs
run_protein_tools(['mNY','IAA'], procedure='search_for_motifs', motif='NY')
#Sequence: mNY
#Motif: NY
#Motif is present in protein sequence starting at positions: 1

#Sequence: IAA
#Motif: NY
#Motif is not present in protein sequence

{'mNY': [1], 'IAA': []}


# search_for_alt_frames
run_protein_tools(['mNYQTMSPYYDMId'], procedure='search_for_alt_frames')  # {'mNYQTMSPYYDMId': ['MSPYYDMId']}
run_protein_tools(['mNYTQTSP'], procedure='search_for_alt_frames', alt_start_aa='T')  # {'mNYTQTSP': ['TQTSP']}


# convert_to_nucl_acids
run_protein_tools(['MNY'], procedure='convert_to_nucl_acids', nucl_acids = 'RNA')  # {'RNA': ['AUGAACUAU']}
run_protein_tools(['MNY'], procedure='convert_to_nucl_acids', nucl_acids = 'DNA')  # {'DNA': ['TACTTGATA']}
run_protein_tools(['MNY'], procedure='convert_to_nucl_acids', nucl_acids = 'both') # {'RNA': ['AUGAACUAU'], 'DNA': ['TACTTGATA']}

```

### Troubleshooting

|  Type of the problem                                             |  Probable cause
| ------------------------------------------------------------ |--------------------
| Output does not correspond the expected resultes             | The name of procedure is wrong. You see the results of another procedure
| ValueError: No sequences provided                            | A list of sequences are not inputed
| ValueError: Wrong procedure                                  | The procedure does not exist in this program
| TypeError: takes from 0 to 1 positional arguments but n were given  | Sequences are not collected into the list type
| ValueError: Invalid sequence given                           | The sequences do not correspond to standard amino acid code
| ValueError: Please provide desired motif                     | There are no an additional argument *motif* in `search_for_motifs`
| ValueError: Invalid start AA                                 | There is more than one letter in an additional argument *alt_start_aa* in `search_for_alt_frames`
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
