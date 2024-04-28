# Bioinformatics_starter_pack
Bioinformatic_starter_pack contains programs that were written while doing home tasks in Python course in [Bioinformatics Institute](https://bioinf.me/en).
<br/><br/> 
![Image alt](https://github.com/EkaterinShitik/the-second-repository/raw/main/image_python_rep.png)
<br/><br/> 
The examples of running some programs that are highlighted with * are presented in `Showcases.ipynb`

## `bioinf_starter_pack.py`
The program `bioinf_starter_pack.py` provides basic tools to work with different types of bioinformatic data.
- `RNASequence`/`DNASequence`/`AminoAcidSequence` classes *
  
  To manipulate with three main bioinformatic data types - *DNA*, *RNA* and *proteins*. The implementation of these classes embraced four basic principles of OOP - encapsulation, inheritance, polymorphism, and abstraction.
- `filter_fastq` function

  To filter fastq files under the GC-content, the length of sequences and the quality scores
- `run_genscan` function *
  
  To make requests on [Genscan website](http://hollywood.mit.edu/GENSCAN.html) to predict possible CDS, exons and introns in the sequence of interest. The output is also implemented as a class.
  
- `telegram_logger` decorator
  
  To send the user information on the function of the main interest. The information includes the result of run, time of running, stdout and stderr output in log file. This function is implemented based on [Telegram API](https://core.telegram.org/) without specific libraries.
      
## `bio_files_processor.py`
The program `bio_files_processor.py` facilitates parsing FASTA and GenBank files.
- `convert_multiline_fasta_to_oneline` function
  
  To convert any number of *DNA*/*RNA*/*protein* sequences in FASTA file from multiline format to one line.
- `select_genes_from_gbk_to_fasta` function
  
  To select nearest neighbors of the gene of main interest from GenBank file and output them in FASTA format.
- `OpenFasta` context manager *

  To open FASTA file and return separate FASTA records including id, description and sequence. The implementation of the context manager is similar to `open` built-in function.
## `custom_random_forest.py`
The program `custom_random_forest.py` contains custom implementation of Random forest classifier as `RandomForestClassifierCustom` class *. 

In terms of this course this class was refined using parallel programming. Using `n_jobs` argument enables speeding up the running 
## `test_bioinf_starter_pack.py`
The program `test_bioinf_starter_pack.py` includes classes to test the programs mentioned above.

Enjoy your use! ✨✨
