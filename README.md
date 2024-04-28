# Bioinformatics_starter_pack
Bioinformatic_starter_pack contains programs that were written while doing home tasks in Python course in [Bioinformatics Institute](https://bioinf.me/en).
<br/><br/> 
![Image alt](https://github.com/EkaterinShitik/the-second-repository/raw/main/image_python_rep.png)
<br/><br/> 
The examples of running some programs that are highlighted with * are presented in `Showcases.ipynb`

## `bioinf_starter_pack.py`
The program `bioinf_starter_pack.py` provides basic tools to work with different types of bioinformatic data:
- `RNASequence`/`DNASequence`/`AminoAcidSequence` classes*
  
  To manipulate with three main bioinformatic data types - *DNA*, *RNA* and *proteins*. The implementation of these classes embraced four basic principles of OOP - encapsulation, inheritance, polymorphism, and abstraction.
- `filter_fastq` function

  To filter fastq files under the GC-content, the length of sequences and the quality scores
- `run_genscan` function*
  
  To make requests on [Genscan website](http://hollywood.mit.edu/GENSCAN.html) to predict possible CDS, exons and introns in the sequence of interest. The output is also implemented as a class.
  
- `telegram_logger` decorator
  
  To send the user information on the function of the main interest. The information includes the result of run, time of running, stdout and stderr output in log file. This function is implemented based on [Telegram API](https://core.telegram.org/) without specific libraries.
      
## `bio_files_processor.py`
## `custom_random_forest.py`
## `test_bioinf_starter_pack.py`
Enjoy your use! ✨✨
