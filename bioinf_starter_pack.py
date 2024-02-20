import os
from typing import Optional
import src.dna_rna_tools as dna_rna_tools
import src.protein_tools as protein_tools
import src.fastq_tools as fastq_tools


def run_dna_rna_tools(*seqs: str, func: str) -> list[str]:
    """
    Main function to process nucleotide sequences by one of the developed tools
    
    Run one procedure at a time:
    - Reverse the sequences from 5'-3' direction to 3'-5' and vice versa
    - Output the complementary sequences
    - Output the complementary sequences in reverse direction
    - Translate RNA sequences to DNA and vice versa
    - Count the GC-content in the sequences

    Arguments:
    - *seqs(str): sequences to process
    - func(str): the operation of the interest:
        - 'transcribe'
        - 'reverse'
        - 'complement'
        - 'reverse_complement'
        - 'gc_count'

    All of the described func arguments correspond to following functions:
    - transcribe(seq)
    - reverse(seq)
    - complement(seq)
    - reverse_complement(seq)
    - gc_count(seq: str, gc_counter=0)

    All functions except *gc_count* are letter case sensitive

    Return:
    - list[str]: the result of processing several sequences
    - str: the result of processing one sequence

    For more details see README
    """
    result = []
    if func not in dna_rna_tools.FUNCTIONS_NA:
        raise ValueError('Invalid operation!')
    for seq in seqs:
        if not (dna_rna_tools.is_dna(seq)) and not (dna_rna_tools.is_rna(seq)):
            number_seq = seqs.index(seq) + 1
            error = 'Sequence of number ' + str(number_seq) + ' is incorrect!'
            raise ValueError(error)
        result.append(dna_rna_tools.FUNCTIONS_NA[func](seq))
    if len(result) == 1:
        return result[0]
    return result


def run_protein_tools(*sequences: str, **kwargs: str):
    """
    Main function to process protein sequence by one of the developed tools
    
    Run one procedure at a time:
    - Search for conserved amino acids residues in protein sequence
    - Search for alternative frames in a protein sequences
    - Convert protein sequences to RNA or DNA sequences
    - Reverse the protein sequences from one-letter to three-letter format and vice-versa
    - Define molecular weight of the protein sequences

    All functions except *search_for_alt_frames* are letter case sensitive
    
    Provide protein sequence in one letter code
    
    You can obtain one letter code from three letter code with *three_one_letter_code*
    
    If more information needed please see README or desired docstring

    Arguments:
    - sequences (str): sequences to process
    - procedure (str): desired procedure:
        - 'search_for_motifs'
        - 'search_for_alt_frames'
        - 'convert_to_nucl_acids'
        - 'three_one_letter_code'
        - 'define_molecular_weight'

    For 'search_for_motif' procedure provide:
    - motif (str]: desired motif to check presense in every given sequence
            Example: motif='GA'
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)
    
            Example: overlapping =False

    For 'search_for_alt_frames' procedure provide:
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start codon (Optional)
    
            Example: alt_start_aa='I'

    For 'convert_to_nucl_acids' procedure provide:
    - nucl_acids (str]: the nucleic acid to convert to
    
            Example: nucl_acids='RNA'
                           nucl_acids='DNA'
                           nucl_acids='both'

    Return:
    - dict: Dictionary with processed sequences. Depends on desired tool
    
            Please see Readme or desired docstring
    """
    procedure_arguments, procedure = protein_tools.check_and_parse_user_input(*sequences, **kwargs)
    return protein_tools.PROCEDURES_TO_FUNCTIONS[procedure](**procedure_arguments)


def select_fastq(input_path: str,
                 output_filename: Optional[str] = None,
                 gc_bounds=(0, 100),
                 length_bounds=(0, 2**32),
                 quality_threshold=0) -> None:
    """
    Main function to select fragmnets in FASTQ format according to three main requirements:\n
    -fall in the range of GC-content bounds
    The range of GC-content bounds is determined with gc_bounds argument
   
    -falls in the range of length bounds
    The range of length bounds is determined with length_bounds argument

    -exceeds the quality threshold of the interest
    The quality threshold is determined with quality_threshold argument

    Function takes the path to the file in input_path argument. Use files with fastq extension only.

    Function output the result of checking in a file that is named according to output_filename(Optional).

    The output file also has fatq extension.

    The output file is saved in the 'fastq_filtrator_results' directory.

    If 'fastq_filtrator_results' directory doesn't exist the program creates it in a current directory.

    Without full names use arguments in a certain order
    Example: select_fastq(input_path, output_filename, (0,100), (0,200), 0)
           # select_fastq(input_path, output_filename, gc_bounds=(0,100), length_bounds=(0,200), quality_threshold=0)
             
             
    In case of changing only one argument, provide its full name!
    Example: select_fastq(input_path, output_filename, length_bounds=(50, 100))

             
    Arguments:
    
    - input_path(str): the path to the file

    - output_filename(str): the name for output file with obtained result
    By default output_filename=None
    Without output_filename argument the output file is named as input file
    Name without fastq extention is acceptible.
    Example: output_filename='result'  # 'result.fastq'
             output_filename='result.fastq'

    - gc_bounds(tuple or int or float): contain minimal and maximum GC-content bounds
    By default gc_bounds=(0,100)
    If input contains one number the function accepts it as a maximum bound
    Examples: gc_bounds=(20,40)
              gc_bounds=40  # (0,40)

    - length_bounds(tuple or int or float): contain minimal and maximum length bounds
    By default length_bounds=(0,4294967296)
    If input contains one number the function accepts it as a maximum bound
    Examples: length_bounds=(10,90)
              length_bounds=90  # (0,90)

    - quality_threshold(int or float): the quality threshold of the interest
    By default quality_threshold=0
    Examples: quality_threshold=10

    There are three functions that are used in the main function:
    
    - is_in_gc_bounds(seq, gc_bounds): 
    Check if the sequence falls in the range of GC-content bounds
    
    - is_in_length_bounds(seq, length_bounds): 
    Check if the sequence falls in the range of length bounds
    
    - is_above_quality_threshold(quality_scores, quality_threshold): 
    Check if the mean of quality values exceeds the quality threshold
    
    Return:
    - file: file with fastq extension containing selected fragments.
    
    For more information please see README
    
    """
    seqs = {}
    with open(input_path) as f:
        name = f.readline().strip()
        seqs[name] = []
        for line in f:
            line = line.strip()
            if line.startswith('@'):
                if len(seqs[name]) == 2:
                    seqs[name].append(line)
                else:
                    name = line
                    seqs[name] = []
            else:
                seqs[name].append(line)
    if type(length_bounds) == int:
        length_bounds = 0, length_bounds
    if type(gc_bounds) == int or type(gc_bounds) == float: 
        gc_bounds = 0, gc_bounds
    result = {}
    for name in seqs.keys():
        seq = seqs[name][0]
        quality_scores = seqs[name][2]
        if (fastq_tools.is_in_gc_bounds(seq, gc_bounds) and
            fastq_tools.is_in_length_bounds(seq, length_bounds) and
            fastq_tools.is_above_quality_threshold(quality_scores, quality_threshold)):
            result[name] = seqs[name]
    if len(result) == 0:
        return 'There are no sequences suited to requirements'
    file_output = []
    for name, [seq, comment, quality] in result.items():
        (file_output.append(name) or 
         file_output.append(seq) or
         file_output.append(comment) or
         file_output.append(quality))  # make list
    current_directory = os.getcwd()
    output_path = os.path.join(current_directory, 'fastq_filtrator_results')  # determine the path to files
    if not(os.path.exists(output_path)):
        os.mkdir(output_path)
    if output_filename is None:
        input_filename = os.path.split(input_path)[-1]  # process output_filename
        output_filename = input_filename
    if not(output_filename.endswith('.fastq')):
        output_filename = output_filename + '.fastq'
    if os.path.exists(os.path.join(output_path, output_filename)):
        raise ValueError('File with such name exists! Change output_filename arg!')
    with open(os.path.join(output_path, output_filename), mode='w') as file:  # write in a new file
        for line in file_output:
            file.write(line + '\n') 
