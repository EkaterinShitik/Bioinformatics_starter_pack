import src.dna_rna_tools as dna_rna_tools
import src.protein_tools as protein_tools
import src.fastq_tools as fastq_tools

def run_dna_rna_tools(*seqs: str, func: str) -> list[str]:
    """
    Main function to process nucleotide sequences by one of the developed tools.\n
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
        if not(dna_rna_tools.is_dna(seq)) and not(dna_rna_tools.is_rna(seq)):
            number_seq = seqs.index(seq) + 1
            error = 'Sequence of number ' + str(number_seq) + ' is incorrect!'
            raise ValueError(error)
        result.append(dna_rna_tools.FUNCTIONS_NA[func](seq))
    if len(result) == 1:
        return result[0]
    return result


def run_protein_tools(*sequences: str, **kwargs: str):
    """
    Main function to process protein sequence by one of the developed tools.\n
    Run one procedure at a time:
    - Search for conserved amino acids residues in protein sequence
    - Search for alternative frames in a protein sequences
    - Convert protein sequences to RNA or DNA sequences
    - Reverse the protein sequences from one-letter to three-letter format and vice-versa
    - Define molecular weight of the protein sequences

    All functions except *search_for_alt_frames* are letter case sensitive\n
    Provide protein sequence in one letter code.\n
    You can obtain one letter code from three letter code with *three_one_letter_code*\n
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
    - motif (str]: desired motif to check presense in every given sequence\n
            Example: motif='GA'
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
            Example: overlapping =False

    For 'search_for_alt_frames' procedure provide:
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start codon (Optional)\n
            Example: alt_start_aa='I'

    For 'convert_to_nucl_acids' procedure provide:
    - nucl_acids (str]: the nucleic acid to convert to\n
            Example: nucl_acids='RNA'\n
                           nucl_acids='DNA'\n
                           nucl_acids='both'

    Return:
    - dict: Dictionary with processed sequences. Depends on desired tool\n
            Please see Readme or desired docstring
    """
    procedure_arguments, procedure = protein_tools.check_and_parse_user_input(*sequences, **kwargs)
    return protein_tools.PROCEDURES_TO_FUNCTIONS[procedure](**procedure_arguments)


def select_fastq(seqs: dict, gc_bounds=(0,100), length_bounds=(0,2**32), quality_threshold=0) -> dict:
    """
    Main function to select fragmnets in FASTQ format according to three main requirements:\n
    -fall in the range of GC-content bounds
   The range of GC-content bounds is determined with gc_bounds argument
   
    -falls in the range of length bounds
    The range of length bounds is determined with length_bounds argument

    -exceeds the quality threshold of the interest
    The quality threshold is determined with quality_threshold argument

    Without full names use arguments in a certain order
    Example: select_fastq(seqs, (0,100), (0,200), 0)
           # select_fastq(seqs, gc_bounds=(0,100), length_bounds=(0,200), quality_threshold=0)
             
             
    In case of changing only one argument, provide its full name!
    Example: select_fastq(seqs, length_bounds=(50, 100))

             
    Arguments:
    
    - seqs (dict): the set of fragments to be analyzed
    Example: {'name' : ('sequence', 'quality values')}

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
    - dict: Dictionary with selected fragments
    
    For more information please see README
    
    """
    if type(length_bounds) == int or type(length_bounds) == float:
        length_bounds = 0, length_bounds
    if type(gc_bounds) == int or type(gc_bounds) == float: 
        gc_bounds = 0, gc_bounds
    result = {}
    for name in seqs.keys():
        seq = seqs[name][0]
        quality_scores = seqs[name][1]
        if (fastq_tools.is_in_gc_bounds(seq, gc_bounds) and
            fastq_tools.is_in_length_bounds(seq, length_bounds) and
            fastq_tools.is_above_quality_threshold(quality_scores, quality_threshold)):
            result[name] = seqs[name]
    if len(result) == 0:
        return 'There are no sequences suited to requirements'
    return result
