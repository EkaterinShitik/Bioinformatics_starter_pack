import dictionaries

def reverse(seq: str) -> str:
    """
    Reverse the sequence from 5'-3' direction to 3'-5' and vice versa
    
    Arguments:
    - seq(str): the sequence to change
    
    Return:
    - str: the changed sequence
    """
    reverse_seq = seq[::-1]
    return reverse_seq


def complement(seq: str) -> str:
    """
    Output the complementary sequence
    The complementarity rule could be found here:
    https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)

    Arguments:
    - seq(str): the sequence to change
    
    Return:
    - str: the changed sequence
    """
    complement_seq = seq.translate(dictionaries.COMPLEMENT_RULE)
    return complement_seq


def reverse_complement(seq: str) -> str:
    """
    Output the complementary sequence in reverse direction

    This function uses two other functions:
    - complement(seq): output the complementary sequence
    - reverse(seq): reverse the sequence from 5'-3' direction to 3'-5' and vice versa
    
    Arguments:
    - seq(str): the sequence to change
    
    Return:
    - str: the changed sequence
    """
    reverse_complement_seq = reverse(complement(seq))
    return reverse_complement_seq


def transcribe(seq: str) -> str:
    """
    Translate RNA sequence to DNA and vice versa
    
    Arguments:
    - seq(str): the sequence to change
    
    Return:
    - str: the changed sequence
    """
    transcribe_seq = seq.translate(dictionaries.RULE_OF_TRANSCRIPTION)
    return transcribe_seq


def gc_count(seq: str, gc_counter=0) -> str:
    """
    Count the GC-content in the sequence
    
    Arguments:
    - seq(str): the sequence to count GC-content
    - gc_counter(int): an additional argument to count GC-content that by default is 0
    
    Do not change the last gc_counter argument!
    
    Return:
    - str: the percentage of GC-content in the sequence
    """
    for nucl in seq:
        if nucl in ('G', 'C', 'g', 'c'):
            gc_counter += 1
    gc_share = gc_counter/len(seq)*100
    return str(gc_share) + '%'


def is_dna(seq: str) -> bool:
    """
    Check if the sequence is DNA
    
    Arguments:
    - seq(str): the sequence to check

    Return:
    - bool: the result of the check
    """
    return set(seq) <= set('AaGgCcTt')


def is_rna(seq: str) -> bool:
    """
    Check if the sequence is RNA
    
    Arguments:
    - seq(str): the sequence to check

    Return:
    - bool: the result of the check
    """
    return set(seq) <= set('AaGgCcUu')


FUNCTIONS = {
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement,
        'gc_count': gc_count
        }


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
    if func not in FUNCTIONS:
        raise ValueError('Invalid operation!')
    for seq in seqs:
        if not(is_dna(seq)) and not(is_rna(seq)):
            number_seq = seqs.index(seq) + 1
            error = 'Sequence of number ' + str(number_seq) + ' is incorrect!'
            raise ValueError(error)
        result.append(FUNCTIONS[func](seq))
    if len(result) == 1:
        return result[0]
    return result
