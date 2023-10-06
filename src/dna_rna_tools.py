COMPLEMENT_RULE = 'ATUCG'.maketrans('AaTtUuCcGg', 'TtAaAaGgCc')
RULE_OF_TRANSCRIPTION = 'AUCG'.maketrans('UuTt', 'TtUu')


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
    complement_seq = seq.translate(COMPLEMENT_RULE)
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
    transcribe_seq = seq.translate(RULE_OF_TRANSCRIPTION)
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

is_rna('AGCT')