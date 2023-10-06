import src.dictionaries
import src.dna_rna_tools as dna_rna_tools

FUNCTIONS = {
        'transcribe': dna_rna_tools.transcribe,
        'reverse': dna_rna_tools.reverse,
        'complement': dna_rna_tools.complement,
        'reverse_complement': dna_rna_tools.reverse_complement,
        'gc_count': dna_rna_tools.gc_count
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
