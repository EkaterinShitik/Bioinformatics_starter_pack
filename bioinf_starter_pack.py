import src.dictionaries
import src.dna_rna_tools as dna_rna_tools
import src.protein_tools as protein_tools


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


PROCEDURES_TO_FUNCTIONS = {
    'search_for_motifs': protein_tools.search_for_motifs,
    'search_for_alt_frames': protein_tools.search_for_alt_frames,
    'convert_to_nucl_acids': protein_tools.convert_to_nucl_acids,
    'three_one_letter_code': protein_tools.three_one_letter_code,
    'define_molecular_weight': protein_tools.define_molecular_weight,
}


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
    return PROCEDURES_TO_FUNCTIONS[procedure](**procedure_arguments)
