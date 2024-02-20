import os
from typing import Optional

COMPLEMENT_RULE = 'ATUCG'.maketrans('AaTtUuCcGg', 'TtAaAaGgCc')
RULE_OF_TRANSCRIPTION = 'AUCG'.maketrans('UuTt', 'TtUu')
AMINO_ACIDS = {
    'A': 'Ala',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr',
    'a': 'ala',
    'c': 'cys',
    'd': 'asp',
    'e': 'glu',
    'f': 'phe',
    'g': 'gly',
    'h': 'his',
    'i': 'ile',
    'k': 'lys',
    'l': 'leu',
    'm': 'met',
    'n': 'asn',
    'p': 'pro',
    'q': 'gln',
    'r': 'arg',
    's': 'ser',
    't': 'thr',
    'v': 'val',
    'w': 'trp',
    'y': 'tyr',
}
TRANSLATION_CODE = {
    'F': 'UUU',
    'f': 'uuu',
    'L': 'CUG',
    'l': 'cug',
    'I': 'AUU',
    'i': 'auu',
    'M': 'AUG',
    'm': 'aug',
    'V': 'GUG',
    'v': 'gug',
    'P': 'CCG',
    'p': 'ccg',
    'T': 'ACC',
    't': 'acc',
    'A': 'GCG',
    'a': 'gcg',
    'Y': 'UAU',
    'y': 'uau',
    'H': 'CAU',
    'h': 'cau',
    'Q': 'CAG',
    'q': 'cag',
    'N': 'AAC',
    'n': 'aac',
    'K': 'AAA',
    'k': 'aaa',
    'D': 'GAU',
    'd': 'gau',
    'E': 'GAA',
    'e': 'gaa',
    'C': 'UGC',
    'c': 'ugc',
    'W': 'UGG',
    'w': 'ugg',
    'R': 'CGU',
    'r': 'cgu',
    'S': 'AGC',
    's': 'agc',
    'G': 'GGC',
    'g': 'ggc',
}
AMINO_ACID_WEIGHT = {
    'A': 89.09,
    'C': 121.16,
    'D': 133.10,
    'E': 147.13,
    'F': 165.19,
    'G': 75.07,
    'H': 155.16,
    'I': 131.17,
    'K': 146.19,
    'L': 131.17,
    'M': 149.21,
    'N': 132.12,
    'P': 115.13,
    'Q': 146.15,
    'R': 174.20,
    'S': 105.09,
    'T': 119.12,
    'V': 117.15,
    'W': 204.23,
    'Y': 181.19,
}
RULE_OF_TRANSLATION = 'MTW'.maketrans(TRANSLATION_CODE)


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
    gc_share = gc_counter / len(seq) * 100
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


FUNCTIONS_NA = {
    'transcribe': transcribe,
    'reverse': reverse,
    'complement': complement,
    'reverse_complement': reverse_complement,
    'gc_count': gc_count
}


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
    if func not in FUNCTIONS_NA:
        raise ValueError('Invalid operation!')
    for seq in seqs:
        if not (is_dna(seq)) and not (is_rna(seq)):
            number_seq = seqs.index(seq) + 1
            error = 'Sequence of number ' + str(number_seq) + ' is incorrect!'
            raise ValueError(error)
        result.append(FUNCTIONS_NA[func](seq))
    if len(result) == 1:
        return result[0]
    return result


def three_one_letter_code(sequences: str) -> list:
    """
    Reverse the protein sequences from one-letter to three-letter format and vice-versa

    Case 1: get three-letter sequence\n
    Use one-letter amino-acids sequences of any letter case

    Case 2: get one-letter sequence\n
    Use three-letter amino-acid separated by "-" sequences.
    Please note that sequences without "-" are parsed as one-letter code sequences\n
    Example of mistake: for sequence "Ala" function will return "Ala-leu-ala"

    Arguments:
    - sequences (tuple[str] or list[str]): protein sequences to convert\n
    Example: ['WAG', 'MkqRe', 'msrlk', 'Met-Ala-Gly', 'Met-arg-asn-Trp-Ala-Gly', 'arg-asn-trp']

    Return:
    - list: one-letter/three-letter protein sequences\n
    Example: ['Met-Ala-Gly', 'Met-arg-asn-Trp-Ala-Gly', 'arg-asn-trp', 'WAG', 'MkqRe', 'rlk']
    """
    inversed_sequences = []
    for sequence in sequences:
        if "-" not in sequence:
            inversed_sequence = []
            for letter in sequence:
                inversed_letter = AMINO_ACIDS[letter]
                inversed_sequence += [inversed_letter]
            inversed_sequence = '-'.join(inversed_sequence)
            inversed_sequences.append(inversed_sequence)
        else:
            inversed_sequence = ""
            aa_splitted = sequence.split("-")
            for aa in aa_splitted:
                dict_three_letter_aa = AMINO_ACIDS.values()
                dict_one_letter_aa = AMINO_ACIDS.keys()
                number_aa_in_dict = list(dict_three_letter_aa).index(aa)
                inversed_sequence += list(dict_one_letter_aa)[number_aa_in_dict]
            inversed_sequences.append(inversed_sequence)
    return inversed_sequences


def define_molecular_weight(sequences: str) -> dict:
    """
    Define molecular weight of the protein sequences

    Use one-letter amino-acids sequences of any letter case
    The molecular weight is:
    - a sum of masses of each atom constituting a molecule
    - expressed in units called daltons (Da)
    - rounded to hundredths

    Arguments:
    - sequences (tuple[str] or list[str]): protein sequences to convert

    Return:
    - dictionary: protein sequences as keys and molecular masses as values\n
    Example: {'WAG': 332.39, 'MkqRe': 690.88, 'msrlk': 633.86}
    """
    sequences_weights = {}
    for sequence in sequences:
        sequence_weight = 0
        for letter in sequence:
            sequence_weight += AMINO_ACID_WEIGHT[letter.upper()]
        sequence_weight -= (len(sequence) - 1) * 18  # deduct water from peptide bond
        sequences_weights[sequence] = round(sequence_weight, 2)
    return sequences_weights


def search_for_motifs(
        sequences: (tuple[str] or list[str]), motif: str, overlapping: bool
) -> dict:
    """
    Search for motifs - conserved amino acids residues in protein sequence

    Search for one motif at a time\n
    Search is letter case sensitive\n
    Use one-letter aminoacids code for desired sequences and motifs\n
    Positions of AA in sequences are counted from 0\n
    By default, overlapping matches are counted

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to check for given motif within\n
        Example: sequences = ['AMGAGW', 'GAWSGRAGA']
    - motif (str]: desired motif to check presense in every given sequence\n
        Example: motif='GA'
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
        Example: overlapping=False
    Return:
    - dictionary: sequences (str] as keys , starting positions for presented motif (list) as values\n
        Example: {'AMGAGW': [2], 'GAWSGRAGA': [0, 7]}
    """
    new_line = '\n'
    all_positions = {}
    for sequence in sequences:
        start = 0
        positions = []
        print(f"Sequence: {sequence}")
        print(f"Motif: {motif}")
        if motif in sequence:
            while True:
                start = sequence.find(motif, start)
                if start == -1:
                    break
                positions.append(start)
                if overlapping:
                    start += 1
                else:
                    start += len(motif)
            print_pos = ", ".join(str(x) for x in positions)
            print_pos = f"{print_pos}{new_line}"
            print(
                f"Motif is present in protein sequence starting at positions: {print_pos}"
            )
        else:
            print(f"Motif is not present in protein sequence{new_line}")
        all_positions[sequence] = positions
    return all_positions


def search_for_alt_frames(sequences: str, alt_start_aa: str) -> dict:
    """
    Search for alternative frames in a protein sequences

    Search is not letter case sensitive\n
    Without an alt_start_aa argument search for frames that start with methionine ('M')
    To search frames with alternative start codon add alt_start_aa argument\n
    In alt_start_aa argument use one-letter code

    The function ignores the last three amino acids in sequences

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to check
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start AA (Optional)\n
    Example: alt_start_aa='I'

    Return:
    - dictionary: the original sequences and a collection of alternative frames
    """
    alternative_frames = {}
    num_position = 0
    for sequence in sequences:
        alternative_frames[sequence] = []
        for amino_acid in sequence[1:-3]:
            alt_frame = ''
            num_position += 1
            if amino_acid == alt_start_aa or amino_acid == alt_start_aa.swapcase():
                alt_frame += sequence[num_position:]
                alternative_frames[sequence].append(alt_frame)
        num_position = 0
    return alternative_frames


def convert_to_nucl_acids(sequences: list, nucl_acids: str) -> dict:
    """
    Convert protein sequences to RNA or DNA sequences.

    Use the most frequent codons in human. The source - https://www.genscript.com/tools/codon-frequency-table\n
    All nucleic acids (DNA and RNA) are showed in 5'-3' direction
    The function outputs the coding DNA strand

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to convert
    - nucl_acids (str): the nucleic acid that is prefered\n
    Example: nucl_acids = 'RNA' - convert to RNA\n
             nucl_acids = 'DNA' - convert to DNA\n
             nucl_acids = 'both' - convert to RNA and DNA
    Return:
    - dictionary: 'DNA' or 'RNA'(str) as keys, collection of the obtained sequences (list) as values
    """
    nucl_acid_seqs = {'RNA': [], 'DNA': []}
    for sequence in sequences:
        rna_seq = sequence.translate(RULE_OF_TRANSLATION)
        dna_seq = rna_seq.translate(RULE_OF_TRANSCRIPTION)
        if nucl_acids == 'RNA':
            nucl_acid_seqs['RNA'].append(rna_seq)
        if nucl_acids == 'DNA':
            nucl_acid_seqs['DNA'].append(dna_seq)
        if nucl_acids == 'both':
            nucl_acid_seqs['RNA'].append(rna_seq)
            nucl_acid_seqs['DNA'].append(dna_seq)
    return nucl_acid_seqs


def check_and_parse_user_input(
        *sequences: str, **kwargs
) -> dict and str:
    """
    Check if user input can be correctly processed\n
    Parse sequences and arguments for desired procedure

    Arguments:
    - sequences (list[str] or tuple[str]): sequences to process
    - **kwargs - needed arguments for completion of desired procedure

    Return:
    - string: procedure name
    - dictionary: a collection of procedure arguments and their values
    """
    if len(sequences) == 0:
        raise ValueError('No sequences provided')
    procedure = kwargs['procedure']
    if procedure not in PROCEDURES_TO_FUNCTIONS.keys():
        raise ValueError('Wrong procedure')
    allowed_inputs = set(AMINO_ACIDS.keys()).union(
        set(AMINO_ACIDS.values())
    )
    allowed_inputs.add('-')
    for sequence in sequences:
        if '-' not in sequence:
            if not set(sequence).issubset(allowed_inputs):
                raise ValueError('Invalid sequence given')
        else:
            if procedure != 'three_one_letter_code':
                raise ValueError('Use three-letter aa only for "three_one_letter_code" procedure!')
            if not set(sequence.split('-')).issubset(allowed_inputs):
                raise ValueError('Invalid sequence given')
    procedure_arguments = {}
    if procedure == 'search_for_motifs':
        if 'motif' not in kwargs.keys():
            raise ValueError('Please provide desired motif')
        procedure_arguments['motif'] = kwargs['motif']
        if 'overlapping' not in kwargs.keys():
            procedure_arguments['overlapping'] = True
        else:
            procedure_arguments['overlapping'] = kwargs['overlapping']
    elif procedure == 'search_for_alt_frames':
        if 'alt_start_aa' not in kwargs.keys():
            procedure_arguments['alt_start_aa'] = 'M'
        else:
            if len(kwargs['alt_start_aa']) > 1:
                raise ValueError('Invalid alternative start AA')
            procedure_arguments['alt_start_aa'] = kwargs['alt_start_aa']
    elif procedure == 'convert_to_nucl_acids':
        if 'nucl_acids' not in kwargs.keys():
            raise ValueError('Please provide desired type of nucl_acids')
        if kwargs['nucl_acids'] not in {'DNA', 'RNA', 'both'}:
            raise ValueError('Invalid nucl_acids argument')
        procedure_arguments['nucl_acids'] = kwargs['nucl_acids']
    procedure_arguments['sequences'] = sequences
    return procedure_arguments, procedure


PROCEDURES_TO_FUNCTIONS = {
    'search_for_motifs': search_for_motifs,
    'search_for_alt_frames': search_for_alt_frames,
    'convert_to_nucl_acids': convert_to_nucl_acids,
    'three_one_letter_code': three_one_letter_code,
    'define_molecular_weight': define_molecular_weight,
}


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
    procedure_arguments, procedure = check_and_parse_user_input(*sequences, **kwargs)
    return PROCEDURES_TO_FUNCTIONS[procedure](**procedure_arguments)


def is_in_gc_bounds(seq: str, gc_bounds: tuple, gc_counter=0) -> bool:
    """
    Check if the sequence falls in the range of GC-content bounds
    The range of GC-content bounds is determined with gc_bounds argument

    Arguments:
    - seq(str): the sequence to check
    - gc_bounds(tuple): contain minimal and maximum GC-content bounds of the main interest
    - gc_counter(int): an additional argument to count GC-content that by default is 0 \n
    Do not change the last gc_counter argument!

    Example: is_in_gc_bounds('AgCC', gc_bounds=(10,90))

    Return:
    - bool: the result of the check
    """
    gc_min, gc_max = gc_bounds[0], gc_bounds[1]
    for nucl in seq:
        if nucl in ('G', 'C', 'g', 'c'):
            gc_counter += 1
    gc_share = gc_counter / len(seq) * 100
    return gc_min <= gc_share <= gc_max


def is_in_length_bounds(seq: str, length_bounds: tuple) -> bool:
    """
    Check if the sequence falls in the range of length bounds
    The range of length bounds is determined with length_bounds argument

    Arguments:
    - seq(str): the sequence to check
    - length_bounds(tuple): contain minimal and maximum length bounds of the main interest

    Example: is_in_length_bounds('AgCC', length_bounds=(0,12000))

    Return:
    - bool: the result of the check
    """
    length_min, length_max = length_bounds[0], length_bounds[1]
    return length_min <= len(seq) <= length_max


def is_above_quality_threshold(quality_scores: str, quality_threshold: float, sum_phred=0) -> bool:
    """
    Check if the mean of the sequence quality values in FASTQ exceeds the quality threshold of the interest
    The quality threshold is determined with quality_threshold argument
    To convert quality values, the function uses phred+33

    Arguments:
    - quality_scores(str): field 4 of FASTQ that encodes the quality values of the sequence
    - quality_threshold(int or float): the quality threshold of the interest
    - sum_phred(int): an additional argument to count phred+33 values that by default is 0\n
    Do not change the last sum_phred argument!

    Example: is_above_quality_threshold('BFFFFFFFB@B@A<@D', quality_threshold=20)

    Return:
    - bool: the result of the check
    """
    length_code = len(quality_scores)
    for symbol in quality_scores:
        sum_phred += (ord(symbol) - 33)
        mean_quality = sum_phred / length_code
    return mean_quality >= quality_threshold


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
    if type(length_bounds) is int:
        length_bounds = 0, length_bounds
    if type(gc_bounds) is int or type(gc_bounds) is float:
        gc_bounds = 0, gc_bounds
    result = {}
    for name in seqs.keys():
        seq = seqs[name][0]
        quality_scores = seqs[name][2]
        if (is_in_gc_bounds(seq, gc_bounds) and
            is_in_length_bounds(seq, length_bounds) and
            is_above_quality_threshold(quality_scores, quality_threshold)):
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
    if not (os.path.exists(output_path)):
        os.mkdir(output_path)
    if output_filename is None:
        input_filename = os.path.split(input_path)[-1]  # process output_filename
        output_filename = input_filename
    if not (output_filename.endswith('.fastq')):
        output_filename = output_filename + '.fastq'
    if os.path.exists(os.path.join(output_path, output_filename)):
        raise ValueError('File with such name exists! Change output_filename arg!')
    with open(os.path.join(output_path, output_filename), mode='w') as file:  # write in a new file
        for line in file_output:
            file.write(line + '\n') 
