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
RULE_OF_TRANSCRIPTION = 'AUCG'.maketrans('UuTt', 'TtUu')
RULE_OF_TRANSLATION = 'MTW'.maketrans(TRANSLATION_CODE)


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
