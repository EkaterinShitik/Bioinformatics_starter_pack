import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta='result') -> str:
    """
    The function converts DNA/RNA/protein sequences from multiline format to one line

    The function takes fasta file as input_fasta argument

    The function output processed sequences in fasta file named as output_fasta argument

    It saves output file in a current directory

    Arguments:

    - input_fasta(str): the name of fasta file to be processed
    
    - output_fasta(str): the name for output file with obtained result (Optional)
    By default output_fasta = 'result'
    Argument output_fasta must be inputed without fasta extension
    Example: output_filename='processed_seqs'  # processed_seqs.fasta
    
    Return:
    - file: file with fasta extension containing processed sequences
    """
    seqs = {}
    with open(input_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line
                seqs[name] = ''
            else:
                seqs[name] += line
    file_output = []
    for name, seq in seqs.items():
        (file_output.append(name) or 
        file_output.append(seq))
    current_directory = os.getcwd()
    with open(os.path.join(current_directory, output_fasta + '.fasta'), mode='w') as file:
        for line in file_output:
            file.write(line + '\n')
