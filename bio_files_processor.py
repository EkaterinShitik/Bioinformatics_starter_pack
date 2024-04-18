import os
from dataclasses import dataclass


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta='result') -> None:
    """
    The function converts DNA/RNA/protein sequences from
    multiline format to one line
    The function takes fasta file as input_fasta argument
    The function output processed sequences in fasta file
    named as output_fasta argument
    It saves output file in a current directory

    Arguments:

    - input_fasta(str): the name of fasta file to be processed

    - output_fasta(str): the name for output file
    with obtained result (Optional)
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
    name_fasta = output_fasta + '.fasta'
    if os.path.exists(os.path.join(current_directory, name_fasta)):
        error = 'File with such name exists! Change output_fasta arg!'
        raise ValueError(error)
    with open(os.path.join(current_directory, name_fasta), mode='w') as file:
        for line in file_output:
            file.write(line + '\n')


def select_genes_from_gbk_to_fasta(input_gbk: str,
                                   *genes: str,
                                   n_before=1,
                                   n_after=1,
                                   output_fasta='result') -> None:
    """
    The main aim of function is to select genes that
    flanking the genes of main interest(*genes argument)
    The function takes a file with gbk extention as input_gbk argument
    The function output selected genes in fasta file
    named as output_fasta argument
    It saves output file in a current directory
    If selected genes are unknown they are named as unknownN
    N is a number of this gene in gbk file

    Arguments:

    - input_gbk(str): the name of gbk file to be processed

    - *genes(str): names of any number of genes to be studied

    - n_before(int): a number of upstream flanking genes (Optional)
    By default n_before=1
    To find more flanking genes provide full name of argument
    Example: n_before=2

    - n_after(int): a number of downstream flanking genes (Optional)
    By default n_after=1
    To find more flanking genes provide full name of argument
    Example: n_after=2

    - output_fasta(str): the name for output file with
    obtained result (Optional)
    By default output_fasta = 'result'
    Argument output_fasta must be inputed without fasta extension
    Use full name of argument
    Example: output_filename='processed_seqs'  # processed_seqs.fasta

    Return:
    - file: file with fasta extension containing names
    and protein sequences of genes flanking the genes of interest
    """
    annot_genes = []
    genes_information = []
    num = 1
    with open(input_gbk) as f:
        for line in f:
            line = line.strip()
            if line.startswith('CDS'):
                genes_information = ''.join(genes_information)
                genes_information = genes_information.split('"')
                if not any('gene' in _ for _ in genes_information):
                    name_annot_gene = 'unknown' + str(num)
                    num += 1
                    sequence_annot_gene = genes_information[-2]
                else:
                    name_annot_gene = genes_information[1]
                    sequence_annot_gene = genes_information[-2]
                annot_genes.append([name_annot_gene, sequence_annot_gene])
                genes_information = []
            else:
                genes_information.append(line)
    result = []
    for gene_in_annot in annot_genes:
        name_gene_in_annot = gene_in_annot[0]
        for gene in genes:
            if name_gene_in_annot == gene:
                posit_cur_gene = annot_genes.index(gene_in_annot)
                for previous_posit in range(n_before, 0, -1):
                    delta_position = posit_cur_gene - previous_posit
                    name_previous_gene = annot_genes[delta_position][0]
                    name_previous_for_fasta = '>gene ' + name_previous_gene
                    seq_previous_gene = annot_genes[delta_position][1]
                    result.append(name_previous_for_fasta)
                    result.append(seq_previous_gene)
                for next_posit in range(1, n_after+1):
                    sum_position = posit_cur_gene + next_posit
                    name_next_gene = annot_genes[sum_position][0]
                    name_next_for_fasta = '>gene ' + name_next_gene
                    seq_next_gene = annot_genes[sum_position][1]
                    result.append(name_next_for_fasta)
                    result.append(seq_next_gene)
    if len(result) == 0:
        return 'There are no provided genes in gbk file'
    current_directory = os.getcwd()
    name_fasta = output_fasta + '.fasta'
    if os.path.exists(os.path.join(current_directory, name_fasta)):
        error = 'File with such name exists! Change output_fasta arg!'
        raise ValueError(error)
    with open(os.path.join(current_directory, name_fasta), mode='w') as file:
        for line in result:
            file.write(line + '\n')


@dataclass
class FastaRecord:
    """
    The class to collect all information
    about a fasta record
    """
    id: str
    description: str
    seq: str

    def __repr__(self):
        return f'{self.id} {self.description} \n{self.seq}'


class OpenFasta:
    """
    The class is a context manager for fasta files
    It processes all reads and returns them separately
    function: read_record - returns each record step
    by step
    function: read_records - returns list of all records
    """
    def __init__(self, file_path: str, mode='r'):
        self.file_path = file_path
        self.mode = mode
        self.file_handler = None
        self.line = None

    def __enter__(self):
        self.file_handler = open(self.file_path, mode=self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self.file_handler:
            self.file_handler.close()

    def __iter__(self):
        return self

    def __next__(self) -> FastaRecord:
        if self.line is None:
            self.line = self.file_handler.readline()
        if self.line == '':
            raise StopIteration
        line_sep = self.line.strip().split(' ')
        seq_id = line_sep[0][1:]
        description = ' '.join(line_sep[1:])
        sequence = ''
        self.line = self.file_handler.readline()
        while not (self.line.startswith('>')):
            sequence += self.line.strip()
            self.line = self.file_handler.readline()
            if self.line == '':
                break
        fasta_record = FastaRecord(seq_id, description, sequence)
        return fasta_record

    def read_record(self) -> FastaRecord:
        """
        The function returns fasta records step by step
        :return: a fasta record
        """
        try:
            return self.__next__()
        except StopIteration:
            return FastaRecord('', '', '')

    def read_records(self) -> list:
        """
        The function collect all fasta records from
        a fasta file
        :return: list of fasta records
        """
        list_records = []
        for record in self.__iter__():
            list_records += [record]
        return list_records
