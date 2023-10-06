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
    gc_share = gc_counter/len(seq)*100 
    return gc_min <= gc_share <=gc_max


def is_in_length_bounds(seq:str, length_bounds: tuple) -> bool:
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
        mean_quality = sum_phred/length_code
    return mean_quality >= quality_threshold


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
        if (is_in_gc_bounds(seq, gc_bounds) and
            is_in_length_bounds(seq, length_bounds) and
            is_above_quality_threshold(quality_scores, quality_threshold)):
            result[name] = seqs[name]
    if len(result) == 0:
        return 'There are no sequences suited to requirements'
    return result
