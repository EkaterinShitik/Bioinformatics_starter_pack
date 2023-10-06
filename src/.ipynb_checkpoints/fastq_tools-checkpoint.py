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
    