import pandas as pd
import numpy as np
from enum import Enum

# define cigar operations
class cigar_op( Enum ):
    # Op BAM Description Consumes query Consumes reference
    M = 0 # alignment match (can be a sequence match or mismatch) yes yes
    I = 1 # insertion to the reference yes no
    D = 2 # deletion from the reference no yes
    N = 3 # skipped region from the reference no yes
    S = 4 # soft clipping (clipped sequences present in SEQ) yes no
    H = 5 # hard clipping (clipped sequences NOT present in SEQ) no no
    P = 6 # padding (silent deletion from padded reference) no no
    EQUAL = 7 # sequence match yes yes
    X = 8 # sequence mismatch yes yes

def Generate_Kmer(seq, window_size):
    """
    # 2021-02-20 15:14:13
    generate a list of Kmers from the sequence with the given window_size"""
    return list(
        seq[i : i + window_size] for i in range(0, len(seq) - window_size + 1, 1)
    )


# In[ ]:


def Reverse_Complement(seq):
    """# 2021-02-04 11:47:19
    Return reverse complement of 'seq'"""
    dict_dna_complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-"}
    return "".join(list(dict_dna_complement[base] for base in seq))[::-1]