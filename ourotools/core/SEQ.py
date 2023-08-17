import pandas as pd
import numpy as np
from enum import Enum
from uuid import uuid4

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

def UUID():
    """return a 128bit universal unique identifier"""
    return uuid4().hex
    
def Generate_Kmer(seq, window_size):
    """
    # 2021-02-20 15:14:13
    generate a list of Kmers from the sequence with the given window_size"""
    return list(
        seq[i : i + window_size] for i in range(0, len(seq) - window_size + 1, 1)
    )

def COUNTER(l_values, dict_counter=None, ignore_float=True):  # 2020-07-29 23:49:51
    """Count values in l_values and return a dictionary containing count values. if 'dict_counter' is given, countinue counting by using the 'dict_counter'. if 'ignore_float' is True, ignore float values, including np.nan"""
    if dict_counter is None:
        dict_counter = dict()
    if ignore_float:  # if 'ignore_float' is True, ignore float values, including np.nan
        for value in l_values:
            if isinstance(value, float):
                continue  # ignore float values
            if value in dict_counter:
                dict_counter[value] += 1
            else:
                dict_counter[value] = 1
    else:  # faster counting by not checking type of value
        for value in l_values:
            if value in dict_counter:
                dict_counter[value] += 1
            else:
                dict_counter[value] = 1
    return dict_counter

def Reverse_Complement(seq):
    """# 2021-02-04 11:47:19
    Return reverse complement of 'seq'"""
    dict_dna_complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-"}
    return "".join(list(dict_dna_complement[base] for base in seq))[::-1]


def Cluster_with_Kmer(
    dict_seq_count,
    int_min_n_overlap_kmer,
    len_kmer,
    float_min_proportion_read_to_select_kmer_representative,
):
    """# 2021-08-02 23:14:12
    cluster given sequences (given as a dictionary containing counts of unique sequences) using kmer of a given length
    """
    dict_cluster = dict()
    for seq in dict_seq_count:
        int_seq_count = dict_seq_count[seq]
        flag_assigned_to_cluster = False
        set_kmer = set(Generate_Kmer(seq, len_kmer))
        for name_cluster in dict_cluster:
            c = dict_cluster[name_cluster]
            n_overlap_kmer = len(c["set_kmer"].intersection(set_kmer))
            if int_min_n_overlap_kmer <= n_overlap_kmer:
                c["seq_count"][seq] = int_seq_count
                c["n_seq"] += int_seq_count
                n_seq = c["n_seq"]
                counter_kmer = COUNTER(
                    list(set_kmer) * int_seq_count, c["counter_kmer"]
                )  # update kmer count
                c["counter_kmer"] = counter_kmer
                c["set_kmer"] = set(
                    kmer
                    for kmer in counter_kmer
                    if counter_kmer[kmer] / n_seq
                    >= float_min_proportion_read_to_select_kmer_representative
                )
                flag_assigned_to_cluster = True
                break
        if not flag_assigned_to_cluster:
            c = dict()
            c["set_kmer"] = set_kmer
            c["n_seq"] = int_seq_count
            c["seq_count"] = dict()
            c["seq_count"][seq] = int_seq_count
            c["counter_kmer"] = COUNTER(list(set_kmer) * int_seq_count)
            dict_cluster[UUID()] = c
    return dict_cluster