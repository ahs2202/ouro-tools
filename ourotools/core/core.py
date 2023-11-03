# # %% CORE %%
# # import internal modules

"""
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
This part should be uncommented in core.py
"""

from . import SC
from . import SAM
from . import MAP
from . import SEQ
from . import STR
from . import biobookshelf as bk


"""
||||||||||||||||||||||||||||||||
"""

# import biobookshelf as bk
# import biobookshelf.STR as STR
# import biobookshelf.SEQ as SEQ
# import biobookshelf.SC as SC
# import biobookshelf.SAM as SAM
# import biobookshelf.MAP as MAP
# bk.Wide( 100 )

"""
This part should be uncommented in jupyter notebook
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
"""

import joblib  # for persistent, reference-counting-free memory
from typing import Union, List, Literal, Dict
import os
import pandas as pd
import numpy as np
import multiprocessing as mp
import math
import logging
from copy import deepcopy
import pickle
import time
import glob
import gzip  # to handle gzip file
import shutil  # for copying file
import base64  # for converting binary to text data (web application)
import json  # to read and write JSON file
import matplotlib.pyplot as plt
import scipy.sparse
import io
import pysam
import intervaltree
import ast

# prepare asynchronous operations
import asyncio
import nest_asyncio

nest_asyncio.apply()

from bitarray import bitarray  ## binary arrays

from tqdm import tqdm as progress_bar  # for progress bar

# from tqdm.autonotebook import tqdm  as progress_bar # for progress bar with jupyter notebook integration # not compatible with multi-processing

import argparse
import traceback
import os, sys, getopt
from io import StringIO, BytesIO
import time
import math
import mappy
import pkg_resources

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

# set logging format
logging.basicConfig(
    format="%(asctime)s [%(name)s] <%(levelname)s> (%(funcName)s) - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger("ouro-count")

# define version
_version_ = "0.0.2"
_ourotools_version_ = _version_
_last_modified_time_ = "2023-11-02 23:01:10"

str_release_note = [
    """
    # %% RELEASE NOTE %%
    # 2023-08-10 16:41:24 
    draft version of 'LongFilterNSplit' function completed
    
    # 2023-08-24 14:40:47 
    draft version of 'LongExtractBarcodeFromBAM' and 'LongCreateReferenceSizeDistribution' function completed.
    
    # 2023-08-30 23:51:34 
    draft version of 'LongExportNormalizedCountMatrix' function was completed. chromosome-level multiprocessing is used, and size-distribution-normalized count matrix can be exported for multiple size ranges of interest, along with the raw count matrix.

    # 2023-09-19 12:38:25 
    utility functions 'DeduplicateBAM' and 'SplitBAMs' were implemented. Currently, 'DeduplicateBAM' function selects the longest read with the same CB-UMI pair for each bucket (unique CB-UMI attachment location on the genome).
    
    # 2023-09-27 23:12:53 
    'LongExportNormalizedCountMatrix' function was modified for more accurate identification of specific transcript using re-alignment (TES matching, filtering alignment with excessive soft-clipping, and strand-specific assignment)
    Also, a typo in graph labels were fixed.

    # 2023-10-11 01:08:56 
    'LongExportNormalizedCountMatrix' function was modified for exporting count features excluding internal poly-A primed reads, which can be helpful for alternative splicing analysis. Most features, excluding variant information features, will be exported with/without internal polyA primed reads if possible.
    
    # 2023-11-02 23:00:23 
    'LongFilterNSplit' : Ouro-Enrich applied samples can be now processed with 'LongFilterNSplit' with 'flag_recover_original_molecule_before_self_circularization_and_digestion' flag turned on.
    ##### Future implementations #####

    """
]
"""
  
  .oooooo.   ooooo     ooo ooooooooo.     .oooooo.           ooooooooooooo   .oooooo.     .oooooo.   ooooo         .oooooo..o 
 d8P'  `Y8b  `888'     `8' `888   `Y88.  d8P'  `Y8b          8'   888   `8  d8P'  `Y8b   d8P'  `Y8b  `888'        d8P'    `Y8 
888      888  888       8   888   .d88' 888      888              888      888      888 888      888  888         Y88bo.      
888      888  888       8   888ooo88P'  888      888              888      888      888 888      888  888          `"Y8888o.  
888      888  888       8   888`88b.    888      888 8888888      888      888      888 888      888  888              `"Y88b 
`88b    d88'  `88.    .8'   888  `88b.  `88b    d88'              888      `88b    d88' `88b    d88'  888       o oo     .d8P 
 `Y8bood8P'     `YbodP'    o888o  o888o  `Y8bood8P'              o888o      `Y8bood8P'   `Y8bood8P'  o888ooooood8 8""88888P'  
 
 
 _______  __   __  ______    _______         _______  _______  _______  ___      _______ 
|       ||  | |  ||    _ |  |       |       |       ||       ||       ||   |    |       |
|   _   ||  | |  ||   | ||  |   _   | ____  |_     _||   _   ||   _   ||   |    |  _____|
|  | |  ||  |_|  ||   |_||_ |  | |  ||____|   |   |  |  | |  ||  | |  ||   |    | |_____ 
|  |_|  ||       ||    __  ||  |_|  |         |   |  |  |_|  ||  |_|  ||   |___ |_____  |
|       ||       ||   |  | ||       |         |   |  |       ||       ||       | _____| |
|_______||_______||___|  |_||_______|         |___|  |_______||_______||_______||_______|


 ______     __  __     ______     ______     ______   ______     ______     __         ______    
/\  __ \   /\ \/\ \   /\  == \   /\  __ \   /\__  _\ /\  __ \   /\  __ \   /\ \       /\  ___\   
\ \ \/\ \  \ \ \_\ \  \ \  __<   \ \ \/\ \  \/_/\ \/ \ \ \/\ \  \ \ \/\ \  \ \ \____  \ \___  \  
 \ \_____\  \ \_____\  \ \_\ \_\  \ \_____\    \ \_\  \ \_____\  \ \_____\  \ \_____\  \/\_____\ 
  \/_____/   \/_____/   \/_/ /_/   \/_____/     \/_/   \/_____/   \/_____/   \/_____/   \/_____/  


"""


str_description = """

 _______  __   __  ______    _______         _______  _______  _______  ___      _______ 
|       ||  | |  ||    _ |  |       |       |       ||       ||       ||   |    |       |
|   _   ||  | |  ||   | ||  |   _   | ____  |_     _||   _   ||   _   ||   |    |  _____|
|  | |  ||  |_|  ||   |_||_ |  | |  ||____|   |   |  |  | |  ||  | |  ||   |    | |_____ 
|  |_|  ||       ||    __  ||  |_|  |         |   |  |  |_|  ||  |_|  ||   |___ |_____  |
|       ||       ||   |  | ||       |         |   |  |       ||       ||       | _____| |
|_______||_______||___|  |_||_______|         |___|  |_______||_______||_______||_______|

                                                                                                 
A comprehensive toolkit for quality control and analysis of single-cell long-read RNA-seq data
"""

str_documentation = """


documentation_date = '2023-07-30 15:57:49'
"""

def _get_random_integer(int_num_possible_integers: int):
    """# 2023-01-09 18:18:37
    int_num_possible_integers : int # the number of integers to create. for example, when 'int_num_possible_integers' =  10, 0 ~ 9 will be randomly generated.
    """
    return int(np.floor(np.random.random() * int_num_possible_integers))

def __chromosome_name_remove_chr__(str_name_chrom):
    """# 2023-01-06 00:38:25
    remove 'chr' prefix from the chromosome name.
    handle chromosome name system with 'chr' prefix and those without 'chr' prefix
    """
    if "chr" == str_name_chrom[:3]:
        str_name_chrom = str_name_chrom[3:]
        if str_name_chrom == "M":  # handle mitochondrial genome
            return "MT"
        else:
            return str_name_chrom
    else:
        return str_name_chrom
    
def Call_Mutation_from_Read(
    path_file_input,
    path_file_bam,
    path_folder_output,
    path_folder_temp,
    path_file_filtered_mutation,
    name_ref,
    int_min_mapq_unique_mapped,
):
    """# 2021-08-29 16:55:23
    retrieve bases at the previously identified interesting sites (given by 'path_file_filtered_mutation') for each aligned read
    """
    str_uuid = bk.UUID()

    # define interger representation of the CIGAR operations used in BAM files
    int_cigarop_S = 4
    int_cigarop_H = 5

    df_mut = pd.read_csv(
        path_file_filtered_mutation, sep="\t"
    )  # read filtered mutation info
    df_mut.refname = df_mut.refname.astype(str).astype(
        object
    )  # convert refname datatype to string
    df_mut.refpos = df_mut.refpos - 1  # 1-based > 0-based coordinates
    dict_ref_site_to_refbase = dict(
        (tuple(arr[:2]), arr[-1])
        for arr in df_mut[["refname", "refpos", "refbase"]].values
    )  # compose a dictionary mapping interesting sites in the reference sequences from the called mutation information to the reference bases

    newfile = gzip.open(
        f"{path_folder_temp}{str_uuid}.call_mutation_from_read.{name_ref}.tsv.gz", "wb"
    )

    for (
        id_read_start,
        id_refname_start,
        int_refstart_start,
        int_index_chunk_start,
        id_read_end,
        id_refname_end,
        int_refstart_end,
        int_index_chunk_end,
    ) in pd.read_csv(
        path_file_input, sep="\t"
    ).values:  # retrieve start and end points in the given bam file
        flag_single_threaded = (
            id_refname_start == "single_thread"
        )  # flag indicating the program is run in a single-threaded mode
        id_refname_start = str(id_refname_start)
        id_refname_end = str(id_refname_end)
        id_read_start = str(id_read_start)
        id_read_end = str(id_read_end)
        flag_start_point_entered = False  # a flag indicating whether the starting point has been entered during the sequential reading of a given bam file
        if flag_single_threaded:
            flag_start_point_entered = True

        with pysam.AlignmentFile(path_file_bam, "rb") as samfile:
            for r in (
                samfile.fetch()
                if flag_single_threaded
                else samfile.fetch(contig=id_refname_start, start=int_refstart_start)
            ):
                refname = r.reference_name
                if (
                    r.mapq < int_min_mapq_unique_mapped
                ):  # skip read whose mapq is below 'int_min_mapq_unique_mapped'
                    continue
                seq = r.seq
                if (
                    seq is None
                ):  # skip multi-mapped reads (minimap2 skip sequence information for multi-mapped reads)
                    continue
                len_seq = len(seq)
                cigartuples = r.cigartuples
                if (
                    int_cigarop_H == cigartuples[0][0]
                    or int_cigarop_H == cigartuples[-1][0]
                ):  # skip hard-clipped reads
                    continue
                qname = r.qname
                # check whether the start point has been reached
                if not flag_start_point_entered and qname == id_read_start:
                    flag_start_point_entered = True
                if not flag_start_point_entered:
                    continue
                # exit if the end point is reached
                refstart, refend = r.reference_start, r.reference_end
                if (
                    refname == id_refname_end
                    and refstart == int_refstart_end
                    and qname == id_read_end
                ):
                    break

                # retrieve soft-clipped sequences
                flag_left_softclipped = int_cigarop_S == cigartuples[0][0]
                flag_right_softclipped = int_cigarop_S == cigartuples[-1][0]
                if not (
                    flag_left_softclipped and flag_right_softclipped
                ):  # skip reads that does not contain soft-clipped reads at both ends (adaptors not detected at least one end)
                    continue

                """ retrieve called base of the current read at the interesting genomic sites (where mutations were previously detected at the bulk level) """
                l_called_base = (
                    []
                )  # retrieve the list of called bast at the interesting genomic sites for the current read
                for int_pos, str_base, str_qual in SAM.Generate_Base_and_Qual(
                    r
                ):  # iterate each match/mismatched position and retrieve base and quality score # 0-based coordinates
                    t_ref_site = (refname, int_pos)
                    if (
                        t_ref_site in dict_ref_site_to_refbase
                    ):  # if the current loci is interesting loci
                        str_refbase = dict_ref_site_to_refbase[
                            t_ref_site
                        ]  # retrieve reference base at the reference site
                        l_called_base.append(
                            refname
                            + ":"
                            + str(int_pos + 1)
                            + ":"
                            + (
                                str_refbase
                                if str_refbase == str_base
                                else str_refbase + ">" + str_base
                            )
                        )  # 0-based coordinate > 1-based coordinate

                # write record when read spans an interesting reference site
                if len(l_called_base) > 0:
                    newfile.write(
                        (
                            "\t".join(
                                [
                                    qname,
                                    name_ref,
                                    str(refstart),
                                    ",".join(l_called_base),
                                ]
                            )
                            + "\n"
                        ).encode()
                    )

    newfile.close()

def Preprocess_and_Load_Annotations(
    path_folder_ref,
    path_file_gtf_genome: Union[str, None] = None,
    str_name_gtf_attr_for_id_gene: Union[str, None] = None,
    str_name_gtf_attr_for_name_gene: Union[str, None] = None,
    str_name_gtf_attr_for_id_transcript: Union[str, None] = None,
    str_name_gtf_attr_for_name_transcript: Union[str, None] = None,
    flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation: Union[
        bool, None
    ] = None,
    flag_does_not_make_gene_names_unique: Union[bool, None] = None,
    int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error: Union[
        int, None
    ] = None,
    path_file_tsv_repeatmasker_ucsc: Union[str, None] = None,
    l_repClass_repeatmasker_ucsc: Union[str, None] = None,
    int_min_length_repeatmasker_ucsc: Union[int, None] = None,
    path_file_gff_regulatory_element: Union[str, None] = None,
    str_name_gff_attr_id_regulatory_element: Union[str, None] = None,
    int_min_length_regulatory_element: Union[int, None] = None,
    int_bp_padding_regulatory_element_anno: Union[int, None] = None,
    path_file_fa_transcriptome: Union[str, None] = None,
    path_file_fa_genome: Union[str, None] = None,
    int_bp_padding_for_defining_promoter_from_transcript_start: Union[int, None] = None,
    flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome: Union[
        bool, None
    ] = False,
    dict_num_manager_processes_for_each_data_object: dict = {
        'dict_it_promoter' : 0,
        'dict_t_splice_junc_to_info_genome' : 0,
        'dict_it_exon' : 0,
        'dict_it_exon_transcriptome' : 3,
        'dict_it_splice_junc_transcriptome' : 3,
        'dict_it_splice_donor_and_acceptor_genome' : 3,
        'dict_it_rpmk' : 5,
        'dict_it_reg' : 3,
        'dict_fa_transcriptome' : 2,
    },
):
    """# 2023-07-26 19:52:00 
    Preprocess and Load annotation prior to the analysis from the given folder 'path_folder_ref'


    dict_num_manager_processes_for_each_data_object : dict = {
        'dict_it_promoter' : 0,
        'dict_t_splice_junc_to_info_genome' : 0,
        'dict_it_exon' : 0,
        'dict_it_exon_transcriptome' : 3,
        'dict_it_splice_junc_transcriptome' : 3,
        'dict_it_splice_donor_and_acceptor_genome' : 3,
        'dict_it_rpmk' : 5,
        'dict_it_reg' : 3,
        'dict_fa_transcriptome' : 2,
    }
    # the number of manager processes to use for each data object that will be shared across the forked processes. If 0 is given, no manager process will be used. Instead, the object will be directly accessed in the forked process, incurring memory bloating.
    # generally, it is better to use more number of manager processes for data object that are more frequently accessed. If increasing the number of manager processes does not improve performance, considering not using the manager process and accessing the object directly.
    # the expected size of bloated memory per process for each data object is given below.
    #
    #   'object name'                                       'the size of bloated memory per process'     'Manager Class Type'
    #   dict_it_exon_transcriptome                          1.617437 GB per process                      HostedDictIntervalTree
    #   dict_it_rpmk                                        1.452151 GB per process                      HostedDictIntervalTree
    #   dict_it_splice_junc_transcriptome                   1.381314 GB per process                      HostedDictIntervalTree
    #   dict_it_splice_donor_and_acceptor_genome            ???????? GB per process (not measured)       HostedDictIntervalTree
    #   dict_fa_transcriptome                               0.460438 GB per process                      HostedDict
    #   dict_it_exon                                        0.271540 GB per process                      HostedDictIntervalTree
    #   dict_it_reg                                         0.271540 GB per process                      HostedDictIntervalTree
    #   dict_t_splice_junc_to_info_genome                   0.188898 GB per process                      HostedDict
    #   dict_it_promoter                                    0.141673 GB per process                      HostedDictIntervalTree
    #   dict_fa_genome                                      0.082643 GB per process                      -
    #   dict_id_tx_to_id_gene                               0.070837 GB per process                      -
    #   dict_id_tx_to_name_tx                               0.070837 GB per process                      -
    #   dict_it_gene                                        0.059031 GB per process                      -
    #   dict_id_gene_to_l_id_tx                             0.059031 GB per process                      -
    #   dict_index_df_gtf_gene                              0.047224 GB per process                      -
    #   arr_data_df_gtf_gene                                0.047224 GB per process                      -
    #   dict_seqname_to_mask_gtf_reg                        0.035418 GB per process                      -
    #   dict_seqname_to_mask_gtf_intron_near_splice_site    0.035418 GB per process                      -
    #   dict_seqname_to_mask_gtf_rpmk_unfiltered            0.035418 GB per process                      -
    #   dict_seqname_to_mask_gtf_rpmk_filtered              0.035418 GB per process                      -
    #   dict_seqname_to_mask_gtf_exon                       0.035418 GB per process                      -
    #
    # if pre-loaded 'scidx' is given, this argument will be ignored.

    """
    # load ouro-count index for persistent access
    path_file_scidx = f"{path_folder_ref}scidx.pickle"
    if os.path.exists(path_file_scidx):
        logger.info(f"loading an existing index '{path_file_scidx}'")

        # load ouro-count index for persistent access
        scidx = joblib.load(path_file_scidx)  # load the scidx
    else:
        # create an output folder
        os.makedirs(path_folder_ref, exist_ok=True)

        # initialize ouro-index
        scidx = dict()

        """ load genome sequences """
        path_file_flag = f"{path_folder_ref}genome.fa.processing_completed.flag"
        path_file_pickle_dict_fa_genome = f"{path_folder_ref}dict_fa_genome.pickle"
        if not os.path.exists(path_file_flag):
            # read genome sequences
            scidx["dict_fa_genome"] = bk.FASTA_Read(
                path_file_fa_genome, header_split_at_space=True
            )
            bk.PICKLE_Write(path_file_pickle_dict_fa_genome, scidx["dict_fa_genome"])

            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            scidx["dict_fa_genome"] = bk.PICKLE_Read(
                path_file_pickle_dict_fa_genome
            )  # load genome
        logger.info("[Completed] the genome was loaded.")

        # retrieve sequence length of reference sequences from the BAM header
        dict_seqname_to_len_seq = dict(
            (e, len(scidx["dict_fa_genome"][e])) for e in scidx["dict_fa_genome"]
        )

        """
        Report annotations settings used in building the reference
        """

        """ load or export setting for reference annotations """
        path_file_json_ref_setting = f"{path_folder_ref}ref_setting.json"
        if os.path.exists(path_file_json_ref_setting):
            """load setting for reference annotations"""
            with open(path_file_json_ref_setting, "r") as file:
                dict_setting_ref = json.load(
                    file
                )  # override current program setting with previous program setting
            # parse settings update values in the local scope
            dict_seqname_to_len_seq = dict_setting_ref["dict_seqname_to_len_seq"]
        else:
            """export setting for reference annotations"""
            for e in [
                path_file_gtf_genome,
                str_name_gtf_attr_for_id_gene,
                str_name_gtf_attr_for_name_gene,
                str_name_gtf_attr_for_id_transcript,
                str_name_gtf_attr_for_name_transcript,
                flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation,
                flag_does_not_make_gene_names_unique,
                int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
                path_file_tsv_repeatmasker_ucsc,
                l_repClass_repeatmasker_ucsc,
                int_min_length_repeatmasker_ucsc,
                path_file_gff_regulatory_element,
                str_name_gff_attr_id_regulatory_element,
                int_min_length_regulatory_element,
                int_bp_padding_regulatory_element_anno,
                path_file_fa_transcriptome,
                path_file_fa_genome,
            ]:
                if e is None:
                    logger.error("required arguments were not given")
                    # return -1

            # record arguments used for the program (metadata)
            dict_setting_ref = {
                "dict_seqname_to_len_seq": dict_seqname_to_len_seq,
                "path_file_gtf_genome": path_file_gtf_genome,
                "str_name_gtf_attr_for_id_gene": str_name_gtf_attr_for_id_gene,
                "str_name_gtf_attr_for_name_gene": str_name_gtf_attr_for_name_gene,
                "str_name_gtf_attr_for_id_transcript": str_name_gtf_attr_for_id_transcript,
                "str_name_gtf_attr_for_name_transcript": str_name_gtf_attr_for_name_transcript,
                "flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation": flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation,
                "flag_does_not_make_gene_names_unique": flag_does_not_make_gene_names_unique,
                "int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error": int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
                "path_file_tsv_repeatmasker_ucsc": path_file_tsv_repeatmasker_ucsc,
                "l_repClass_repeatmasker_ucsc": l_repClass_repeatmasker_ucsc,
                "int_min_length_repeatmasker_ucsc": int_min_length_repeatmasker_ucsc,
                "path_file_gff_regulatory_element": path_file_gff_regulatory_element,
                "str_name_gff_attr_id_regulatory_element": str_name_gff_attr_id_regulatory_element,
                "int_min_length_regulatory_element": int_min_length_regulatory_element,
                "int_bp_padding_regulatory_element_anno": int_bp_padding_regulatory_element_anno,
                "path_file_fa_transcriptome": path_file_fa_transcriptome,
                "path_file_fa_genome": path_file_fa_genome,
                "int_bp_padding_for_defining_promoter_from_transcript_start": int_bp_padding_for_defining_promoter_from_transcript_start,
                "flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome": flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome,
            }
            with open(path_file_json_ref_setting, "w") as newfile:
                json.dump(dict_setting_ref, newfile)
        """ preprocess UCSC repeatmasker annotation """
        path_file_flag = f"{path_folder_ref}repeatmasker.processing_completed.flag"
        path_file_pickle_df_rpmk = f"{path_folder_ref}df_rpmk.pickle"
        path_file_pickle_df_rpmk_unfiltered = (
            f"{path_folder_ref}df_rpmk_unfiltered.pickle"
        )
        if not os.path.exists(path_file_flag):

            def __Preprocess_UCSC_repeatmasker_anno__(
                path_file_tsv_repeatmasker_ucsc,
                l_repClass_repeatmasker_ucsc,
                int_min_length_repeatmasker_ucsc,
            ):
                """preprocess UCSC repeatmasker annotation"""
                # compose empty rmpk dataframe
                df_rmpk_empty = pd.DataFrame(
                    columns=[
                        "seqname",
                        "source",
                        "feature",
                        "start",
                        "end",
                        "score",
                        "strand",
                        "frame",
                        "attribute",
                        "repName",
                        "id_repeat",
                    ]
                )
                # if 'path_file_tsv_repeatmasker_ucsc' is given
                if path_file_tsv_repeatmasker_ucsc is not None:
                    # read repeatmasker table
                    df_rpmk = pd.read_csv(path_file_tsv_repeatmasker_ucsc, sep="\t")

                    # if 'genoName' contains 'chr', remove 'chr'
                    df_rpmk["genoName"] = (
                        list(e[3:] for e in df_rpmk.genoName.values)
                        if df_rpmk.genoName.values[0][:3] == "chr"
                        else df_rpmk.genoName
                    )

                    # rename columns to match that of df_gtf
                    df_rpmk.rename(
                        columns={
                            "genoName": "seqname",
                            "genoStart": "start",
                            "genoEnd": "end",
                        },
                        inplace=True,
                    )
                    # create a column to match that of df_gtf
                    df_rpmk["score"] = df_rpmk["swScore"]
                    df_rpmk["feature"] = "gene"
                    df_rpmk["source"] = "repeatmasker_ucsc"
                    df_rpmk["frame"] = "."
                    df_rpmk["attribute"] = ""

                    df_rpmk_unfiltered = df_rpmk[
                        [
                            "seqname",
                            "source",
                            "feature",
                            "start",
                            "end",
                            "score",
                            "strand",
                            "frame",
                            "attribute",
                            "repName",
                            "repClass",
                        ]
                    ]
                    """ filtering """
                    # select only given repClass entries from the UCSC repeatmasker table
                    if (
                        len(l_repClass_repeatmasker_ucsc) > 0
                    ):  # if valid list of repClass is given, select only the given repClass entries
                        df_rpmk = bk.PD_Select(
                            df_rpmk, repClass=l_repClass_repeatmasker_ucsc
                        )

                    # filtering out repeat element with its length shorter than the minimum threshold
                    df_rpmk["int_len_repeat"] = (
                        df_rpmk.end - df_rpmk.start + 1
                    )  # retrieve length of repeat element # 1-based
                    df_rpmk = df_rpmk[
                        df_rpmk.int_len_repeat >= int_min_length_repeatmasker_ucsc
                    ]

                    # discard unnecessary columns
                    df_rpmk = df_rpmk[
                        [
                            "seqname",
                            "source",
                            "feature",
                            "start",
                            "end",
                            "score",
                            "strand",
                            "frame",
                            "attribute",
                            "repName",
                            "repClass",
                        ]
                    ]

                    # compose identifier of repeatmasker entries
                    def __get_string_object_series__(s):
                        return s.astype(str).astype(object)

                    df_rpmk["id_repeat"] = (
                        "repeatmasker_ucsc|repClass="
                        + df_rpmk.repClass
                        + "|repName="
                        + df_rpmk.repName
                        + "|pos="
                        + __get_string_object_series__(df_rpmk.seqname)
                        + ":"
                        + __get_string_object_series__(df_rpmk.start)
                        + "-"
                        + __get_string_object_series__(df_rpmk.end)
                    )
                else:
                    df_rpmk, df_rpmk_unfiltered = df_rmpk_empty, df_rmpk_empty
                return df_rpmk, df_rpmk_unfiltered

            df_rpmk, df_rpmk_unfiltered = __Preprocess_UCSC_repeatmasker_anno__(
                path_file_tsv_repeatmasker_ucsc,
                l_repClass_repeatmasker_ucsc,
                int_min_length_repeatmasker_ucsc,
            )
            # write filtered repeatmasker GTF as files
            bk.GTF_Write(df_rpmk, f"{path_folder_ref}repeatmasker_ucsc.filtered.gtf.gz")
            bk.GTF_Write(
                df_rpmk_unfiltered,
                f"{path_folder_ref}repeatmasker_ucsc.unfiltered.gtf.gz",
            )
            bk.PICKLE_Write(path_file_pickle_df_rpmk, df_rpmk)
            bk.PICKLE_Write(path_file_pickle_df_rpmk_unfiltered, df_rpmk_unfiltered)
            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            df_rpmk = bk.PICKLE_Read(path_file_pickle_df_rpmk)
            df_rpmk_unfiltered = bk.PICKLE_Read(path_file_pickle_df_rpmk_unfiltered)
        logger.info("[Completed] Repeatmasker annotations were loaded.")

        """ preprocess Regulatory Build annotations """
        path_file_flag = (
            f"{path_folder_ref}regulatory_element.processing_completed.flag"
        )
        path_file_pickle_df_gtf_reg = f"{path_folder_ref}df_gtf_reg.pickle"
        if not os.path.exists(path_file_flag):

            def __Preprocess_Regulatory_Annotations__(
                path_file_gff_regulatory_element,
                str_name_gff_attr_id_regulatory_element,
                int_min_length_regulatory_element,
                dict_seqname_to_len_seq,
                int_bp_padding_regulatory_element_anno=None,
            ):
                """preprocess regulatory element annotations"""
                # compose empty gtf dataframe
                df_gtf_reg_empty = pd.DataFrame(
                    columns=[
                        "seqname",
                        "source",
                        "feature",
                        "start",
                        "end",
                        "score",
                        "strand",
                        "frame",
                        "attribute",
                        str_name_gff_attr_id_regulatory_element,
                        "id_regulatory_element",
                    ]
                )
                # if 'path_file_tsv_repeatmasker_ucsc' is given
                if path_file_gff_regulatory_element is not None:
                    # read repeatmasker table
                    df_gtf_reg = bk.GTF_Read(
                        path_file_gff_regulatory_element,
                        flag_gtf_format=False,
                        remove_chr_from_seqname=True,
                    )

                    # filtering out repeat element with its length shorter than the minimum threshold
                    df_gtf_reg["int_len_regulatory_element"] = (
                        df_gtf_reg.end - df_gtf_reg.start + 1
                    )  # retrieve length of repeat element
                    df_gtf_reg = df_gtf_reg[
                        df_gtf_reg.int_len_regulatory_element
                        >= int_min_length_regulatory_element
                    ]

                    """ apply padding to regulatory element annotations """
                    if int_bp_padding_regulatory_element_anno is not None:
                        df_gtf_reg.start = (
                            df_gtf_reg.start - int_bp_padding_regulatory_element_anno
                        )
                        df_gtf_reg.start[
                            df_gtf_reg.start < 0
                        ] = 0  # handle invalid coordinates outside boundaries after applying padding
                        df_gtf_reg.end = (
                            df_gtf_reg.end + int_bp_padding_regulatory_element_anno
                        )
                        # handle invalid coordinates outside boundaries after applying padding
                        l_end = []
                        for seqname, end in df_gtf_reg[["seqname", "end"]].values:
                            if seqname in dict_seqname_to_len_seq:
                                end = min(end, dict_seqname_to_len_seq[seqname])
                            l_end.append(end)
                        df_gtf_reg["end"] = l_end

                    # compose identifier of repeatmasker entries
                    def __get_string_object_series__(s):
                        return s.astype(str).astype(object)

                    df_gtf_reg["id_regulatory_element"] = (
                        "regulatory_element|ID="
                        + df_gtf_reg[str_name_gff_attr_id_regulatory_element]
                        + "|pos="
                        + __get_string_object_series__(df_gtf_reg.seqname)
                        + ":"
                        + __get_string_object_series__(df_gtf_reg.start)
                        + "-"
                        + __get_string_object_series__(df_gtf_reg.end)
                    )
                else:
                    df_gtf_reg = df_gtf_reg_empty
                return df_gtf_reg

            df_gtf_reg = __Preprocess_Regulatory_Annotations__(
                path_file_gff_regulatory_element,
                str_name_gff_attr_id_regulatory_element,
                int_min_length_regulatory_element,
                dict_seqname_to_len_seq,
                int_bp_padding_regulatory_element_anno,
            )
            # write filtered regulatory elements as a file
            bk.GTF_Write(df_gtf_reg, f"{path_folder_ref}regulatory_element.gtf.gz")
            bk.PICKLE_Write(path_file_pickle_df_gtf_reg, df_gtf_reg)
            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            df_gtf_reg = bk.PICKLE_Read(path_file_pickle_df_gtf_reg)
        logger.info("[Completed] Regulatory annotations were loaded.")

        """ 
        pre-process gene annotations 
        """
        path_file_flag = f"{path_folder_ref}gtf_gene.processing_completed.flag"
        path_file_pickle_df_gtf = f"{path_folder_ref}df_gtf.pickle"
        path_file_pickle_df_gtf_gene = f"{path_folder_ref}df_gtf_gene.pickle"
        path_file_pickle_df_gtf_transcript = (
            f"{path_folder_ref}df_gtf_transcript.pickle"
        )
        path_file_pickle_df_gtf_promoter = f"{path_folder_ref}df_gtf_promoter.pickle"
        if not os.path.exists(path_file_flag):
            # only save gtf records of genes & transcripts for faster GTF file loading
            df_gtf = bk.GTF_Read(path_file_gtf_genome, parse_attr=True)
            df_gtf_gene = bk.PD_Select(df_gtf, feature="gene")
            df_gtf_gene.dropna(
                subset=[str_name_gtf_attr_for_id_gene], inplace=True
            )  # drop entries without 'id_gene'

            """
            merge & make name_gene unique for genes whose name is not unique
            """
            # retrieve the list of gene_name of genes whose name is not unique
            l_name_gene_not_unique = bk.LIST_COUNT(
                df_gtf_gene[str_name_gtf_attr_for_name_gene], duplicate_filter=2
            ).index.values

            def __combine_id_gene__(e_reduced, e_new):
                return f"{e_reduced};{e_new}"

            l_l = []
            for name_gene in l_name_gene_not_unique:
                df = bk.PD_Select(
                    df_gtf_gene, **{str_name_gtf_attr_for_name_gene: name_gene}
                )
                df = df.sort_values(
                    ["seqname", "start", "end", "strand"]
                )  # sort annotations in the order to provide consistently (non-randomly) modified annotations
                dict_it_for_a_gene = dict()
                for seqname, start, end, strand, id_gene in df[
                    ["seqname", "start", "end", "strand", str_name_gtf_attr_for_id_gene]
                ].values:  # 1-based coordinates
                    start -= 1  # 1->0 based coordinates
                    if (seqname, strand) not in dict_it_for_a_gene:
                        dict_it_for_a_gene[
                            seqname, strand
                        ] = intervaltree.IntervalTree()
                    dict_it_for_a_gene[seqname, strand][start:end] = id_gene
                """ merge overlapping gene annotations with the same gene_name """
                if (
                    not flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation
                ):
                    for seqname, strand in dict_it_for_a_gene:
                        dict_it_for_a_gene[seqname, strand].merge_overlaps(
                            data_reducer=__combine_id_gene__
                        )

                """ rename genes with the same gene_name so that the gene_name becomes unique  """
                int_n_duplicated_gene_names = 0
                for seqname, strand in dict_it_for_a_gene:
                    for start, end, id_gene in dict_it_for_a_gene[
                        seqname, strand
                    ]:  # 0-based
                        l_l.append(
                            [
                                seqname,
                                start,
                                end,
                                strand,
                                id_gene,
                                name_gene
                                if flag_does_not_make_gene_names_unique
                                or int_n_duplicated_gene_names == 0
                                else f"{name_gene}_dup{int_n_duplicated_gene_names}",
                            ]
                        )
                        int_n_duplicated_gene_names += 1
            df_gtf_gene_whose_name_was_not_unique = pd.DataFrame(
                l_l,
                columns=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    str_name_gtf_attr_for_id_gene,
                    str_name_gtf_attr_for_name_gene,
                ],
            )
            # convert some columns' dtype to string
            for name_col in [
                "seqname",
                str_name_gtf_attr_for_id_gene,
                str_name_gtf_attr_for_name_gene,
            ]:
                df_gtf_gene_whose_name_was_not_unique[name_col] = (
                    df_gtf_gene_whose_name_was_not_unique[name_col]
                    .astype(str)
                    .astype(object)
                )
            # fill out uniform values
            for name_col, uniform_value in zip(
                ["source", "score", "frame", "feature"],
                ["processed_by_ouro", ".", ".", "gene"],
            ):
                df_gtf_gene_whose_name_was_not_unique[name_col] = uniform_value
            # combine and construct new gene annotations
            df_gtf_gene = pd.concat(
                [
                    bk.PD_Select(
                        df_gtf_gene,
                        **{str_name_gtf_attr_for_name_gene: l_name_gene_not_unique},
                        deselect=True,
                    ),
                    df_gtf_gene_whose_name_was_not_unique,
                ]
            )  # replace gtf record of genes whose gene_name is not unique with processed gtf records whose gene_name was made unique or annotations were merged

            """ 
            pre-process transcript annotations 
            """
            """
            apply and link modified gene annotation to all gtf annotations
            """
            # map original gene_id to new annotation
            # acquire mapping
            dict_id_gene_original_to_new_anno = dict()
            for id_gene, name_gene in df_gtf_gene_whose_name_was_not_unique[
                [str_name_gtf_attr_for_id_gene, str_name_gtf_attr_for_name_gene]
            ].values:
                for id_gene_original in id_gene.split(";"):
                    dict_id_gene_original_to_new_anno[id_gene_original] = (
                        id_gene,
                        name_gene,
                    )
            # perform mapping original gene_id to new annotation for the entire GTF records
            l_id_gene, l_name_gene = [], []
            for id_gene_original, name_gene_original in df_gtf[
                [str_name_gtf_attr_for_id_gene, str_name_gtf_attr_for_name_gene]
            ].values:
                if id_gene_original in dict_id_gene_original_to_new_anno:
                    id_gene, name_gene = dict_id_gene_original_to_new_anno[
                        id_gene_original
                    ]
                    l_id_gene.append(id_gene), l_name_gene.append(name_gene)
                else:
                    l_id_gene.append(id_gene_original), l_name_gene.append(
                        name_gene_original
                    )
            df_gtf[str_name_gtf_attr_for_id_gene] = l_id_gene
            df_gtf[str_name_gtf_attr_for_name_gene] = l_name_gene

            # retrieve GTF file for transcripts
            df_gtf_transcript = bk.PD_Select(df_gtf, feature="transcript")

            """ write data as GTF files and pickle files """
            df_gtf = df_gtf[
                [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    str_name_gtf_attr_for_id_gene,
                    str_name_gtf_attr_for_name_gene,
                    str_name_gtf_attr_for_id_transcript,
                    str_name_gtf_attr_for_name_transcript,
                ]
            ]  # drop unnecessary columns
            path_file_gtf_genome_processed = f"{path_folder_ref}gtf.genome.gtf.gz"
            bk.GTF_Write(
                df_gtf, path_file_gtf_genome_processed, flag_update_attribute=True
            )
            bk.PICKLE_Write(path_file_pickle_df_gtf, df_gtf)

            path_file_gtf_gene = f"{path_folder_ref}gtf.gene.genome.gtf"
            bk.GTF_Write(df_gtf_gene, path_file_gtf_gene)
            bk.PICKLE_Write(path_file_pickle_df_gtf_gene, df_gtf_gene)

            path_file_gtf_transcript = f"{path_folder_ref}gtf.transcript.genome.gtf"
            bk.GTF_Write(df_gtf_transcript, path_file_gtf_transcript)
            bk.PICKLE_Write(path_file_pickle_df_gtf_transcript, df_gtf_transcript)
            """
            Prepare promoter annotations from processed transcript annotations
            """
            # make sure each record in 'df_gtf_transcript' contains exactly one record for each unique transcript
            assert len(df_gtf_transcript) == len(
                df_gtf_transcript.transcript_id.unique()
            )
            # process transcript with + strand
            s_tx_positive_strand_promoter_end_in_1basedcoord = (
                bk.PD_Select(df_gtf_transcript, strand="+")[
                    ["start", str_name_gtf_attr_for_id_transcript]
                ]
                .groupby(str_name_gtf_attr_for_id_transcript)
                .min()["start"]
                - 1
            )  # retrieve promoter end sites for '+' strand transcripts
            s_tx_positive_strand_promoter_start_in_1basedcoord = (
                s_tx_positive_strand_promoter_end_in_1basedcoord
                - int_bp_padding_for_defining_promoter_from_transcript_start
                + 1
            )
            # process transcript with - strand
            s_tx_negative_strand_promoter_start_in_1basedcoord = (
                bk.PD_Select(df_gtf_transcript, strand="-")[
                    ["end", str_name_gtf_attr_for_id_transcript]
                ]
                .groupby(str_name_gtf_attr_for_id_transcript)
                .max()["end"]
                + 1
            )  # retrieve promoter start sites for '-' strand transcripts
            s_tx_negative_strand_promoter_end_in_1basedcoord = (
                s_tx_negative_strand_promoter_start_in_1basedcoord
                + int_bp_padding_for_defining_promoter_from_transcript_start
                - 1
            )
            # combine information
            s_tx_promoter_end_in_1basedcoord = pd.concat(
                [
                    s_tx_positive_strand_promoter_end_in_1basedcoord,
                    s_tx_negative_strand_promoter_end_in_1basedcoord,
                ]
            )
            s_tx_promoter_start_in_1basedcoord = pd.concat(
                [
                    s_tx_positive_strand_promoter_start_in_1basedcoord,
                    s_tx_negative_strand_promoter_start_in_1basedcoord,
                ]
            )
            # correct invalid (0 or negative) promoter start position
            s_tx_promoter_start_in_1basedcoord[
                s_tx_promoter_start_in_1basedcoord < 1
            ] = 1

            # compose 'df_gtf_promoter'
            df_gtf_promoter = pd.DataFrame(
                {
                    "start": s_tx_promoter_start_in_1basedcoord,
                    "end": s_tx_promoter_end_in_1basedcoord,
                }
            )
            df = deepcopy(df_gtf_transcript)
            df.drop(
                columns=["source", "feature", "start", "end", "attribute"], inplace=True
            )
            df.set_index(str_name_gtf_attr_for_id_transcript, inplace=True)
            df_gtf_promoter = df_gtf_promoter.join(
                df, how="left"
            )  # retrieve transcript annotation from 'df_gtf_transcript'
            df_gtf_promoter["source"] = "processed_by_ouro"
            df_gtf_promoter["feature"] = "promoter"
            df_gtf_promoter.reset_index(drop=False, inplace=True)

            # reorder columns to be compatible with downstream applications
            l_col_essential_for_gtf = [
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
            ]
            df_gtf_promoter = df_gtf_promoter[
                l_col_essential_for_gtf
                + list(
                    col
                    for col in df_gtf_promoter.columns.values
                    if col not in l_col_essential_for_gtf
                )
            ]

            # retrieve chromosome length
            df_gtf_promoter["chrom_length"] = df_gtf_promoter.seqname.apply(
                MAP.Map(dict_seqname_to_len_seq).a2b
            )
            # correct invalid (larger than chromosome length) promoter end position
            df_gtf_promoter.loc[
                df_gtf_promoter.end > df_gtf_promoter.chrom_length, "end"
            ] = df_gtf_promoter.loc[
                df_gtf_promoter.end > df_gtf_promoter.chrom_length, "chrom_length"
            ]
            df_gtf_promoter.drop(columns=["chrom_length"], inplace=True)
            df_gtf_promoter.loc[
                df_gtf_promoter.end < df_gtf_promoter.start, "end"
            ] = df_gtf_promoter.start  # handle the invalid intervals

            # write an output file
            path_file_gtf_promoter = f"{path_folder_ref}gtf.promoter.genome.gtf"
            bk.GTF_Write(
                df_gtf_promoter,
                path_file_gtf_promoter,
                flag_update_attribute=True,
                flag_filetype_is_gff3=False,
            )
            bk.PICKLE_Write(path_file_pickle_df_gtf_promoter, df_gtf_promoter)
            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            """read data"""
            path_file_gtf_genome_processed = f"{path_folder_ref}gtf.genome.gtf.gz"
            path_file_gtf_gene = f"{path_folder_ref}gtf.gene.genome.gtf"
            path_file_gtf_transcript = f"{path_folder_ref}gtf.transcript.genome.gtf"
            path_file_gtf_promoter = f"{path_folder_ref}gtf.promoter.genome.gtf"

            df_gtf = bk.PICKLE_Read(path_file_pickle_df_gtf)
            df_gtf_gene = bk.PICKLE_Read(path_file_pickle_df_gtf_gene)
            df_gtf_transcript = bk.PICKLE_Read(path_file_pickle_df_gtf_transcript)
            df_gtf_promoter = bk.PICKLE_Read(path_file_pickle_df_gtf_promoter)
        logger.info("[Completed] Gene annotation pre-processing was completed.")

        """ build interval trees of the annotations """
        path_file_flag = f"{path_folder_ref}interval_tree.export_completed.flag"
        path_file_pickle_dict_it_gene = f"{path_folder_ref}dict_it_gene.pickle"
        path_file_pickle_dict_it_exon = f"{path_folder_ref}dict_it_exon.pickle"
        path_file_pickle_dict_it_promoter = f"{path_folder_ref}dict_it_promoter.pickle"
        path_file_pickle_dict_it_rpmk = f"{path_folder_ref}dict_it_rpmk.pickle"
        path_file_pickle_dict_it_reg = f"{path_folder_ref}dict_it_reg.pickle"
        path_file_pickle_df_gtf_exon = f"{path_folder_ref}df_gtf_exon.pickle"
        if not os.path.exists(path_file_flag):
            """load gene and exon annotation as an interval tree"""
            # at gene body level
            scidx["dict_it_gene"] = bk.GTF_Interval_Tree(
                df_gtf_gene, feature="gene", value=str_name_gtf_attr_for_id_gene
            )  # load gene annotations
            bk.PICKLE_Write(path_file_pickle_dict_it_gene, scidx["dict_it_gene"])

            # at exon level (drop duplicated exon annotation for each unique id_gene)
            df_gtf_exon = bk.PD_Select(df_gtf, feature="exon")
            df_gtf_exon.sort_values(
                ["seqname", str_name_gtf_attr_for_id_transcript, "start"], inplace=True
            )  # sort exons by start site position for downstream analysis
            bk.PICKLE_Write(path_file_pickle_df_gtf_exon, df_gtf_exon)
            dict_it_exon = bk.GTF_Interval_Tree(
                df_gtf_exon.drop_duplicates(
                    subset=["seqname", "start", "end", str_name_gtf_attr_for_id_gene]
                ),
                feature="exon",
                value=str_name_gtf_attr_for_id_gene,
            )  # load drop duplicates across transcripts, but does not drop duplicates across genes so that duplicated exons belonging to two genes can be labeled for each gene separately
            bk.PICKLE_Write(path_file_pickle_dict_it_exon, dict_it_exon)

            """ load promoter annotations as an interval tree """
            dict_it_promoter = bk.GTF_Interval_Tree(
                df_gtf_promoter.drop_duplicates(
                    subset=[
                        "seqname",
                        "start",
                        "end",
                        "strand",
                        str_name_gtf_attr_for_id_gene,
                    ]
                ),
                feature="promoter",
                value=[str_name_gtf_attr_for_id_gene, "strand"],
            )  # load gene annotations # Some promoter sequences are shared between different genes (ensembl)
            bk.PICKLE_Write(path_file_pickle_dict_it_promoter, dict_it_promoter)

            """ load repeatmasker annotations as an interval tree """
            dict_it_rpmk = bk.GTF_Interval_Tree(
                df_rpmk, feature="gene", value="id_repeat"
            )  # load repeatmasker annotations # filtered repeatmasker annotations
            bk.PICKLE_Write(path_file_pickle_dict_it_rpmk, dict_it_rpmk)

            """ load regulatory element as an interval tree """
            dict_it_reg = bk.GTF_Interval_Tree(
                df_gtf_reg, feature=None, value="id_regulatory_element"
            )  # load regulatory annotations # use all features
            bk.PICKLE_Write(path_file_pickle_dict_it_reg, dict_it_reg)

            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            # read interval trees
            scidx["dict_it_gene"] = bk.PICKLE_Read(path_file_pickle_dict_it_gene)
            df_gtf_exon = bk.PICKLE_Read(path_file_pickle_df_gtf_exon)
        logger.info(
            "[Completed] the construction of the interval tree of the gene bodies, exons, promoters, repeatmasker annotations, and regulatory annotations was completed."
        )

        """ export masks """
        path_file_flag = f"{path_folder_ref}mask.export_completed.flag"
        path_file_gtf_rpmk_filtered = (
            f"{path_folder_ref}repeatmasker_ucsc.filtered.gtf.gz"
        )
        path_file_gtf_rpmk_unfiltered = (
            f"{path_folder_ref}repeatmasker_ucsc.unfiltered.gtf.gz"
        )
        path_file_gtf_reg = f"{path_folder_ref}regulatory_element.gtf.gz"
        path_folder_mask_gtf_exon = f"{path_folder_ref}mask.gtf.exon.genome/"
        path_folder_mask_gtf_rpmk_filtered = (
            f"{path_folder_ref}mask.gtf.repeatmasker_ucsc.filtered/"
        )
        path_folder_mask_gtf_rpmk_unfiltered = (
            f"{path_folder_ref}mask.gtf.repeatmasker_ucsc.unfiltered/"
        )
        path_folder_mask_gtf_reg = f"{path_folder_ref}mask.gtf.regulatory_element/"
        if not os.path.exists(path_file_flag):
            """build masks of exons for filtering intronic reads"""
            scidx["dict_seqname_to_mask_gtf_exon"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                df_gtf=df_gtf,
                str_feature="exon",
                remove_chr_from_seqname=True,
                path_folder_output=path_folder_mask_gtf_exon,
            )

            """ build masks of repeatmasker elements (filtered & unfiltered) for flagging reads significantly overlaps with a repeat element """
            scidx["dict_seqname_to_mask_gtf_rpmk_filtered"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                df_gtf=df_rpmk,
                str_feature="gene",
                remove_chr_from_seqname=True,
                path_folder_output=path_folder_mask_gtf_rpmk_filtered,
            )  # filtered
            scidx["dict_seqname_to_mask_gtf_rpmk_unfiltered"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                df_gtf=df_rpmk_unfiltered,
                str_feature="gene",
                remove_chr_from_seqname=True,
                path_folder_output=path_folder_mask_gtf_rpmk_unfiltered,
            )  # unfiltered

            """ build masks of regulatory element for flagging reads significantly overlaps with a regulatory element """
            scidx["dict_seqname_to_mask_gtf_reg"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                df_gtf=df_gtf_reg,
                str_feature=None,
                remove_chr_from_seqname=True,
                path_folder_output=path_folder_mask_gtf_reg,
            )  # use all features

            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            # load masks
            scidx["dict_seqname_to_mask_gtf_exon"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                path_folder_output=path_folder_mask_gtf_exon,
            )
            scidx["dict_seqname_to_mask_gtf_rpmk_filtered"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                path_folder_output=path_folder_mask_gtf_rpmk_filtered,
            )
            scidx["dict_seqname_to_mask_gtf_rpmk_unfiltered"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                path_folder_output=path_folder_mask_gtf_rpmk_unfiltered,
            )
            scidx["dict_seqname_to_mask_gtf_reg"] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                path_folder_output=path_folder_mask_gtf_reg,
            )

        logger.info(
            "[Completed] the construction of masks of exons, filtered repeatmasker regions, unfiltered repeatmasker regions, and regulatory annotations was completed."
        )

        """ load transcriptome sequences and retrieve information about transcript from the input transcriptome fasta header """
        path_file_flag = f"{path_folder_ref}transcriptome.fa.processing_completed.flag"
        path_file_pickle_dict_fa_transcriptome = (
            f"{path_folder_ref}dict_fa_transcriptome.pickle"
        )
        if not os.path.exists(path_file_flag):
            dict_fa_transcriptome = bk.FASTA_Read(path_file_fa_transcriptome)
            """ remove the version info. from id_transcript if 'flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome' is False """
            if (
                not flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome
            ):
                dict_fa_transcriptome = dict(
                    (k.split(" ", 1)[0].rsplit(".", 1)[0], dict_fa_transcriptome[k]) # remove version info from transcriptome
                    for k in dict_fa_transcriptome
                )
            else :
                dict_fa_transcriptome = dict(
                    (k.split(" ", 1)[0], dict_fa_transcriptome[k]) # does not remove the version info from transcriptome
                    for k in dict_fa_transcriptome
                )
            bk.FASTA_Write(
                f"{path_folder_ref}transcriptome.fa.gz",
                dict_fasta=dict_fa_transcriptome,
            )
            bk.PICKLE_Write(
                path_file_pickle_dict_fa_transcriptome, dict_fa_transcriptome
            )
            with open(path_file_flag, "w") as file:
                file.write("completed")
        logger.info("[Completed] the transcriptome was loaded.")

        """ build masks of intronic regions near exons for filtering false positive mutations in the intron near the splice site from splice-site detection errors """
        """ retrieve splice junctions """
        path_file_flag = f"{path_folder_ref}mask.gtf.intron_near_splice_site.genome.processing_completed.flag"
        path_folder_mask_gtf_intron_near_splice_site = (
            f"{path_folder_ref}mask.gtf.intron_near_splice_site.genome/"
        )
        path_file_pickle_dict_t_splice_junc_to_info_genome = (
            f"{path_folder_ref}dict_t_splice_junc_to_info_genome.pickle"
        )
        path_file_gtf_exon_transcriptome = (
            f"{path_folder_ref}gtf.exon.transcriptome.gtf"
        )
        path_file_gtf_splice_junc_transcriptome = (
            f"{path_folder_ref}gtf.splice_junc.transcriptome.gtf"
        )
        path_file_gtf_splice_donor_and_acceptor_genome = (
            f"{path_folder_ref}gtf.splice_donor_and_acceptor.genome.gtf"
        )
        path_file_pickle_dict_it_exon_transcriptome = (
            f"{path_folder_ref}dict_it_exon_transcriptome.pickle"
        )
        path_file_pickle_dict_it_splice_junc_transcriptome = (
            f"{path_folder_ref}dict_it_splice_junc_transcriptome.pickle"
        )
        path_file_pickle_dict_it_splice_donor_and_acceptor_genome = (
            f"{path_folder_ref}dict_it_splice_donor_and_acceptor_genome.pickle"
        )
        if not os.path.exists(path_file_flag):
            arr_df = df_gtf_exon[
                [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    str_name_gtf_attr_for_id_gene,
                ]
            ].values
            dict_index = bk.DF_Build_Index_Using_Dictionary(
                df_gtf_exon, l_col_for_index=str_name_gtf_attr_for_id_transcript
            )

            # label intronic regions where potential splice site detection error can occur based on the given window size 'int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error'
            l_l_intronic_region_near_splice_site = []
            l_l_df_gtf_exon_transcriptome = (
                []
            )  # a list of list that will contain records of GTF of transcriptome
            l_l_df_gtf_splice_junc_transcriptome = (
                []
            )  # a list of list that will contain splicing junction records of GTF of transcriptome
            l_df_gtf_splice_donor_and_acceptor_genome = (
                []
            )  # a list of list that will contain splicing donor and acceptor records of GTF of genome
            # collect splice site information
            dict_t_splice_junc_to_info_genome = dict()
            for id_tx in dict_index:
                arr_df_subset = arr_df[dict_index[id_tx]]
                """ identify single-exon transcript """
                int_n_exon = arr_df_subset.shape[0]

                """
                Retrieve splice junction from genomic information
                """
                """ retrieve records of intronic regions near exon-exon splice sites """
                seqname = arr_df_subset[0][0]  # retrieve seqname
                strand = arr_df_subset[0][6]  # retrieve strand
                id_gene = arr_df_subset[0][8]  # retrieve id_gene
                if (
                    int_n_exon > 1
                ):  # skip identification of splice junctions for single-exon transcripts
                    for left_exon_end, right_exon_start in (
                        arr_df_subset[:, 3:5].ravel()[1:-1].reshape((int_n_exon - 1, 2))
                    ):  # retrieve exon-exon junctions
                        l_l_intronic_region_near_splice_site.append(
                            [
                                seqname,
                                left_exon_end + 1,
                                left_exon_end
                                + int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
                            ]
                        )  # 1-based coordinates
                        l_l_intronic_region_near_splice_site.append(
                            [
                                seqname,
                                right_exon_start
                                - int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
                                right_exon_start - 1,
                            ]
                        )  # 1-based coordinates

                        """ retrieve a tuple representing a splice junction and collect information about the splicing junc """
                        t_splice_junction = (
                            seqname,
                            left_exon_end, # 1>0-based
                            right_exon_start - 1, # 1>0-based
                        )  # define a tuple representing a splice junction (a region of spliced out intron) # 0-based coordinates
                        
                        
                        ''' add a record to 'dict_t_splice_junc_to_info_genome' '''
                        if t_splice_junction not in dict_t_splice_junc_to_info_genome:
                            dict_t_splice_junc_to_info_genome[t_splice_junction] = []
                        dict_t_splice_junc_to_info_genome[t_splice_junction].append(
                            (id_tx, strand)
                        )
                        
                        ''' add a record to 'l_df_gtf_splice_donor_and_acceptor_genome' '''
                        l_df_gtf_splice_donor_and_acceptor_genome.extend( [
                            [
                                seqname, # add seqname of the gene
                                left_exon_end, # 1-based SJ donor (last base of the left exon)
                                left_exon_end, # 1-based SJ donor
                                strand, # add strand of the gene
                                id_gene, # add id_gene
                                'L', # add the relative location of exon to intron
                            ],
                            [
                                seqname, # add seqname of the gene
                                right_exon_start, # 1-based SJ acceptor (first base of the right exon)
                                right_exon_start, # 1-based SJ acceptor
                                strand, # add strand of the gene
                                id_gene, # add id_gene
                                'R', # add the relative location of exon to intron
                            ],
                        ] )
                        
                """
                Retrieve splice junction from transcriptomic information
                """
                """ retrieve exons of the transcript in the order """
                l_l_df_gtf_exon_transcriptome_for_the_current_transcript = (
                    []
                )  # initialize a list of list that will contain records for the current transcript
                int_pos_transcript = (
                    0  # initialize the position on the transcript # 0-based coordinates
                )
                for exon_start, exon_end in (
                    arr_df_subset[:, 3:5] if strand == "+" else arr_df_subset[::-1, 3:5]
                ):  # retrieve exons # flip the order of exons if transcript's direction is '-'
                    len_exon = (
                        exon_end - exon_start + 1
                    )  # retrieve the length of the current exon
                    l_l_df_gtf_exon_transcriptome_for_the_current_transcript.append(
                        [
                            id_tx,
                            int_pos_transcript + 1,
                            int_pos_transcript + len_exon,
                            seqname,
                            exon_start,
                            exon_end,
                            strand,
                        ]
                    )  # 1-based coordinates # retain the information of the original exon
                    int_pos_transcript += len_exon
                l_l_df_gtf_exon_transcriptome.extend(
                    l_l_df_gtf_exon_transcriptome_for_the_current_transcript
                )  # append the records to the GTF

                if (
                    int_n_exon > 1
                ):  # skip identification of splice junctions for single-exon transcripts
                    arr_df_gtf_exon_transcriptome_for_the_current_transcript = np.array(
                        l_l_df_gtf_exon_transcriptome_for_the_current_transcript,
                        dtype=object,
                    )  # convert the list of list to array
                    for exon_coord_in_transcriptome, exon_coord_in_genome in zip(
                        arr_df_gtf_exon_transcriptome_for_the_current_transcript[:, 1:3]
                        .ravel()[1:-1]
                        .reshape((int_n_exon - 1, 2)),
                        arr_df_gtf_exon_transcriptome_for_the_current_transcript[
                            :, [4, 5] if strand == "+" else [5, 4]
                        ]
                        .ravel()[1:-1]
                        .reshape((int_n_exon - 1, 2)),
                    ):  # retrieve exon-exon junctions # also use exon coordinates on the genomes # if strand is '-', change the exon start and exon end position
                        left_exon_end, right_exon_start = exon_coord_in_transcriptome
                        left_exon_end_genome, right_exon_start_genome = (
                            exon_coord_in_genome
                            if strand == "+"
                            else exon_coord_in_genome[::-1]
                        )  # reverse the order of splicing junction start and ends (end and start positions of downstream-exon and upstream-exon, respectively)
                        """ retrieve a tuple representing a splice junction and collect information about the splicing junc for each transcript """
                        l_l_df_gtf_splice_junc_transcriptome.append(
                            [
                                id_tx,
                                left_exon_end,
                                left_exon_end,
                                seqname,
                                left_exon_end_genome,
                                right_exon_start_genome,
                                strand,
                            ]
                        )  # left_exon_end, left_exon_end, 1-based coordinates # since the transcript sequence lack intron, the splice junction will be identified by using the coordinate of the end of the previous exon (exon on the left side) # left_exon_end_genome, right_exon_start_genome: 1-based coordinates # collect the corresponding splicing junctiuon on the genome

            """ write a pickle file """
            bk.PICKLE_Write(
                path_file_pickle_dict_t_splice_junc_to_info_genome,
                dict_t_splice_junc_to_info_genome,
            )

            """ compose a dataframe """
            df_intronic_region_near_splice_site = pd.DataFrame(
                l_l_intronic_region_near_splice_site,
                columns=["seqname", "start", "end"],
            )
            df_intronic_region_near_splice_site[
                "feature"
            ] = "intron_near_splice_site"  # name the records
            df_intronic_region_near_splice_site.to_csv(
                f"{path_folder_ref}df_intronic_region_near_splice_site.tsv.gz",
                sep="\t",
                index=False,
            )
            
            """ build mask """
            scidx[
                "dict_seqname_to_mask_gtf_intron_near_splice_site"
            ] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                df_gtf=df_intronic_region_near_splice_site,
                str_feature="intron_near_splice_site",
                remove_chr_from_seqname=True,
                path_folder_output=path_folder_mask_gtf_intron_near_splice_site,
            )

            """ compose a GTF of exons of transcript sequences """ 
            df_gtf_exon_transcriptome = pd.DataFrame(
                l_l_df_gtf_exon_transcriptome,
                columns=[
                    "seqname",
                    "start",
                    "end",
                    "seqname_genome",
                    "start_genome",
                    "end_genome",
                    "strand_genome",
                ],
            )  # compose a dataframe
            for name_col, val in zip(
                ["source", "feature", "score", "strand", "frame"],
                ["processed_by_ouro", "exon", ".", "+", "."],
            ):
                df_gtf_exon_transcriptome[name_col] = val
            df_gtf_exon_transcriptome = df_gtf_exon_transcriptome[
                [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "seqname_genome",
                    "start_genome",
                    "end_genome",
                    "strand_genome",
                ]
            ]  # reorder columns to match that of typical GTF dataframe
            bk.GTF_Write(
                df_gtf_exon_transcriptome,
                path_file_gtf_exon_transcriptome,
                flag_update_attribute=True,
                flag_filetype_is_gff3=False,
            )
            """ load exons of transcript sequences as an interval tree data structure """
            dict_it_exon_transcriptome = bk.GTF_Interval_Tree(
                df_gtf_exon_transcriptome,
                feature=None,
                value=["seqname_genome", "start_genome", "end_genome", "strand_genome"],
            )  # load regulatory annotations # use all features
            bk.PICKLE_Write(
                path_file_pickle_dict_it_exon_transcriptome, dict_it_exon_transcriptome
            )

            """ compose a GTF of splice junction of transcript sequences """
            df_gtf_splice_junc_transcriptome = pd.DataFrame(
                l_l_df_gtf_splice_junc_transcriptome,
                columns=[
                    "seqname",
                    "start",
                    "end",
                    "seqname_genome",
                    "left_exon_end_genome",
                    "right_exon_start_genome",
                    "strand_genome",
                ],
            )  # compose a dataframe
            for name_col, val in zip(
                ["source", "feature", "score", "strand", "frame"],
                ["processed_by_ouro", "splice_junction", ".", "+", "."],
            ):
                df_gtf_splice_junc_transcriptome[name_col] = val
            df_gtf_splice_junc_transcriptome = df_gtf_splice_junc_transcriptome[
                [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "seqname_genome",
                    "left_exon_end_genome",
                    "right_exon_start_genome",
                    "strand_genome",
                ]
            ]  # reorder columns to match that of typical GTF dataframe
            bk.GTF_Write(
                df_gtf_splice_junc_transcriptome,
                path_file_gtf_splice_junc_transcriptome,
                flag_update_attribute=True,
                flag_filetype_is_gff3=False,
            )
            """ load data as an interval tree data structure """
            dict_it_splice_junc_transcriptome = bk.GTF_Interval_Tree(
                df_gtf_splice_junc_transcriptome,
                feature=None,
                value=[
                    "seqname_genome",
                    "left_exon_end_genome",
                    "right_exon_start_genome",
                    "strand_genome",
                ],
            ) # build interval tree
            bk.PICKLE_Write(
                path_file_pickle_dict_it_splice_junc_transcriptome,
                dict_it_splice_junc_transcriptome,
            )
            
            """ compose a GTF of splice donor and acceptor sites of genome """ 
            df_gtf_splice_donor_and_acceptor_genome = pd.DataFrame(
                l_df_gtf_splice_donor_and_acceptor_genome,
                columns=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    'id_gene',
                    'relative_location_of_exon_to_intron'
                ],
            ).drop_duplicates( )  # compose a dataframe # drop duplicate values
            for name_col, val in zip(
                ["source", "feature", "score", "frame"],
                ["processed_by_ouro", "splice_donor_and_acceptor", ".", "."],
            ):
                df_gtf_splice_donor_and_acceptor_genome[name_col] = val # fill required columns
            df_gtf_splice_donor_and_acceptor_genome = df_gtf_splice_donor_and_acceptor_genome[ 
                [
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    'id_gene',
                    'relative_location_of_exon_to_intron',
                ]
            ]  # reorder columns to match that of typical GTF dataframe
            bk.GTF_Write(
                df_gtf_splice_donor_and_acceptor_genome,
                path_file_gtf_splice_donor_and_acceptor_genome,
                flag_update_attribute=True,
                flag_filetype_is_gff3=False,
            ) # write the GTF file
            """ load data as an interval tree data structure """
            dict_it_splice_donor_and_acceptor = bk.GTF_Interval_Tree(
                df_gtf_splice_donor_and_acceptor_genome,
                feature=None,
                value = [ "id_gene", 'relative_location_of_exon_to_intron' ],
            ) # build interval tree
            bk.PICKLE_Write(
                path_file_pickle_dict_it_splice_donor_and_acceptor_genome,
                dict_it_splice_donor_and_acceptor,
            )
            # write a log 
            with open(path_file_flag, "w") as file:
                file.write("completed")
        else:
            scidx[
                "dict_seqname_to_mask_gtf_intron_near_splice_site"
            ] = bk.GTF_Build_Mask(
                dict_seqname_to_len_seq=dict_seqname_to_len_seq,
                path_folder_output=path_folder_mask_gtf_intron_near_splice_site,
            )
        logger.info(
            "[Completed] the construction of the masks of introns near splice junctions was completed."
        )

        """ map id_tx to id_gene and vice versa """
        """ map id_tx to name_tx """
        (
            scidx["dict_id_gene_to_l_id_tx"],
            scidx["dict_id_tx_to_id_gene"],
            scidx["dict_id_tx_to_name_tx"],
        ) = (dict(), dict(), dict())
        for id_gene, id_tx, name_tx in df_gtf_transcript[
            [
                str_name_gtf_attr_for_id_gene,
                str_name_gtf_attr_for_id_transcript,
                str_name_gtf_attr_for_name_transcript,
            ]
        ].values:
            scidx["dict_id_tx_to_id_gene"][id_tx] = id_gene
            scidx["dict_id_tx_to_name_tx"][id_tx] = (
                name_tx if not isinstance(name_tx, float) else id_tx
            )  # if 'name_tx' is available, use 'name_tx'. if 'name_tx' is not available, use id_tx as 'name_tx'
            if id_gene not in scidx["dict_id_gene_to_l_id_tx"]:
                scidx["dict_id_gene_to_l_id_tx"][id_gene] = []
            scidx["dict_id_gene_to_l_id_tx"][id_gene].append(id_tx)

        """ prepare efficient access of 'df_gtf_gene' using 'id_gene' """
        scidx["dict_index_df_gtf_gene"] = bk.DF_Build_Index_Using_Dictionary(
            df_gtf_gene, l_col_for_index=str_name_gtf_attr_for_id_gene
        )
        scidx["arr_data_df_gtf_gene"] = df_gtf_gene[
            [
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
                str_name_gtf_attr_for_name_gene,
            ]
        ].values
        logger.info(
            "[Completed] completed the loading of miscellaneous data structures were completed."
        )

        joblib.dump(
            scidx, path_file_scidx
        )  # dump the sparse matrix for paralleled access

    # initialize manager proxy objects to avoid memory bloating
    import asyncio

    lp = asyncio.get_event_loop()

    async def __load_managers_of_a_data_object(
        name_data: str,
        type_managed_data: Literal["HostedDictIntervalTree", "HostedDict"],
        path_file_pickle: str,
        int_num_manager_processes: int,
    ):
        """# 2023-01-09 01:03:29
        initialize managed data

        int_num_manager_processes : int = 0 # the number of managed processes to start.
        """

        async def __load_a_manager_of_a_data_object(
            name_data: str,
            type_managed_data: Literal["HostedDictIntervalTree", "HostedDict"],
            path_file_pickle: str,
        ):
            """# 2023-01-09 16:23:09
            load and start a manager, and return a manager and the proxy object from the manager
            """
            logger.info(f"loading a manager for '{name_data}' data object started.")
            # initialize the manager
            manager = bk.ManagerReadOnly()
            manager.start()  # start the manager
            managed_data = getattr(manager, type_managed_data)(path_file_pickle)
            logger.info(f"loading a manager for '{name_data}' data object completed.")
            return {
                "name_data": name_data,
                "type_managed_data": type_managed_data,
                "manager": manager,
                "managed_data": managed_data,
            }

        if int_num_manager_processes > 0:
            return list(
                __load_a_manager_of_a_data_object(
                    name_data, type_managed_data, path_file_pickle
                )
                for i in range(int_num_manager_processes)
            )  # load data into manager processes

    async def __load_managers():
        """# 2023-01-09 13:57:09
        load managers of data objects asynchronously
        """
        l_cor = []  # initialize the list of coroutin objects

        # for each data object, either gather futures of manager objects or directly load data in the main process
        for name_data, type_managed_data, path_file_pickle in zip(
            [
                "dict_it_exon_transcriptome",
                "dict_it_rpmk",
                "dict_it_splice_junc_transcriptome",
                "dict_it_splice_donor_and_acceptor_genome",
                "dict_it_exon",
                "dict_it_reg",
                "dict_it_promoter",
                "dict_fa_transcriptome",
                "dict_t_splice_junc_to_info_genome",
            ],
            [
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDictIntervalTree",
                "HostedDict",
                "HostedDict",
            ],
            [
                f"{path_folder_ref}dict_it_exon_transcriptome.pickle",
                f"{path_folder_ref}dict_it_rpmk.pickle",
                f"{path_folder_ref}dict_it_splice_junc_transcriptome.pickle",
                f"{path_folder_ref}dict_it_splice_donor_and_acceptor_genome.pickle",
                f"{path_folder_ref}dict_it_exon.pickle",
                f"{path_folder_ref}dict_it_reg.pickle",
                f"{path_folder_ref}dict_it_promoter.pickle",
                f"{path_folder_ref}dict_fa_transcriptome.pickle",
                f"{path_folder_ref}dict_t_splice_junc_to_info_genome.pickle",
            ],
        ):
            int_num_manager_processes = dict_num_manager_processes_for_each_data_object[
                name_data
            ]  # retrieve 'int_num_manager_processes'
            if int_num_manager_processes == 0:
                # if 'int_num_manager_processes' == 0, directly load data into the current process
                scidx[name_data] = bk.PICKLE_Read(path_file_pickle)
            else:
                # load data into manager processes (gather future object of the manager processes)
                l_cor.extend(
                    lp.run_until_complete(
                        __load_managers_of_a_data_object(
                            name_data,
                            type_managed_data,
                            path_file_pickle,
                            int_num_manager_processes,
                        )
                    )
                )

        # retrieve managed data objects and manager objects and save these object in scidx data object
        if len(l_cor) > 0:
            for e in await asyncio.gather(*l_cor):
                name_data = e["name_data"]
                # initialize
                if f"l_managed_{name_data}" not in scidx:
                    scidx[f"l_managed_{name_data}"] = []
                    scidx[f"l_manager_of_{name_data}"] = []
                scidx[f"l_managed_{name_data}"].append(e["managed_data"])
                scidx[f"l_manager_of_{name_data}"].append(e["manager"])

    lp.run_until_complete(__load_managers())  # load managers

    return scidx

def Convert_df_count_to_MTX_10X(
    path_file_df_count: str,
    path_folder_mtx_10x_output: str,
    path_folder_mtx_10x_filtered_output: str,
    chunksize: int = 1000000,
    int_min_count_features_for_filtering_barcodes: int = 50,
    flag_output_dtype_is_integer : bool = True,
):
    """# 2023-01-06 23:46:20
    convert df_count (ouro output) to 10X MTX (matrix market) format in a memory-efficient manner.

    path_file_df_count : str, # file path to 'df_count'
    path_folder_mtx_10x_output : str, # a folder containing 10x output matrix (unfiltered)
    path_folder_mtx_10x_filtered_output : str, # a folder containing 10x output matrix (filtered)
    chunksize : int = 500000,
    int_min_count_features_for_filtering_barcodes : int = 50, # the minimum number of features in a barcode to be included in the filtered output
    flag_output_dtype_is_integer : bool = True, # a boolean flag indicating the output dtype is integer dtype. Set this flag to False if the output dtype is float
    """
    # create output folders
    os.makedirs(path_folder_mtx_10x_output, exist_ok=True)
    os.makedirs(path_folder_mtx_10x_filtered_output, exist_ok=True)

    # retrieve unique feature/barcode information from df_count without loading entire data in the memory
    bk.DF_Deduplicate_without_loading_in_memory(
        path_file_df_count,
        f"{path_folder_mtx_10x_output}_features.tsv.gz",
        l_col_for_identifying_duplicates=["feature", "id_feature"],
        str_delimiter="\t",
    )
    res = bk.DF_Deduplicate_without_loading_in_memory(
        path_file_df_count,
        f"{path_folder_mtx_10x_output}_barcodes.tsv.gz",
        l_col_for_identifying_duplicates="barcode",
        str_delimiter="\t",
    )  # collect the number of records
    int_num_lines = res["int_num_lines"]
    s_num_records_for_each_barcode = pd.Series(
        res["dict_t_val_count"]
    )  # retrieve the number of records for each barcode
    del res
    s_num_records_for_each_barcode = s_num_records_for_each_barcode[
        s_num_records_for_each_barcode >= int_min_count_features_for_filtering_barcodes
    ]  # filter barcodes using the given setting
    df_barcode_filtered = pd.DataFrame(
        {"barcode": s_num_records_for_each_barcode.index.values}
    )  # compose a dataframe containing filtered barcodes

    # read features and barcode information
    df_barcode = pd.read_csv(
        f"{path_folder_mtx_10x_output}_barcodes.tsv.gz", sep="\t", usecols=["barcode"]
    )
    df_feature = pd.read_csv(
        f"{path_folder_mtx_10x_output}_features.tsv.gz",
        sep="\t",
        usecols=["feature", "id_feature"],
    )
    df_feature = df_feature.loc[:, ["id_feature", "feature"]]
    df_feature["10X_type"] = "Gene Expression"
    # save feature/cell metadata
    df_barcode.to_csv(
        f"{path_folder_mtx_10x_output}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_barcode_filtered.to_csv(
        f"{path_folder_mtx_10x_filtered_output}barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_feature.to_csv(
        f"{path_folder_mtx_10x_output}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )
    df_feature.to_csv(
        f"{path_folder_mtx_10x_filtered_output}features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
    )

    # retrieve barcode/feature to integer representation of barcode/feature mapping
    dict_to_int_barcode = dict(
        (e, i + 1) for i, e in enumerate(df_barcode.iloc[:, 0].values)
    )
    dict_to_int_barcode_filtered = dict(
        (e, i + 1) for i, e in enumerate(df_barcode_filtered.iloc[:, 0].values)
    )
    dict_to_int_feature = dict(
        (e, i + 1) for i, e in enumerate(df_feature.iloc[:, 0].values)
    )

    (
        int_num_features,
        int_num_barcodes,
        int_num_barcodes_filtered,
        int_num_records,
        int_num_records_filtered,
    ) = (
        len(df_feature),
        len(df_barcode),
        len(df_barcode_filtered),
        int_num_lines,
        s_num_records_for_each_barcode.sum(),
    )  # retrieve metadata of the output matrix
    del (
        df_feature,
        df_barcode,
        df_barcode_filtered,
        s_num_records_for_each_barcode,
    )  # delete objects

    # write mtx file
    str_dtype_mtx_header = 'integer' if flag_output_dtype_is_integer else 'real' # retrieve dtype name for the matrix header
    with gzip.open(f"{path_folder_mtx_10x_output}matrix.mtx.gz", "wb") as newfile:
        newfile.write(
            f"""%%MatrixMarket matrix coordinate {str_dtype_mtx_header} general\n%\n{int_num_features} {int_num_barcodes} {int_num_records}\n""".encode()
        )
        with gzip.open(
            f"{path_folder_mtx_10x_filtered_output}matrix.mtx.gz", "wb"
        ) as newfile_filtered:
            newfile_filtered.write(
                f"""%%MatrixMarket matrix coordinate {str_dtype_mtx_header} general\n%\n{int_num_features} {int_num_barcodes_filtered} {int_num_records_filtered}\n""".encode()
            )

            # iterate through each chunk
            for df_chunk in pd.read_csv(
                path_file_df_count,
                iterator=True,
                header=0,
                chunksize=chunksize,
                sep="\t",
                usecols=["id_feature", "barcode", "read_count"],
            ):
                df_chunk = df_chunk[
                    ["id_feature", "barcode", "read_count"]
                ]  # reorder columns
                mask_filtered = np.array(
                    list(
                        e in dict_to_int_barcode_filtered
                        for e in df_chunk.barcode.values
                    ),
                    dtype=bool,
                )  # retrieve a mask for filtering records
                df_chunk["id_feature"] = df_chunk.id_feature.apply(
                    MAP.Map(dict_to_int_feature).a2b
                )
                df_chunk_filtered = df_chunk[mask_filtered]  # filter records
                df_chunk["barcode"] = df_chunk.barcode.apply(
                    MAP.Map(dict_to_int_barcode).a2b
                )
                df_chunk_filtered["barcode"] = df_chunk_filtered.barcode.apply(
                    MAP.Map(dict_to_int_barcode_filtered).a2b
                )
                df_chunk.to_csv(newfile, sep=" ", header=None, index=False)
                df_chunk_filtered.to_csv(
                    newfile_filtered, sep=" ", header=None, index=False
                )  # write filtered records
                del mask_filtered, df_chunk, df_chunk_filtered

    # delete temporary files
    os.remove(f"{path_folder_mtx_10x_output}_barcodes.tsv.gz")
    os.remove(f"{path_folder_mtx_10x_output}_features.tsv.gz")

def _update_size_distribution( new_size : int, arr_dist = None, int_default_size : int = 10000 ) :
    """ # 2023-08-05 09:05:32 
    update size distribution using the new size
    """
    if arr_dist is None :
        arr_dist = np.zeros( max( int( new_size * 2 ), int_default_size ), dtype = int ) # initialize an array containing the size distribution # the array will accomodate sizes upto 2x of the initial input
    if len( arr_dist ) <= new_size : # when the new molecule size exceed that of the container
        arr_dist_new = np.zeros( int( new_size * 2 ), dtype = int ) # initialize an array containing the size distribution # the array will accomodate sizes upto 2x of the initial input
        arr_dist_new[ : len( arr_dist ) ] = arr_dist[ : ] # copy 'arr_dist' to 'arr_dist_new'
        arr_dist = arr_dist_new # replace existing 'arr_dist'
    arr_dist[ new_size ] += 1 # update the size distribution
    return arr_dist

def _batch_update_size_distribution( l_new_size, arr_dist = None, int_default_size : int = 10000 ) :
    """ # 2023-08-05 09:05:32 
    batch update size distribution using the list of new size (a list of integers)
    """
    n = len( l_new_size ) # retrieve the number of new sizes
    if arr_dist is None :
        arr_dist = np.zeros( max( int( np.max( l_new_size ) if n > 30 else max( l_new_size ) ) * 2, int_default_size ), dtype = int ) # initialize an array containing the size distribution # the array will accomodate sizes upto 2x of the initial input
    for new_size in l_new_size : # for each new size 
        if len( arr_dist ) <= new_size : # when the new molecule size exceed that of the container
            arr_dist_new = np.zeros( int( new_size * 2 ), dtype = int ) # initialize an array containing the size distribution # the array will accomodate sizes upto 2x of the initial input
            arr_dist_new[ : len( arr_dist ) ] = arr_dist[ : ] # copy 'arr_dist' to 'arr_dist_new'
            arr_dist = arr_dist_new # replace existing 'arr_dist'
        arr_dist[ new_size ] += 1 # update the size distribution
    return arr_dist

def _combine_size_distribution( arr_dist_1, arr_dist_2 ) :
    """ # 2023-08-03 11:41:13 
    combine two size distributions. one of the distributions will be modified in-place
    """
    ''' handle cases when one if the distribution is empty '''
    if arr_dist_1 is None :
        return arr_dist_2
    if arr_dist_2 is None :
        return arr_dist_1

    ''' when both of the distributions are not empty '''
    len_arr_dist_1, len_arr_dist_2 = len( arr_dist_1 ), len( arr_dist_2 )
    if len_arr_dist_1 >= len_arr_dist_2 : # if the length of distribution 1 is larger than the distribution 2
        arr_dist_1[ : len_arr_dist_2 ] += arr_dist_2 # add distribution 2 to distribution 1
        return arr_dist_1
    else :
        arr_dist_2[ : len_arr_dist_1 ] += arr_dist_1
        return arr_dist_2

def _index_array( l_index : list ) :
    """ # 2023-08-12 20:46:39 
    return a dictionary where key = unique value of 'l_index' and value = list of integer indices of the entries that are equal to the unique value.
    Of note, ignore 'float' type values, including np.nan values.
    """
    dict_index = dict( )
    for i, index in enumerate( l_index ) :
        if isinstance( index, float ) :
            continue
        if index not in dict_index :
            dict_index[ index ] = [ ]
        dict_index[ index ].append( i )
    return dict_index

def _combine_dictionary_of_size_distributions( dict_arr_dist_existing : dict, dict_arr_dist_new : dict ) :
    """ # 2023-08-14 18:27:36 
    combine two 'dict_arr_dist' object
    """
    for e in set( dict_arr_dist_existing ).intersection( dict_arr_dist_new ) : # for existing keys in the container
        dict_arr_dist_existing[ e ] = _combine_size_distribution( dict_arr_dist_existing[ e ], dict_arr_dist_new[ e ] ) # update the distribution
    for e in set( dict_arr_dist_new ).difference( dict_arr_dist_existing ) : # new keys that do not exist in the container
        dict_arr_dist_existing[ e ] = dict_arr_dist_new[ e ] # copy distribution of a new key to existing container
    return dict_arr_dist_existing

def _batch_update_dictionary_of_size_distributions( dict_arr_dist : dict, dict_l_len : dict ) :
    """ # 2023-08-14 18:27:36 
    dict_arr_dist : dict 
    dict_l_len : dict

    update size distributions in 'dict_arr_dist' using the list of sizes in 'dict_l_len'
    """
    for e in set( dict_arr_dist ).intersection( dict_l_len ) : # for existing keys in the container
        dict_arr_dist[ e ] = _batch_update_size_distribution( l_new_size = dict_l_len[ e ], arr_dist = dict_arr_dist[ e ] ) # update the distribution
    for e in set( dict_l_len ).difference( dict_arr_dist ) : # new keys that do not exist in the container
        dict_arr_dist[ e ] = _batch_update_size_distribution( l_new_size = dict_l_len[ e ], arr_dist = None ) # create a new distribution
    return dict_arr_dist
    
def _get_df_bar( dict_arr_dist : dict, l_name_type_dist : Union[ list, None ] = None, int_size_bin_in_base_pairs : int = 50, int_max_size_in_base_pairs : int = 7500 ) :
    """ # 2023-08-15 15:42:39 
    compose a dataframe containing summarized distributions
    """
    # set default values
    if l_name_type_dist is None : # if 'l_name_type_dist' is None, use all key values of the given 'dict_arr_dist'
        l_name_type_dist = list( dict_arr_dist )
    
    int_max_length = math.ceil( max( len( dict_arr_dist[ e ] ) for e in l_name_type_dist if dict_arr_dist[ e ] is not None ) / int_size_bin_in_base_pairs ) * int_size_bin_in_base_pairs # retrieve max length size of the distributions (should be multiple of 'int_size_bin_in_base_pairs')
    l_arr = [ ]
    for e in l_name_type_dist :
        ''' compose array of distribution '''
        arr = np.zeros( int_max_length, dtype = int )
        arr_dist = dict_arr_dist[ e ]
        if arr_dist is not None : # if 'arr_dist' is not empty copy the 'arr_dist' to the array
            arr[ : len( arr_dist ) ] = arr_dist

        ''' summarize and clip the distribution '''
        arr = arr.reshape( ( int( int_max_length / int_size_bin_in_base_pairs ), int_size_bin_in_base_pairs ) ).sum( axis = 1 ) # summarize distribution
        arr = arr[ : math.ceil( int_max_size_in_base_pairs / int_size_bin_in_base_pairs ) ] # clip distribution

        l_arr.append( arr )
    df_bar = pd.DataFrame( np.vstack( l_arr ), index = l_name_type_dist, columns = np.arange( 1, len( l_arr[ 0 ] ) + 1, dtype = int ) * int_size_bin_in_base_pairs )
    return df_bar

def _draw_bar_plot( 
    df_bar, 
    l_status : list, 
    title : str = '', 
    y_format : str = ':.3f', 
    flag_save_figure : bool = False, 
    path_folder_graph : Union[ str, None ] = None,
    l_color : Union[ List[ str ], None ] = None,
    flag_use_proportion : bool = False,
    float_ratio_padding_y_axis : float = 0.05,
    xaxis_range : Union[ List[ float ], None ] = None,
    int_min_total_counts_bin_for_proportion_calculation : int = 10,
) :
    ''' # 2023-10-25 21:16:00 
    draw barplot for the given 'df_bar'
    Return plotly fig
    maximum number of l_status (unique color will be assigned to each) is 56
    l_color : Union[ List[ str ], None ] = None, # list of color for each status
    flag_use_proportion : bool = False, # if True, draw proportion graph
    float_ratio_padding_y_axis : float = 0.05, # padding added to the max value of the y-axis
    xaxis_range : Union[ List[ float ], None ] = None, # if given, set the xaxis range using the given start and end positions
    int_min_total_counts_bin_for_proportion_calculation : int = 10, # min total counts for a bin to plot the proportions of categories for the bin. bins below this total counts will be shown as '0'
    '''
    import plotly.express as px
    import plotly.graph_objects as go
    x = df_bar.columns.values # retrieve x categories
    if l_color is None :
        l_color = px.colors.qualitative.Dark24 + px.colors.qualitative.Pastel2 + px.colors.qualitative.Light24  # retrieve unique colors for each status # plotly express can draw barplot more easily, but for annotating each status with the same color annotation, for loop with go.Bar will be used instead
    set_valid_status = set( df_bar.index.values ) # retrieve a set of valid status
    df_bar = df_bar.copy( ) # copy the data
    df_bar = df_bar.loc[ list( e for e in l_status if e in set_valid_status ) ] # drop invalid categories
    if flag_use_proportion : # if use proportions, calculate proportions
        s_sum = df_bar.sum( axis = 0 )
        df_bar = df_bar / s_sum
        df_bar.loc[ :, s_sum.values < int_min_total_counts_bin_for_proportion_calculation ] = 0 # does not plot bins below 'int_min_total_counts_bin_for_proportion_calculation'
        df_bar.fillna( 0, inplace = True )
    
    l_go = list( go.Bar( x = x, y = df_bar.loc[ str_status ].values, name = str_status, marker_color = str_color, hovertemplate = 'size_bin: <b>%{x}</b><br><br>proportion: <b>%{y' + y_format + '}</b>', width = df_bar.columns[ 1 ] - df_bar.columns[ 0 ] ) for str_status, str_color in zip( l_status, l_color ) if str_status in set_valid_status ) # retrieve graph object for each valid str_status # infer width from the input 'df_bar'

    # compose a bar plot
    fig = go.Figure( l_go[ 0 ] )
    for go_bar in l_go[ 1 : ] :
        fig.add_trace( go_bar )
    fig.update_traces(marker=dict( line = dict(width=0) )) 
    fig.update_layout( barmode = 'stack', xaxis = { 'categoryorder' : 'category ascending' }, title_text = title, plot_bgcolor='white' )
    fig.update_layout( yaxis_range = [ 0, 1 + float_ratio_padding_y_axis ] if flag_use_proportion else [ 0, df_bar.sum( axis = 0 ).max( ) * ( 1 + float_ratio_padding_y_axis ) ] ) # update y-axis range        
    if xaxis_range is not None : # update xaxis_range
        fig.update_layout( xaxis_range = [ xaxis_range[ 0 ], xaxis_range[ 1 ] ] ) 
    fig.update_xaxes( mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='#f8f8f8' )
    fig.update_yaxes( mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='#f8f8f8' )
    if flag_save_figure and path_folder_graph is not None :
        fig.write_html( f"{path_folder_graph}distribution.bar.{'proportion' if flag_use_proportion else 'read_count'}.{title}.html" ) # write HTML file
    else :
        return fig
    
def _check_binary_flags( flags : int, int_bit_flag_position : int ) :
    """ # 2023-08-08 22:47:02 
    check a flag in the binary flags at the given position
    """
    return ( flags & ( 1 << int_bit_flag_position ) ) > 0 

def _argmin(seq):
    """ # 2023-09-20 19:47:25 
    reference: https://gist.github.com/mscheltienne/4b02ff0095e2d847e1e572d0aae216c6
    Get the index of the minimum item in a list or dict.
    
    Parameters
    ----------
    seq : list | dict
        The list or dict on which to find its min.
    
    Returns:
    --------
    int
        The minimum value of the list or dict.
    """
    if isinstance(seq, (list, tuple)):
        return min(range(len(seq)), key=seq.__getitem__)

    elif isinstance(seq, dict):
        return min(seq, key=seq.__getitem__)

    else:
        raise TypeError

def _argmax(seq):
    """ # 2023-09-20 19:47:19 
    reference: https://gist.github.com/mscheltienne/4b02ff0095e2d847e1e572d0aae216c6
    Get the index of the maximum item in a list or dict
    
    Parameters
    ----------
    seq : list | dict
        The list or dict on which to find its max.
    
    Returns:
    --------
    int
        The maximum value of the list or dict.
    """
    if isinstance(seq, (list, tuple)):
        return max(range(len(seq)), key=seq.__getitem__)

    elif isinstance(seq, dict):
        return max(seq, key=seq.__getitem__)

    else:
        raise TypeError
    
def LongFilterNSplit(
    flag_usage_from_command_line_interface: bool = False,
    path_file_minimap_index_genome: Union[str, None] = None,
    l_path_file_minimap_index_unwanted: List[str] = [ ],
    l_path_file_fastq_input: Union[list, None] = None,
    l_path_folder_output: [list[str], None] = None,
    n_threads: int = 32,
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_base_pairs_in_a_batch : int = 2_500_000, # the number of base pairs in a batch
    float_memory_in_GiB: float = 50,
    verbose: bool = True,
    str_minimap_aligner_preset : str = 'splice', # preset of the minimap2 aligner
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome
    int_size_window_for_searching_poly_a_tail : int = 15, # the size of the window from the end of the alignment to search for poly A tail.
    int_max_size_intervening_sequence_between_alignment_end_and_poly_A : int = 20, # max size of the intervening sequence between alignment end position and poly A tract. it will be applied for both internal poly A or external (enzymatically attached) poly A.
    float_min_A_frequency_for_identifying_poly_A : float = 0.75, # the minimum frequency to determine a sequence contains a poly A tract
    int_min_size_intervening_sequence_for_splitting : int = 150, # the minimum length of intervening sequence between alignments for splitting the reads
    int_max_intron_size_for_determining_chimeric_molecule : int = 200_000, # the maximum allowed intron size for classifying considering the molecule as a intra-chromosomal chimeric read
    int_max_read_length : int = 20_000, # the maximum read length to analyze. If the speed of the analysis seems to be slower than expected, try lowering this parameter to filter out repeat-containing artifact reads present in long-read sequencing data, which takes a long time to align and filter out based on the alignment profile.
    flag_recover_original_molecule_before_self_circularization_and_digestion : bool = False, # by default, this flag is set to False. However, when Ourn-Enrich was applied and the sequenced molecules are rearranged so that the adaptors are located in the middle, please set this flag to True, which enables a recovery mode to reconstruct the original molecules from the rearranged sequences.
    int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion : int = 30, # the max coordinate difference at the predicted cut site for recovering original molecule before self-circularization. this argument is only active when 'flag_recover_original_molecule_before_self_circularization' is True
    am_genome = None, # mappy aligner for genome (optional. if given, will override 'path_file_minimap_index_genome' argument)
    l_am_unwanted : Union[ None, List ] = None, # mappy aligner for unwanted sequences (optional. if given, will override 'l_path_file_minimap_index_unwanted' argument)
) -> None :
    """# 2023-11-02 20:16:31 
    
    flag_usage_from_command_line_interface: bool = False,
    path_file_minimap_index_genome: Union[str, None] = None, # required for identifying valid regions of a read and identify chimeric transcripts
    l_path_file_minimap_index_unwanted: List[str] = [ ], # a list of minimap indices of sequences to which unwanted reads can originate (e.g., mt-dna, rRNA repeat, etc.) in a decreasing order of priority
    l_path_file_fastq_input: Union[list, None] = None,
    l_path_folder_output: [list[str], None] = None,
    n_threads: int = 32,
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_base_pairs_in_a_batch : int = 2_500_000, # the number of base pairs in a batch
    float_memory_in_GiB: float = 50,
    str_minimap_aligner_preset = 'splice', # preset of the minimap2 aligner
    verbose: bool = True,
    int_min_mapq = 1, # minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome
    int_size_window_for_searching_poly_a_tail : int = 16, # the size of the window from the end of the alignment to search for poly A tail.
    int_max_size_intervening_sequence_between_alignment_end_and_poly_A : int = 20, # max size of the intervening sequence between alignment end position and poly A tract. it will be applied for both internal poly A or external (enzymatically attached) poly A.
    float_min_A_frequency_for_identifying_poly_A : float = 0.75, # the minimum frequency to determine a sequence contains a poly A tract
    int_max_intron_size_for_determining_chimeric_molecule : int = 200000, # the maximum allowed intron size for classifying considering the molecule as a intra-chromosomal chimeric read
    int_max_read_length : int = 30_000, # the maximum read length to analyze. If the speed of the analysis seems to be slower than expected, try lowering this parameter to filter out repeat-containing artifact reads present in long-read sequencing data, which takes a long time to align and filter out based on the alignment profile.
    flag_recover_original_molecule_before_self_circularization_and_digestion : bool = False, # by default, this flag is set to False. However, when Ourn-Enrich was applied and the sequenced molecules are rearranged so that the adaptors are located in the middle, please set this flag to True, which enables a recovery mode to reconstruct the original molecules from the rearranged sequences.
    int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion : int = 30, # the max coordinate difference at the predicted cut site for recovering original molecule before self-circularization. this argument is only active when 'flag_recover_original_molecule_before_self_circularization_and_digestion' is True
    am_genome = None, # mappy aligner for genome (optional. if given, will override 'path_file_minimap_index_genome' argument)
    int_min_size_intervening_sequence_for_splitting : int = 150, # the minimum length of intervening sequence between alignments for splitting the reads
    l_am_unwanted : Union[ None, List ] = None, # mappy aligner for unwanted sequences (optional. if given, will override 'l_path_file_minimap_index_unwanted' argument)

    returns
    
    * of note, strand information is not used for identifying chimeric molecules, since hairpin formation during reverse transcription can lead to molecules with hairpin alignment patterns, which are not chimeric molecules. (reference: 10x Genomics' technical document)
    """
    """
    Parse arguments
    """
    if flag_usage_from_command_line_interface:  # parse arguments
        """parse arguments when the function was called from the command-line interface"""
        # {  } # unused arguments
        # command line arguments
        parser = argparse.ArgumentParser(
            description=str_description,
            usage="ourotools LongFilterNSplit",
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument("LongFilterNSplit")

        arg_grp_general = parser.add_argument_group("General")
        arg_grp_general.add_argument(
            "-q",
            "--l_path_file_fastq_input",
            help="",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-o",
            "--l_path_folder_output",
            help="",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-t",
            "--n_threads",
            help="(default: 32) the number of processors to use concurrently.",
            default=32,
            type=int,
        )
        arg_grp_general.add_argument(
            "-b",
            "--int_num_base_pairs_in_a_batch",
            help="(default: 2,500,000) the number of base pairs in a batch.",
            default = 2_500_000,
            type=int,
        )
        arg_grp_general.add_argument(
            "-s",
            "--int_num_samples_analyzed_concurrently",
            help="(default: 2) the number of samples that can be analyzed concurrently.",
            default=2,
            type=int,
        )
        arg_grp_general.add_argument(
            "-m",
            "--float_memory_in_GiB",
            help="(default: 50) the maximum memory usage of the pipeline in GiB",
            default=50,
            type=float,
        )
        arg_grp_general.add_argument(
            "-v", 
            "--verbose", 
            help="turn on verbose mode", 
            action="store_true"
        )
        arg_grp_alignment = parser.add_argument_group("Alignment")
        arg_grp_alignment.add_argument(
            "-i",
            "--path_file_minimap_index_genome",
            help="",
            type=str,
        )
        arg_grp_alignment.add_argument(
            "-u",
            "--l_path_file_minimap_index_unwanted",
            help="",
            nargs="*",
        )
        arg_grp_alignment.add_argument(
            "-Q", 
            "--int_min_mapq", 
            help="(default: 1) minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome", 
            default=1,
            type=int,
        )
        arg_grp_alignment.add_argument(
            "-x",
            "--str_minimap_aligner_preset",
            help="(default: 'splice') preset of the minimap2 aligner",
            default="splice",
            type=str,
        )
        arg_grp_alignment.add_argument(
            "-L", 
            "--int_max_read_length", 
            help="(default: 30,000) the maximum read length to analyze. If the speed of the analysis seems to be slower than expected, try lowering this parameter to filter out repeat-containing artifact reads present in long-read sequencing data, which takes a long time to align and filter out based on the alignment profile.", 
            default=30_000,
            type=int,
        )
        arg_grp_poly_a_tail_detection = parser.add_argument_group("Poly A tail detection")
        arg_grp_poly_a_tail_detection.add_argument(
            "-w",
            "--int_size_window_for_searching_poly_a_tail",
            help="(default: 16) the size of the window from the end of the alignment to search for poly A tail.",
            default=16,
            type=int,
        )
        arg_grp_poly_a_tail_detection.add_argument(
            "-A",
            "--float_min_A_frequency_for_identifying_poly_A",
            help="(default: 0.75) the minimum frequency to determine a sequence contains a poly A tract.",
            default=0.75,
            type=float,
        )
        arg_grp_poly_a_tail_detection.add_argument(
            "-I",
            "--int_max_size_intervening_sequence_between_alignment_end_and_poly_A",
            help="(default: 20) the maximum size of the intervening sequence between alignment end position and poly A tract. it will be applied for both internal poly A or external (enzymatically attached) poly A.",
            default=20,
            type=int,
        )
        
        arg_grp_read_splitting = parser.add_argument_group("Read Splitting ")
        arg_grp_read_splitting.add_argument(
            "-S",
            "--int_min_size_intervening_sequence_for_splitting",
            help="(default: 150) the minimum length of intervening sequence between alignments for splitting the reads.",
            default=150,
            type=int,
        )
        arg_grp_read_splitting.add_argument(
            "-C",
            "--int_max_intron_size_for_determining_chimeric_molecule",
            help="(default: 200,000) the maximum allowed intron size for classifying considering the molecule as a intra-chromosomal chimeric read.",
            default=200000,
            type=int,
        )
        arg_grp_read_splitting.add_argument(
            "-r", 
            "--flag_recover_original_molecule_before_self_circularization_and_digestion", 
            help="by default, this flag is set to False. However, when Ourn-Enrich was applied and the sequenced molecules are rearranged so that the adaptors are located in the middle, please set this flag to True, which enables a recovery mode to reconstruct the original molecules from the rearranged sequences.", 
            action="store_true"
        )
        arg_grp_read_splitting.add_argument(
            "-D",
            "--int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion",
            help="(default: 30) the max coordinate difference at the predicted cut site for recovering original molecule before self-circularization. this argument is only active when 'flag_recover_original_molecule_before_self_circularization_and_digestion' is True",
            default=30,
            type=int,
        )
        
        args = parser.parse_args()

        path_file_minimap_index_genome = args.path_file_minimap_index_genome
        l_path_file_fastq_input = args.l_path_file_fastq_input
        l_path_folder_output = args.l_path_folder_output
        n_threads = args.n_threads
        int_num_samples_analyzed_concurrently = args.int_num_samples_analyzed_concurrently
        float_memory_in_GiB = args.float_memory_in_GiB
        verbose = args.verbose
        int_num_base_pairs_in_a_batch = args.int_num_base_pairs_in_a_batch
        l_path_file_minimap_index_unwanted = args.l_path_file_minimap_index_unwanted
        int_min_mapq = args.int_min_mapq
        str_minimap_aligner_preset = args.str_minimap_aligner_preset
        int_size_window_for_searching_poly_a_tail = args.int_size_window_for_searching_poly_a_tail
        int_max_size_intervening_sequence_between_alignment_end_and_poly_A = args.int_max_size_intervening_sequence_between_alignment_end_and_poly_A
        float_min_A_frequency_for_identifying_poly_A = args.float_min_A_frequency_for_identifying_poly_A
        int_min_size_intervening_sequence_for_splitting = args.int_min_size_intervening_sequence_for_splitting
        int_max_intron_size_for_determining_chimeric_molecule = args.int_max_intron_size_for_determining_chimeric_molecule
        int_max_read_length = args.int_max_read_length
        flag_recover_original_molecule_before_self_circularization_and_digestion = args.flag_recover_original_molecule_before_self_circularization_and_digestion
        int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion = args.int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion

    """
    Start of the pipeline
    """
    logger.info(str_description)
    logger.info(
        "Ouro-Tools LongReadFilterNSplit, a preprocessing pipeline for filtering undesired reads and spliting chimeric reads FASTQ files"
    )
    logger.info(f"Started.")

    """ handle special cases and invalid inputs """
    if l_path_file_fastq_input is None or ( path_file_minimap_index_genome is None and am_genome is None ) : # when both the minimap2 aligner and index path are not given
        logger.error(
            "Required argument(s) is missing. to view help message, type -h or --help"
        )
        return -1

    """ process required input directories """
    path_file_minimap_index_genome = os.path.abspath(path_file_minimap_index_genome)
    l_path_file_minimap_index_unwanted = list( os.path.abspath( e ) for e in l_path_file_minimap_index_unwanted )

    """ process input directory  """
    l_path_file_fastq_input = list(
        os.path.abspath(path_file_fastq_input)
        for path_file_fastq_input in l_path_file_fastq_input
    )
    if l_path_folder_output is not None:
        """# when a valid list of output folders were given # ensure directories of the output folder ends with '/' characters"""
        l_path_folder_output = list(
            os.path.abspath(path_folder) + "/" for path_folder in l_path_folder_output
        )
    else:
        """# compose a list of default 'path_folder_output' values for the given list of input files"""
        l_path_file_fastq_input_reversed = deepcopy(
            l_path_file_fastq_input[::-1]
        )  # reverse the input file paths so that pop operation yield the element located at the front
        l_path_folder_output = []
        for str_mode_ouro_count in l_str_mode_ouro_count:
            path_file_fastq = l_path_file_fastq_input_reversed.pop()
            path_folder_output = (
                f"{path_file_fastq.rsplit( '/', 1 )[ 0 ]}LongFilterNSplit_output/"
            )
            l_path_folder_output.append(path_folder_output)

    """ 
    Fixed Settings
    """
    # internal settings
    int_highest_mapq = 60

    """
    Exit early when no samples is anlayzed
    """
    # if no samples will be analyzed, return
    if len(l_path_folder_output) == 0:
        logger.info(f"no output folders were given, exiting")
        return

    """
    Initiate pipelines for off-loading works
    """
    pipelines = bk.Offload_Works(
        None
    )  # no limit for the number of works that can be submitted.

    int_num_samples_analyzed_concurrently = min(
        len(l_path_folder_output), int_num_samples_analyzed_concurrently
    )  # if the number of samples are smaller than 'int_num_samples_analyzed_concurrently', adjust 'int_num_samples_analyzed_concurrently' so that it matches the number of samples

    n_threads = int(
        np.ceil(n_threads / int_num_samples_analyzed_concurrently)
    )  # divide the number of processes that can be used by each pipeline by the number of pipelines that will be run concurrently.
    
    """
    Load Minimap2 indices
    # 2023-07-31 13:42:38 
    """
    # load minimap2 aligners from index file paths if the aligners were not given.
    if am_genome is None :
        am_genome = mappy.Aligner( path_file_minimap_index_genome, preset = str_minimap_aligner_preset )
    if l_am_unwanted is None :
        l_am_unwanted = list( mappy.Aligner( e, preset = str_minimap_aligner_preset ) for e in l_path_file_minimap_index_unwanted )
        
    """
    Pipeline specific functions and variables
    """
    l_type_molecule = [ 'non_chimeric', 'chimeric' ] 
    if flag_recover_original_molecule_before_self_circularization_and_digestion : # add a 'type_molecule' if 'flag_recover_original_molecule_before_self_circularization_and_digestion' is True
        l_type_molecule += [ 'non_chimeric__before_self_circularization_and_digestion' ]
    l_name_output = [ 'aligned_to_unwanted_sequences', 'cannot_aligned_to_genome' ] # 'aligned_to_genome__non_chimeric__poly_A__plus_strand' : main fastq output file
    l_name_dist = [ 'aligned_to_unwanted_sequence', 'cannot_aligned_to_genome', 'aligned_to_genome' ] 
    for type_molecule in l_type_molecule :
        l_name_output.extend( [ f'aligned_to_genome__{type_molecule}__no_poly_A', f'aligned_to_genome__{type_molecule}__poly_A__plus_strand', ] )
        l_name_dist.extend( [ f'aligned_to_genome__{type_molecule}__no_poly_A', f'aligned_to_genome__{type_molecule}__external_poly_A__unreferenced_G', f'aligned_to_genome__{type_molecule}__internal_poly_A__unreferenced_G', f'aligned_to_genome__{type_molecule}__external_poly_A__no_unreferenced_G', f'aligned_to_genome__{type_molecule}__internal_poly_A__no_unreferenced_G', ] )
    
    def _write_a_fastq_record( newfastqfile, r ) :
        """ # 2023-08-01 12:21:30 
        write a fastq record to a given fastq output file (gzipped).
        """
        newfastqfile.write( ( r[ 0 ] + '\n' + r[ 1 ] + '\n+\n' + r[ 3 ] + '\n' ).encode( ) ) # the given record should contain '@' in front of the qname

    def _calculate_proportion_of_a_base( str_seq : str, str_base : str, int_min_length : int = 1 ) :
        """ # 2023-08-03 21:39:57 
        calculate a proportion of a base in a given sequence
        
        int_min_length : int # if the length of the given sequence is smaller than the given threshold, return np.nan as the proportion of the base. should be larger than 0 to avoid zero division error.
        """
        len_seq = len( str_seq ) # retrieve the length of the given sequence
        return [ len_seq, str_seq.count( str_base ) / len_seq ] if len_seq >= int_min_length else [ len_seq, np.nan ]
    
    def _classify_read( seq : str, q_st : int, q_en : int ) :
        """ # 2023-11-02 19:12:39 
        classlfy read by detecting poly A and unreferenced Gs.

        # possible labels
        'no_poly_A', 'external_poly_A__unreferenced_G', 'internal_poly_A__unreferenced_G', 'external_poly_A__no_unreferenced_G', 'internal_poly_A__no_unreferenced_G'

        return the direction of the read and the classification label
        """
        # internal setting
        float_min_G_frequency_for_identifying_unreferenced_Gs = 0.5
        int_size_window_for_identifying_unreferenced_Gs = 3

        # calculate metrics for classification
        len_left_clipped, left_clipped_T_prop = _calculate_proportion_of_a_base( seq[ max( 0, q_st - int_size_window_for_searching_poly_a_tail ) : q_st ], 'T', int_size_window_for_searching_poly_a_tail )
        len_left_internal, left_internal_T_prop = _calculate_proportion_of_a_base( seq[ q_st : q_st + int_size_window_for_searching_poly_a_tail ], 'T', int_size_window_for_searching_poly_a_tail )
        len_right_clipped, right_clipped_A_prop = _calculate_proportion_of_a_base( seq[ q_en : q_en + int_size_window_for_searching_poly_a_tail ], 'A', int_size_window_for_searching_poly_a_tail )
        len_right_internal, right_internal_A_prop = _calculate_proportion_of_a_base( seq[ max( 0, q_en - int_size_window_for_searching_poly_a_tail ) : q_en ], 'A', int_size_window_for_searching_poly_a_tail )
        len_left_clipped_3bp, left_clipped_3bp_G_prop = _calculate_proportion_of_a_base( seq[ max( 0, q_st - int_size_window_for_identifying_unreferenced_Gs ) : q_st ], 'G', int_size_window_for_identifying_unreferenced_Gs )
        len_right_clipped_3bp, right_clipped_3bp_C_prop = _calculate_proportion_of_a_base( seq[ q_en : q_en + int_size_window_for_identifying_unreferenced_Gs ], 'C', int_size_window_for_identifying_unreferenced_Gs )

        # retrieve flags for poly A
        flag_external_poly_A_flipped = left_clipped_T_prop >= float_min_A_frequency_for_identifying_poly_A # if length = 0, the value will be np.nan, and the comparison will be automatically failed
        flag_internal_poly_A_flipped = left_internal_T_prop >= float_min_A_frequency_for_identifying_poly_A
        flag_external_poly_A = right_clipped_A_prop >= float_min_A_frequency_for_identifying_poly_A
        flag_internal_poly_A = right_internal_A_prop >= float_min_A_frequency_for_identifying_poly_A
        # retrieve flags for unreferenced Gs
        flag_unreferenced_Gs = left_clipped_3bp_G_prop >= float_min_G_frequency_for_identifying_unreferenced_Gs
        flag_unreferenced_Gs_flipped = right_clipped_3bp_C_prop >= float_min_G_frequency_for_identifying_unreferenced_Gs

        ''' handles 'no_poly_A' '''
        if not( flag_external_poly_A_flipped or flag_internal_poly_A_flipped or flag_external_poly_A or flag_internal_poly_A ) :
            """
            [rebound] - search poly A that might be located farther downstream/upstream of the alignment end position using costly but more accurate search algorithm. The exact mechanism by which poly A can be found 5~20bp up/downstram of the alignment end position is currently not known.
            """
            int_min_poly_A_length = int( math.ceil( int_size_window_for_searching_poly_a_tail * float_min_A_frequency_for_identifying_poly_A ) ) # retrieve minimum length of a stretch of 'A' (or 'T') for identification of poly A tail.
            int_size_search_window = int_size_window_for_searching_poly_a_tail + int_max_size_intervening_sequence_between_alignment_end_and_poly_A # retrieve size of the search window

            # retry cacluation of flags for searching poly A
            flag_external_poly_A_flipped = len( STR.Find_stretch_of_a_character( seq[ max( 0, q_st - int_size_search_window ) : q_st ], 'T', int_len_threshold = int_min_poly_A_length ) ) > 0 # True if at least one stretch of 'T' (flipped) or 'A' exists in the search window.
            flag_internal_poly_A_flipped = len( STR.Find_stretch_of_a_character( seq[ q_st : q_st + int_size_search_window ], 'T', int_len_threshold = int_min_poly_A_length ) ) > 0
            flag_external_poly_A = len( STR.Find_stretch_of_a_character( seq[ q_en : q_en + int_size_search_window ], 'A', int_len_threshold = int_min_poly_A_length ) ) > 0
            flag_internal_poly_A = len( STR.Find_stretch_of_a_character( seq[ max( 0, q_en - int_size_search_window ) : q_en ], 'A', int_len_threshold = int_min_poly_A_length ) ) > 0

            ''' handles 'no_poly_A' (after retrying) '''
            if not( flag_external_poly_A_flipped or flag_internal_poly_A_flipped or flag_external_poly_A or flag_internal_poly_A ) :
                return 'no_poly_A', None # no direction for the 'no_poly_A' label

        ''' handles internal poly A '''
        if flag_internal_poly_A_flipped or flag_internal_poly_A :
            if flag_internal_poly_A_flipped :
                if flag_unreferenced_Gs_flipped :
                    return 'internal_poly_A__unreferenced_G', '-'
                else :
                    return 'internal_poly_A__no_unreferenced_G', '-'
            else :
                if flag_unreferenced_Gs :
                    return 'internal_poly_A__unreferenced_G', '+'
                else :
                    return 'internal_poly_A__no_unreferenced_G', '+'

        ''' handles external poly A '''
        if flag_external_poly_A_flipped :
            if flag_unreferenced_Gs_flipped :
                return 'external_poly_A__unreferenced_G', '-'
            else :
                return 'external_poly_A__no_unreferenced_G', '-'
        else :
            if flag_unreferenced_Gs :
                return 'external_poly_A__unreferenced_G', '+'
            else :
                return 'external_poly_A__no_unreferenced_G', '+'
    
    def _initialize_dict_arr_dist( ) :
        """ # 2023-08-03 11:49:26 
        initialize 'dict_arr_dist'
        """
        return dict( ( name_dist, None ) for name_dist in l_name_dist )
        
    def run_pipeline():
        """# 2023-10-03 20:00:57 
        analyze a pipeline for a given list of samples
        """
        # retrieve id of the pipeline
        str_uuid_pipeline = bk.UUID()
        logger.info(
            f"[Pipeline Start] Forked Pipeline (id={str_uuid_pipeline}) Started."
        )

        """
        Initiate workers for off-loading works
        """
        workers = bk.Offload_Works(
            None
        )  # no limit for the number of works that can be submitted.

        """
        Run pipeline for each sample
        """
        for path_file_fastq_input, path_folder_output in zip( l_path_file_fastq_input, l_path_folder_output ) :  # retrieve an output folder for the current sample
            """
            define a function to release a lock
            """
            def release_lock():
                """# 2023-01-14 20:36:17
                release the lock file
                """
                path_file_lock = (
                    f"{path_folder_output}ourotools.lock"
                )

                # check the existence of output files for the output folder of each input file of the current sample
                flag_all_output_files_exist = True  # initialize the flag
                
                if not os.path.exists(
                    f"{path_folder_output}pipeline_completed.txt"
                ):
                    flag_all_output_files_exist = False

                # check the existence of the lock file
                if (
                    os.path.exists(path_file_lock) and flag_all_output_files_exist
                ):  # if all output files exist and the lock file exists
                    # check whether the lock file has been created by the current pipeline
                    with open(path_file_lock, "rt") as file_lock:
                        str_uuid_pipeline_lock = file_lock.read() # retrieve uuid of lock
                        flag_lock_acquired = str_uuid_pipeline_lock == str_uuid_pipeline
                    if (
                        flag_lock_acquired
                    ):  # if the lock file has been created by the current pipeline, delete the lock file
                        os.remove(path_file_lock)
                        # lock has been released
                        if verbose:
                            logger.warning(
                                f"[{path_folder_output}] The forked pipeline (id={str_uuid_pipeline}) released the lock"
                            )
                    else :
                        # lock has been released
                        if verbose:
                            logger.warning(
                                f"[{path_folder_output}] The lock belongs to the forked pipeline (id={str_uuid_pipeline_lock}), and the lock was not released."
                            )
                else:
                    if verbose:
                        logger.warning(
                            f"[{path_folder_output}] The forked pipeline (id={str_uuid_pipeline}) attempted to release the lock, but some output files are missing, and the lock will not be released, yet."
                        )

            """
            Run pipeline for each sample
            """
            """
            create a lock
            """
            os.makedirs(path_folder_output, exist_ok=True)
            path_file_lock = (
                f"{path_folder_output}ourotools.lock"
            )
            # check the existence of the lock file
            if os.path.exists(path_file_lock):
                logger.warning(
                    f"[Output folder unavailable] the output folder {path_folder_output} contains a lock file, which appears to be processed by a different process. Therefore, the output folder will be skipped."
                )
                continue
            flag_lock_acquired = False  # initialize 'flag_lock_acquired'
            try:
                # create the lock file
                with open(path_file_lock, "wt") as newfile_lock:
                    newfile_lock.write(str_uuid_pipeline)
                # check whether the lock file has been created correctly (check for collision).
                with open(path_file_lock, "rt") as file_lock:
                    flag_lock_acquired = file_lock.read() == str_uuid_pipeline
            except Exception as e:
                logger.critical(
                    e, exc_info=True
                )  # if an exception occurs, print the error message
            if not flag_lock_acquired:
                logger.warning(
                    f"[Output folder unavailable] an attempt to acquire a lock for the output folder {path_folder_output} failed, which appears to be processed by a different process. Therefore, the output folder will be skipped."
                )
                continue
            # lock has been acquired

            """
            Run pipeline for each input file
            """
            # define folders and directories
            path_file_fastq_input = os.path.abspath(path_file_fastq_input)
            if path_folder_output is None:  # set default 'path_folder_output'
                path_folder_output = (
                    f"{path_file_fastq.rsplit( '/', 1 )[ 0 ]}LongFilterNSplit_output/"
                )
            path_folder_output = os.path.abspath(path_folder_output)
            path_folder_output += "/"
            path_folder_temp = f"{path_folder_output}temp/"
            path_folder_graph = f"{path_folder_output}graph/"

            """ if the output folder already exists """
            if os.path.exists(path_folder_output):
                """check whether the pipeline has been completed"""
                if os.path.exists( f"{path_folder_output}pipeline_completed.txt" ) :  # intermediate files should not exists, while all output files should exist
                    logger.info(
                        f"[Output folder Already Exists] the output folder {path_folder_output} contains valid output files. Therefore, the output folder will be skipped."
                    )
                    release_lock( ) # release the lock
                    continue  # skip if the pipeline has been completed for the output folder
                else:
                    """if required output files does not exist or the an intermediate file exists, remove the entire output folder, and rerun the pipeline"""
                    if (
                        len(glob.glob(f"{path_folder_output}*/")) > 0
                    ):  # detect a folder inside the output folder and report the presence of the existing folders.
                        logger.info(
                            f"[Output folder Already Exists] the output folder {path_folder_output} does not contain valid output files. The output folder will be cleaned and the pipeline will start anew at the folder."
                        )
                    # delete the folders
                    for path_folder in glob.glob(f"{path_folder_output}*/"):
                        shutil.rmtree(path_folder, ignore_errors = True)
                    # delete the files, excluding the lock file that has been acquired by the current pipeline
                    for path_file in glob.glob(f"{path_folder_output}*"):
                        if (
                            path_file_lock != path_file
                        ):  # does not delete the lock file
                            os.remove(path_file)

            """ create directories """
            for path_folder in [
                path_folder_output,
                path_folder_temp,
                path_folder_graph,
            ]:
                os.makedirs(path_folder, exist_ok=True)

            """
            Report program arguments
            """
            # record arguments used for the program (metadata)
            dict_program_setting = {
                "version": _version_,  # record version
                # external
                "flag_usage_from_command_line_interface" : flag_usage_from_command_line_interface,
                "path_file_minimap_index_genome" : path_file_minimap_index_genome,
                "path_file_fastq_input" : path_file_fastq_input,
                "path_folder_output" : path_folder_output,
                "n_threads" : n_threads,
                "int_num_samples_analyzed_concurrently" : int_num_samples_analyzed_concurrently,
                "float_memory_in_GiB" : float_memory_in_GiB,
                "int_num_base_pairs_in_a_batch" : int_num_base_pairs_in_a_batch,
                'str_minimap_aligner_preset' : str_minimap_aligner_preset,
                'int_min_mapq' : int_min_mapq,
                'int_size_window_for_searching_poly_a_tail' : int_size_window_for_searching_poly_a_tail,
                'float_min_A_frequency_for_identifying_poly_A' : float_min_A_frequency_for_identifying_poly_A,
                'int_max_read_length' : int_max_read_length,
                'flag_recover_original_molecule_before_self_circularization_and_digestion' : flag_recover_original_molecule_before_self_circularization_and_digestion,
                'int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion' : int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion,
                # internal
                "path_folder_temp": path_folder_temp,
                "path_folder_graph": path_folder_graph,
            }
            logger.info(
                f"[Setting] program will be run with the following setting for the input file {path_file_fastq_input} : {str( dict_program_setting )}"
            )

            """ export program setting """
            path_file_json_setting_program = (
                f"{path_folder_output}program_setting.json"
            )
            if os.path.exists(path_file_json_setting_program):
                with open(path_file_json_setting_program, "r") as file:
                    j = json.load(file)
                if j != dict_program_setting:
                    logger.info(
                        f"[Warning] the current program setting is different from the previous program setting recorded in the pipeline folder. The previous setting will be used."
                    )
                    with open(path_file_json_setting_program, "r") as file:
                        dict_program_setting = json.load(
                            file
                        )  # override current program setting with previous program setting
            with open(path_file_json_setting_program, "w") as newfile:
                json.dump(dict_program_setting, newfile)
                
            """
            Define a generator for partitioning input file
            """

            def gen_batch( ):
                """# 2023-07-30 18:37:49 
                create batch from the input fastq file
                """
                int_base_pair_counter = 0 # initialize base pair counter
                l_r_for_a_batch = [ ] # a list of records for a batch
                for r in bk.FASTQ_Iterate( path_file_fastq_input ) : # iterate through the input FASTQ file
                    l_r_for_a_batch.append( r ) # add the record
                    int_base_pair_counter += len( r[ 1 ] ) # increase the counter
                    if int_base_pair_counter >= int_num_base_pairs_in_a_batch : # if the batch is full, yield the batch
                        yield l_r_for_a_batch
                        l_r_for_a_batch = [ ] # initialize the next batch
                        int_base_pair_counter = 0 # initialize base pair counter
                if len( l_r_for_a_batch ) > 0 : # if records are remaining in the list, yield the list as the last batch
                    yield l_r_for_a_batch

            def process_batch(pipe_receiver, pipe_sender):
                """
                # 2022-04-24 01:29:59
                Requires loading several data objects (using copy-on-write method)

                receives a bookmark file (either file directory of a tsv file or a dataframe)
                """
                """
                initialize the worker 
                # 2023-08-01 12:19:06 
                """
                str_uuid = bk.UUID()  # retrieve id
                if verbose:
                    logger.info(f"[Started] start working (worker_id={str_uuid})")
                
                """ open output files """
                str_uuid_for_a_batch = bk.UUID( ) # retrieve id for the specific batch
                dict_newfile_fastq_output = dict( ( name_output, gzip.open( f"{path_folder_temp}{str_uuid}.{name_output}.fastq.gz", "wb", ) ) for name_output in l_name_output )
                
                while True:
                    ins = pipe_receiver.recv()
                    if ins is None:
                        break
                    l_r_for_a_batch = ins  # parse input
                    
                    """
                    Filter reads that were aligned to unwanted sequences.
                    # 2023-07-31 13:40:24 
                    """
                    int_total_num_records_for_a_batch = len( l_r_for_a_batch ) # record the total number of records

                    # initialize summary metrics
                    dict_arr_dist = _initialize_dict_arr_dist( ) # initialize 'dict_arr_dist'

                    """
                    define batch-specific function
                    """
                    def _process_molecule( qname : str, seq : str, qual : str, hit, st : Union[ None, int ] = 0, en : Union[ None, int ] = None, type_molecule : Literal[ 'non_chimeric', 'chimeric', 'non_chimeric__before_self_circularization_and_digestion' ] = 'non_chimeric' ) :
                        """ # 2023-11-02 20:12:33 
                        process a molecule (segment) of a sequence

                        hit, # a Minimap2 mappy alignment record. multiple hits can be given as a list of hits. 
                        type_molecule : Literal[ 'non_chimeric', 'chimeric', 'non_chimeric__before_self_circularization_and_digestion' ] = 'non_chimeric' # type of the molecule
                        """
                        if type_molecule == 'non_chimeric__before_self_circularization_and_digestion' :
                            """ analyze circular molecule """
                            '''
                            recover original molecule before self-circularization and digestion 

                            q_st q_en r_st r_en strand
                            1409 1762 4595 4948 1 
                            1892 2576 3913 4597 1 

                            OR

                            q_st q_en r_st  r_en  strand
                            137  2000 11741 13619 -1 
                            2134 2252 13619 13738 -1 
                            '''
                            hit_1, hit_2 = l_hit_segment # parse 'l_hit_segment'
                            q_st_1, q_en_1, q_st_2, q_en_2 = hit_1.q_st, hit_1.q_en, hit_2.q_st, hit_2.q_en # retrieve properties
                            seq_intervening, qual_intervening = seq[ q_en_1 : q_st_2 ], qual[ q_en_1 : q_st_2 ] # retrieve intervening sequence (cyclized adaptors)
                            len_seq_intervening = len( seq_intervening )
                            seg = seq_intervening + seq[ q_st_2 : q_en_2 ] + seq[ q_st_1 : q_en_1 ] + seq_intervening # composing original sequence 
                            qual_seg = qual_intervening + qual[ q_st_2 : q_en_2 ] + qual[ q_st_1 : q_en_1 ] + qual_intervening
                            qname_suffix = '_circ' + str( q_st_1 ) + 'to' + str( q_en_2 ) # compose suffix
                            q_st, q_en = len_seq_intervening, len( seg ) - len_seq_intervening # compose 'q_st' and 'q_en'
                        else :
                            """ analyze linear molecule """
                            ''' Retrieve the segment and the length of the segment. Also, compose suffix to the qname '''
                            if st == 0 and en is None : # analyze an entire molecule
                                seg, qual_seg, qname_suffix = seq, qual, '' 
                            elif st == 0 :
                                seg, qual_seg, qname_suffix = seq[ : en ], qual[ : en ], '_to' + str( en )
                            elif en is None :
                                seg, qual_seg, qname_suffix = seq[ st : ], qual[ st : ], '_' + str( st ) + 'to'
                            else :
                                seg, qual_seg, qname_suffix = seq[ st : en ], qual[ st : en ], '_' + str( st ) + 'to' + str( en )

                            ''' compose 'q_st' and 'q_en' (simply substract 'st') '''
                            if isinstance( hit, list ) : # if a list of hits were given, 
                                q_st, q_en = min( h.q_st for h in hit ) - st, max( h.q_en for h in hit ) - st # retrieve smallest q_st and largest q_en for q_st and q_en of the segment if multiple hits were given
                            else : # when single hit was given
                                q_st, q_en = hit.q_st - st, hit.q_en - st

                        # retrieve length of the segment
                        len_seg = len( seg )

                        label, direction = _classify_read( seg, q_st, q_en ) # classify the segment 

                        dict_arr_dist[ f'aligned_to_genome__{type_molecule}__{label}' ] = _update_size_distribution( new_size = len_seg, arr_dist = dict_arr_dist[ f'aligned_to_genome__{type_molecule}__{label}' ] ) # update appropriate distribution of reads using the label

                        if label == 'no_poly_A' : # (likely to be not analyzed)
                            _write_a_fastq_record( dict_newfile_fastq_output[ f'aligned_to_genome__{type_molecule}__no_poly_A' ], [ "@" + qname + qname_suffix, seg, '+', qual_seg ] ) # write the current read to the appropriate output fastq file
                        else : # collect all reads with poly A to 'poly_A__plus_strand' output file. (likely to be analyzed together)
                            """
                            if needed, modify the fastq record so that poly A can be located at the right, representing the '+' strand of the original mRNA initially captured by the primer.
                            """
                            if direction == '-' :
                                seg = SEQ.Reverse_Complement( seg ) # reverse complement the sequence
                                qual_seg = qual_seg[ : : -1 ] # flip the quality scores
                            _write_a_fastq_record( dict_newfile_fastq_output[ f'aligned_to_genome__{type_molecule}__poly_A__plus_strand' ], [ "@" + qname + qname_suffix + '_R', seg, '+', qual_seg ] ) # write the current read to the appropriate output fastq file # add additional suffix to show the sequence has been reverse complemented.

                    for r in l_r_for_a_batch :
                        header, seq, _, qual = r # parse fastq record
                        len_seq = len( seq ) # retrieve length of the sequence
                        qname = header.split( ' ', 1 )[ 0 ][ 1 : ] # retrieve qname
                        
                        """
                        handle the case when read length exceed the given limit (will be considered as 'cannot be aligned to genome')
                        """
                        if len_seq > int_max_read_length :
                            dict_arr_dist[ 'cannot_aligned_to_genome' ] = _update_size_distribution( new_size = len_seq, arr_dist = dict_arr_dist[ 'cannot_aligned_to_genome' ] ) # update distribution of reads that cannot be aligned to the genome
                            _write_a_fastq_record( dict_newfile_fastq_output[ 'cannot_aligned_to_genome' ], r ) # write the current read to the appropriate output fastq file
                            continue # skip the remaining operations
                        
                        """
                        align read to the list of unwanted sequences
                        """
                        flag_aligned_to_unwanted_sequence = False
                        for am_unwanted in l_am_unwanted : # for each aligner for unwanted sequences
                            l_hit_unwanted = list( hit for hit in am_unwanted.map( seq ) ) # exhuast the iterator to avoid the memory leak
                            if len( l_hit_unwanted ) > 0 :
                                flag_aligned_to_unwanted_sequence = True
                                break
                            for hit in l_hit_unwanted :
                                l_seq, int_total_aligned_length = bk.SAM.Retrive_List_of_Mapped_Segments( hit.cigar, hit.r_st, flag_is_cigartuples_from_mappy = True )
                                
                        """
                        handle the case when read was aligned to unwanted sequences
                        """
                        if flag_aligned_to_unwanted_sequence :
                            dict_arr_dist[ 'aligned_to_unwanted_sequence' ] = _update_size_distribution( new_size = len_seq, arr_dist = dict_arr_dist[ 'aligned_to_unwanted_sequence' ] ) # update distribution of reads aligned to unwanted sequences
                            _write_a_fastq_record( dict_newfile_fastq_output[ 'aligned_to_unwanted_sequences' ], r ) # write the current read to the appropriate output fastq file
                            continue # skip the remaining operations
                            
                        """
                        align the read to genome
                        """
                        l_hit_genome = list( hit for hit in am_genome.map( seq ) if hit.mapq >= int_min_mapq ) # exhuast the iterator to avoid the memory leak # filter hits using mapping quality

                        """
                        handle the case when read was not aligned to genome
                        """
                        if len( l_hit_genome ) == 0 :
                            dict_arr_dist[ 'cannot_aligned_to_genome' ] = _update_size_distribution( new_size = len_seq, arr_dist = dict_arr_dist[ 'cannot_aligned_to_genome' ] ) # update distribution of reads that cannot be aligned to the genome
                            _write_a_fastq_record( dict_newfile_fastq_output[ 'cannot_aligned_to_genome' ], r ) # write the current read to the appropriate output fastq file
                            continue # skip the remaining operations
                        
                        """
                        analyze the alignments to the genome
                        """
                        dict_arr_dist[ 'aligned_to_genome' ] = _update_size_distribution( new_size = len_seq, arr_dist = dict_arr_dist[ 'aligned_to_genome' ] ) # update distribution of reads aligned to genome

                        """
                        handle the case when read was aligned to genome only once (non-chimeric read, the majority of cases)
                        """
                        if len( l_hit_genome ) == 1 :
                            _process_molecule( qname, seq, qual, l_hit_genome[ 0 ] ) # process non-chimeric segment
                            continue # skip the remaining operations
                            
                        """
                        handle reads with multiple genome alignments
                        """
                        l_l = [ ] # initialize the container # 
                        for hit in l_hit_genome :
                            l_l.append( [ hit.q_st, hit ] )
                            # l_seq, int_total_aligned_length = bk.SAM.Retrive_List_of_Mapped_Segments( hit.cigar, hit.r_st, flag_is_cigartuples_from_mappy = True )
    
                        arr_algn = np.array( l_l, dtype = object ) # create an array of alignments
                        arr_algn = arr_algn[ arr_algn[ :, 0 ].argsort( ) ] # sort alignments using 'q_st'
                        
                
                        """
                        Split a read into multiple segments (multiple output reads)
                        """
                        # initialize the search
                        flag_is_segment_chimeric, q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev, strand_prev = False, 0, None, None, None, None, None # initialize a flag indicating the segment is chimeric or non-chimeric # set 'q_st_segment' as 0
                        l_hit_segment = [ ] # initialize a list of hit of a segment
                        for q_st, hit in arr_algn : # for each q_st and hit
                            q_en = hit.q_en # retrieve q_en
                            ''' initialize the segment with current alignment (if it was not initialized) '''
                            if q_en_segment is None :
                                q_en_segment, ctg_prev, r_st_prev, r_en_prev, strand_prev = hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.strand # use the start of the molecule as 'q_en_segment'
                                l_hit_segment.append( hit )
                                continue # continue to the next hit

                            """
                            recover original molecule before self-circularization by identifying the unique patterns:

                            q_st q_en r_st r_en strand
                            1409 1762 4595 4948 1 
                            1892 2576 3913 4597 1 

                            OR

                            q_st q_en r_st  r_en  strand
                            137  2000 11741 13619 -1 
                            2134 2252 13619 13738 -1 
                            """
                            if flag_recover_original_molecule_before_self_circularization_and_digestion and ( len( l_hit_segment ) == 1 ) and ( ctg_prev == hit.ctg ) and ( strand_prev == hit.strand ) and ( ( abs( r_st_prev - hit.r_en ) if strand_prev > 0 else abs( r_en_prev - hit.r_st ) ) <= int_max_coordinate_difference_at_predicted_cut_site_for_recovering_original_molecule_before_self_circularization_and_digestion ) : # there should be a single previous segment, the strand and contig of the previous segment should be matched with those of the current segment in order to be called as a pair of seperated segments of a self-circularized molecule. # if the coordinate difference between predicted cut sites is below the given threshold, call the pair of segments as a self-circularized molecule.
                                l_hit_segment.append( hit ) # add current alignment to the segment
                                _process_molecule( qname, seq, qual, l_hit_segment, q_st_segment, q_en_segment, 'non_chimeric__before_self_circularization_and_digestion' ) # not including chimeric reads (all recovered molecule are considered non-chimeric) # process 'non_chimeric__before_self_circularization_and_digestion' molecule
                                # flush the data (segmant will start anew)
                                flag_is_segment_chimeric, q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev, strand_prev = False, q_en, None, None, None, None, None # set q_st_segment as 'q_en' to include a flanking sequence
                                l_hit_segment = [ ] 
                                continue            

                            """ split the read """
                            int_size_gap = q_st - q_en_segment # calculate the size of the intervening sequence
                            if int_size_gap >= int_min_size_intervening_sequence_for_splitting : # if the size of the intervening sequences is larger then the threshold, split the read
                                int_size_flanking = min( int_size_gap, int_min_size_intervening_sequence_for_splitting ) # retrieve the size of the flanking sequence to include in the segment. the size of the flanking sequence cannot be larger than the intervening sequence
                                _process_molecule( qname, seq, qual, l_hit_segment, q_st_segment, q_en_segment + int_size_flanking, 'chimeric' if flag_is_segment_chimeric else 'non_chimeric' ) # process chimeric segment # add 'int_size_flanking' to q_en_segment to include a flanking sequence
                                # initialize the next segment
                                flag_is_segment_chimeric, q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev, strand_prev = False, q_st - int_size_flanking, hit.q_en, hit.ctg, hit.r_st, hit.r_en, hit.strand # set q_st_segment as q_st - int_size_flanking to include a flanking sequence
                                l_hit_segment = [ hit ] 
                                continue

                            """ concatenate genomic alignments and determine whether the segment is chimeric or not. (regardless of whether the segment is chimeric or not, the genomic alignment will be concatenated.) """
                            # identify chimeric molecule
                            if ctg_prev != hit.ctg :
                                flag_is_segment_chimeric = True
                            elif max( r_st_prev - hit.r_en, hit.r_st - r_en_prev ) > int_max_intron_size_for_determining_chimeric_molecule : # if the distance between alignment is longer than maximum intron size, consider reads as an intra-chromosomal chimeric molecule
                                flag_is_segment_chimeric = True
                            # extend segment
                            l_hit_segment.append( hit )
                            q_en_segment = max( q_en_segment, q_en ) # update 'q_en_segment' ('q_st_segment' will not change)
                            ctg_prev, r_st_prev, r_en_prev, strand_prev = hit.ctg, hit.r_st, hit.r_en, hit.strand # update the previous alignment (for the next segment, determining whether the segment is also part of the chimeric molecule will depend on the alignment information of the current segment.)

                        if len( l_hit_segment ) > 0 : # if a segment is remaining, process the segment
                            _process_molecule( qname, seq, qual, l_hit_segment, q_st_segment, None, 'chimeric' if flag_is_segment_chimeric else 'non_chimeric' ) # process chimeric segment # use the end of the molecule as 'q_en_segment'

                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : int_total_num_records_for_a_batch,
                        'dict_arr_dist' : dict_arr_dist,
                    } )  # report the number of processed records
                    """ report the worker has completed a batch """
                    if verbose:
                        logger.info(f"[Completed] completed a batch (worker_id={str_uuid})")
                    
                """ close output files """
                for name_type in dict_newfile_fastq_output :
                    dict_newfile_fastq_output[ name_type ].close( )
                    
                """ report the worker has completed all works """
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")
                pipe_sender.send( 'completed' )  

            name_file_input = path_file_fastq_input.rsplit( '/', 1 )[ 1 ] if '/' in path_file_fastq_input else path_file_fastq_input # retrieve name of the input file
            ns = dict()  # define a namespace
            ns[ "int_num_read_currently_processed" ] = 0  # initialize total number of reads processed by the algorithm
            ns[ 'dict_arr_dist' ] = _initialize_dict_arr_dist( ) # initialize 'dict_arr_dist'

            def post_process_batch(res):
                # parse received result
                int_total_num_records_for_a_batch = res[ 'int_total_num_records_for_a_batch' ]
                ns["int_num_read_currently_processed"] += int_total_num_records_for_a_batch
                if verbose :
                    logger.info( f"[{path_file_fastq_input}] a batch has been completed, {0 if res[ 'dict_arr_dist' ][ 'aligned_to_unwanted_sequence' ] is None else res[ 'dict_arr_dist' ][ 'aligned_to_unwanted_sequence' ].sum( )}/{int_total_num_records_for_a_batch} number of reads were aligned to unwanted sequences, {0 if res[ 'dict_arr_dist' ][ 'cannot_aligned_to_genome' ] is None else res[ 'dict_arr_dist' ][ 'cannot_aligned_to_genome' ].sum( )}/{int_total_num_records_for_a_batch} number of reads cannot be aligned to genome" )
                    logger.info( f"[{path_file_fastq_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed." )  # report
                
                # combine distributions
                for name_cat_dist in _initialize_dict_arr_dist( ) : # for each category
                    ns[ 'dict_arr_dist' ][ name_cat_dist ] = _combine_size_distribution( ns[ 'dict_arr_dist' ][ name_cat_dist ], res[ 'dict_arr_dist' ][ name_cat_dist ] ) # combine and update the global distributions
                    
            """
            Analyze an input file
            """
            if verbose:
                logger.info(
                    f"[{path_file_fastq_input}] the analysis pipeline will be run with {n_threads} number of threads"
                )
            bk.Multiprocessing_Batch_Generator_and_Workers( 
                gen_batch=gen_batch(),
                process_batch=process_batch,
                post_process_batch=post_process_batch,
                int_num_threads=n_threads
                + 2,  # one thread for generating batch, another thread for post-processing of the batch
                flag_wait_for_a_response_from_worker_after_sending_termination_signal = True, # wait until all worker exists before resuming works in the main process
            )

            """ 
            post-processing
            """

            def post_processing():  # off-loading a single-core work
                logger.info(
                    f"[{path_file_fastq_input}] post-processing started"
                )
                # combine results into a single output file (initial read analysis)
                for name_output in l_name_output :
                    bk.OS_Run( [ 'cat' ] + glob.glob( f"{path_folder_temp}*.{name_output}.fastq.gz" ), stdout_binary = True, path_file_stdout = f"{path_folder_output}{name_output}.fastq.gz" )
                    
                ''' summarize distributions '''
                dict_arr_dist = ns[ 'dict_arr_dist' ] # retrieve 'dict_arr_dist'
                l_l = [ ]
                for e in dict_arr_dist :
                    arr = dict_arr_dist[ e ]
                    if arr is None :
                        continue
                    arr_bins = np.arange( len( arr ) ) # retrieve bin size of the histograms
                    int_num_reads, int_num_base_pairs = arr.sum( ), ( arr * arr_bins ).sum( )
                    int_avg_length_base_pairs = int_num_base_pairs / int_num_reads
                    float_standard_deviation_length_base_pairs = np.sqrt( np.average((arr_bins - int_avg_length_base_pairs)**2, weights=arr) )
                    l_l.append( [ e, int_num_reads, int_num_base_pairs, int_avg_length_base_pairs, float_standard_deviation_length_base_pairs ] )
                df_summary_of_distributions = pd.DataFrame( l_l, columns = [ 'name_type_distribution', 'int_num_reads', 'int_num_base_pairs', 'int_avg_length_base_pairs', 'float_standard_deviation_length_base_pairs' ] )
                df_summary_of_distributions.to_csv( f"{path_folder_output}df_summary_of_distributions.tsv.gz", sep = '\t', index = False ) # export 'df_summary_of_distributions'
                
                """
                Draw plots of distributions
                """
                # create output folders
                path_folder_graph_noninteractive, path_folder_graph_interactive = f"{path_folder_graph}noninteractive_graph/", f"{path_folder_graph}interactive_graph/"
                for path_folder in [ path_folder_graph_noninteractive, path_folder_graph_interactive ] :
                    os.makedirs( path_folder, exist_ok = True )

                ''' draw simple line plots '''
                # plot settings
                int_max_molecule_size_plot = 6500
                for name_cat_dist in _initialize_dict_arr_dist( ) : # for each category
                    if dict_arr_dist[ name_cat_dist ] is not None :
                        len_max_molecule_size_data = len( dict_arr_dist[ name_cat_dist ] ) # retrieve max molecule size 
                        plt.plot( np.arange( min( int_max_molecule_size_plot, len_max_molecule_size_data ) ), dict_arr_dist[ name_cat_dist ] if len_max_molecule_size_data <= int_max_molecule_size_plot else dict_arr_dist[ name_cat_dist ][ : int_max_molecule_size_plot ] )
                        plt.title( f"{name_cat_dist} ({dict_arr_dist[ name_cat_dist ].sum( )} molecules)" )
                        bk.MPL_SAVE( f"{name_cat_dist}.distribution", folder = path_folder_graph_noninteractive, l_format=['.pdf', '.png'] )
                        
                ''' draw interactive stacked bar graphs '''
                df_bar = _get_df_bar( dict_arr_dist, int_size_bin_in_base_pairs = 50, int_max_size_in_base_pairs = int_max_molecule_size_plot ) # retrieve a dataframe for drawing a bar graph
                
                for flag_use_proportion in [ True, False ] :
                    _draw_bar_plot( 
                        df_bar, 
                        [ 'aligned_to_unwanted_sequence', 'cannot_aligned_to_genome',  'aligned_to_genome', ],
                        title = f"Alignment Results of '{name_file_input}'",
                        flag_use_proportion = flag_use_proportion,
                        flag_save_figure = True, path_folder_graph = path_folder_graph_interactive,
                    )
                    _draw_bar_plot( 
                        df_bar, l_name_dist[ : 2 ] + l_name_dist[ 3 : ], # drop 'aligned_to_genome' from the result
                        title = f"Adaptor-Based Classification Results of '{name_file_input}'",
                        flag_use_proportion = flag_use_proportion,
                        flag_save_figure = True, path_folder_graph = path_folder_graph_interactive,
                    )
                
                ''' export pickle files '''
                # write distribution data as a pickle file
                bk.PICKLE_Write( f"{path_folder_output}dict_arr_dist.pkl", ns[ 'dict_arr_dist' ] )
                    
#                 pd.DataFrame( ns[ "l_l" ], columns = [ 'qname', 'len_seq', 'mapq', 'q_st', 'q_en', 'ref_name', 'ref_st', 'ref_en', 'strand' ] ).to_csv( f"{path_folder_output}output.tsv.gz", sep = '\t', index = False ) # ❤️ 'len_left_clipped', 'left_clipped_T_prop', 'len_left_internal', 'left_internal_T_prop', 'len_right_clipped', 'right_clipped_A_prop', 'len_right_internal', 'right_internal_A_prop', 'len_left_clipped_3bp', 'right_left_clipped_3bp_G_prop', 'len_right_clipped_3bp', 'right_clipped_3bp_C_prop'
                
                # write a flag indicating that the processing has been completed
                with open( f"{path_folder_output}pipeline_completed.txt", 'w' ) as newfile :
                    newfile.write( 'completed' )

                # delete temporary files
                shutil.rmtree( path_folder_temp, ignore_errors = True )
                    
                release_lock()  # release the lock
                logger.info(
                    f"[{path_file_fastq_input}] post-processing completed"
                )

            workers.submit_work(post_processing)

            release_lock()  # release the lock

        # wait all the single-core works offloaded to the workers to be completed.
        workers.wait_all()
        logger.info(
            f"[Pipeline Completion] Forked Pipeline (id={str_uuid_pipeline}) Completed."
        )

    for _ in range(
        int_num_samples_analyzed_concurrently
    ):  # run 'int_num_samples_analyzed_concurrently' number of pipelines
        pipelines.submit_work(run_pipeline)

    # wait all pipelines to be completed
    pipelines.wait_all()
    logger.info(f"Completed.")
    return 

def LongExtractBarcodeFromBAM(
    flag_usage_from_command_line_interface: bool = False, # a flag indicating the usage in the command line
    l_path_file_bam_input: Union[list, None] = None, # list of input BAM files
    l_path_folder_output: [list[str], None] = None, # list of output folders
    n_threads: int = 32, # the number of threads to use
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_base_pairs_in_a_batch : int = 2_500_000, # the number of base pairs in a batch
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to filter read with low alignment quality
    float_memory_in_GiB : float = 50, # expected memory usage of the pipeline
    float_error_rate : float = 0.2, # maximum error rate to consider when searching adaptor sequence in the read
    int_length_cb : int = 16, # the length of the cell barcode
    int_length_umi : int = 12, # the length of the UMI (unique molecular identifier)
    str_seq_r1 : str = 'CTACACGACGCTCTTCCGATCT', # the sequence of R1 adaptor (in 10x GEX v3 kit, located upstream of CB and UMI)
    str_seq_tso : str = 'AAGCAGTGGTATCAACGCAGAG', # the sequence of TSO adaptor (in 10x GEX v3 kit, located at 5' end of the molecule)
    path_file_valid_cb : str = None, # (required argument) the path to tsv file of whitelist barcodes. For more details, please see 10x cellranger references.
    int_max_num_cell_expected : int = 20_000, # the max number of expected cells
    int_len_sliding_window_internal_polyT : int = 10, # the length of sliding window for searching internal poly T (poly A) tract. (When poly-A tailed read is reverse complemented, R1 adaptor become situated in the forward direction
    int_len_window_internal_polyT : int = 30, # the size of window for searching for internal poly T
    float_min_T_fraction : float = 0.8, # the minimum T fraction for identifying the stretch of poly T tract
    int_min_n_overlap_kmer_for_clustering_umi : int = 1, # the minimum number of overlapped kmer for initiating UMI clustering 
    int_len_kmer_for_clustering_umi : int = 7, # the length of kmer for clustering UMI
    float_min_proportion_read_to_select_kmer_representative_for_clustering_umi : float = 0.75, # if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.
    int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level : int = 50, # the size of the bin (in base pairs) for collecting size distributions at the single-cell level
    verbose: bool = True,
) -> None :
    """# 2023-10-03 23:35:54 
    
    flag_usage_from_command_line_interface: bool = False, # a flag indicating the usage in the command line
    l_path_file_bam_input: Union[list, None] = None, # list of input BAM files
    l_path_folder_output: [list[str], None] = None, # list of output folders
    n_threads: int = 32, # the number of threads to use
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_base_pairs_in_a_batch : int = 2_500_000, # the number of base pairs in a batch
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to filter read with low alignment quality
    float_memory_in_GiB: float = 50,
    float_error_rate : float = 0.2, # maximum error rate to consider when searching adaptor sequence in the read
    int_length_cb : int = 16, # the length of the cell barcode
    int_length_umi : int = 12, # the length of the UMI (unique molecular identifier)
    str_seq_r1 : str = 'CTACACGACGCTCTTCCGATCT', # the sequence of R1 adaptor (in 10x GEX v3 kit, located upstream of CB and UMI)
    str_seq_tso : str = 'AAGCAGTGGTATCAACGCAGAG', # the sequence of TSO adaptor (in 10x GEX v3 kit, located at 5' end of the molecule)
    path_file_valid_cb : str = None, # (required argument) the path to tsv file of whitelist barcodes. For more details, please see 10x cellranger references.
    int_max_num_cell_expected : int = 20000, # the max number of expected cells
    int_len_sliding_window_internal_polyT : int = 10, # the length of sliding window for searching internal poly T (poly A) tract. (When R1 adaptor become situated in the forward direction, poly-A tailed read is reverse complemented to become poly-T containing read)
    int_len_window_internal_polyT : int = 30, # the size of window for searching for internal poly T
    float_min_T_fraction : float = 0.8, # the minimum T fraction for identifying the stretch of poly T tract
    int_min_n_overlap_kmer_for_clustering_umi : int = 1, # the minimum number of overlapped kmer for initiating UMI clustering 
    int_len_kmer_for_clustering_umi : int = 7, # the length of kmer for clustering UMI
    float_min_proportion_read_to_select_kmer_representative_for_clustering_umi : float = 0.75, # if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.
    int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level : int = 50, # the size of the bin (in base pairs) for collecting size distributions at the single-cell level
    verbose: bool = True,
    
    returns None

    --------------------------------------------------------------------------------------------------------
    Ouro-Tools
    BAM tags description
    # 2023-08-13 17:59:06 

    'CB', 'Z' : corrected cell barcode sequence (cell barcode correction at sample level, using all reads from a BAM)
    'UB', 'Z' : corrected UMI sequence using UMI clustering (UMI clustering was performed using all reads with the same each poly(A)-tail attached site, considering the direction of the transcript)
    'UR', 'Z' : uncorrected UMI sequence before UMI clustering (a subsequence of the raw sequence from the 'CU' tag).
    'XR', 'i' : the number of errors for identification of R1 adaptor (marks the 3' end of the original cDNA molecule for 10x GEX 3' products, where cell barcode and UMI sequences can be found). -1 indicates that the adaptor was not identified.
    'XT', 'i' : the number of errors for identification of TSO adaptor (marks the 5' end of the original cDNA molecule). -1 indicates that the adaptor was not identified.
    'CU', 'Z' : the uncorrected raw CB-UMI sequence before processing
    'IA', 'i' : the length of detected internal poly(A) priming region in the genomic alignment. 
    'LE', 'i' : the total length of genomic regions that are actually covered by the read, excluding spliced introns (the sum of exons).
    """
    """
    Parse arguments
    """
    if flag_usage_from_command_line_interface:  # parse arguments
        """parse arguments when the function was called from the command-line interface"""
        # {  } # unused arguments
        # command line arguments
        parser = argparse.ArgumentParser(
            description=str_description,
            usage="ourotools LongExtractBarcodeFromBAM",
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument("LongExtractBarcodeFromBAM")

        arg_grp_general = parser.add_argument_group("General")
        arg_grp_general.add_argument(
            "-l",
            "--l_path_file_bam_input",
            help="",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-o",
            "--l_path_folder_output",
            help="",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-t",
            "--n_threads",
            help="(default: 32) the number of processors to use concurrently.",
            default=32,
            type=int,
        )
        arg_grp_general.add_argument(
            "-b",
            "--int_num_base_pairs_in_a_batch",
            help="(default: 2,500,000) the number of base pairs in a batch.",
            default=2_500_000,
            type=int,
        )
        arg_grp_general.add_argument(
            "-s",
            "--int_num_samples_analyzed_concurrently",
            help="(default: 2) the number of samples that can be analyzed concurrently.",
            default=2,
            type=int,
        )
        arg_grp_general.add_argument(
            "-m",
            "--float_memory_in_GiB",
            help="(default: 50) the maximum memory usage of the pipeline in GiB",
            default=50,
            type=float,
        )
        arg_grp_general.add_argument(
            "-v", 
            "--verbose", 
            help="turn on verbose mode", 
            action="store_true"
        )
        arg_grp_alignment = parser.add_argument_group("Alignment")
        arg_grp_alignment.add_argument(
            "-Q", 
            "--int_min_mapq", 
            help="(default: 1) minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome", 
            default=1,
            type=int,
        )
        # define adaptor sequences (10X)
        # define cell barcode and umi length
        arg_grp_barcode_extraction = parser.add_argument_group("Barcode Extraction")
        arg_grp_barcode_extraction.add_argument( "-x", "--int_length_cb", help = "(default: 16) the length of the cell barcode", default = 16, type = int )
        arg_grp_barcode_extraction.add_argument( "-y", "--int_length_umi", help = "(default: 12) the length of the UMI (unique molecular identifier)", default = 12, type = int )
        arg_grp_barcode_extraction.add_argument( "-r", "--str_seq_r1", help = "(default: CTACACGACGCTCTTCCGATCT) the sequence of R1 (Read1) adaptor (in 10x GEX v3 kit, located upstream of CB and UMI)", default = 'CTACACGACGCTCTTCCGATCT' )
        arg_grp_barcode_extraction.add_argument( "-e", "--str_seq_tso", help = "(default: AAGCAGTGGTATCAACGCAGAG) the sequence of TSO (Template Switching Oligo) adaptor (in 10x GEX v3 kit, located at 5' end of the molecule)", default = 'AAGCAGTGGTATCAACGCAGAG' )
        arg_grp_barcode_extraction.add_argument( "-E", "--float_error_rate", help = "(default: 0.2) maximum error rate to consider when searching adaptor sequence in the read", default = 0.2, type = float )
        
        arg_grp_cb_correction = parser.add_argument_group("Cell Barcode Correction")
        arg_grp_cb_correction.add_argument( "-V", "--path_file_valid_cb", help = "(required argument) the path to tsv file of whitelist barcodes. For more details, please see 10x cellranger references." ) # required argument
        arg_grp_cb_correction.add_argument( "-N", "--int_max_num_cell_expected", help = "(default: 20000) the max number of expected cells", default = 1000, type = int )

        arg_grp_internal_polyt = parser.add_argument_group("Internal Poly(A) Tract-Primed Read Identification")
        arg_grp_internal_polyt.add_argument( "-S", "--int_len_sliding_window_internal_polyT", help = "(default: 10) the length of sliding window for searching internal poly T (poly A) tract. (When poly-A tailed read is reverse complemented, R1 adaptor become situated in the forward direction", type = int, default = 10 )
        arg_grp_internal_polyt.add_argument( "-w", "--int_len_window_internal_polyT", help = "(default: 30) the size of window for searching for internal poly T", type = int, default = 30 )
        arg_grp_internal_polyt.add_argument( "-F", "--float_min_T_fraction", help = "(default: 0.8) the minimum T fraction for identifying the stretch of poly T tract", type = float, default = 0.8 )
                    
        arg_grp_umi_clustering = parser.add_argument_group("UMI Clustering")
        arg_grp_umi_clustering.add_argument( "-O", "--int_min_n_overlap_kmer_for_clustering_umi", help = "(default: 1) the minimum number of overlapped kmer for initiating UMI clustering ", type = int, default = 1 )
        arg_grp_umi_clustering.add_argument( "-L", "--int_len_kmer_for_clustering_umi", help = "(default: 7) the length of kmer for clustering UMI", type = int, default = 7 )
        arg_grp_umi_clustering.add_argument( "-P", "--float_min_proportion_read_to_select_kmer_representative_for_clustering_umi", help = "(default: 0.75) if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.", type = float, default = 0.75 )
        
        arg_grp_size_distribution = parser.add_argument_group("Collecting Molecule Length Distributions")
        arg_grp_size_distribution.add_argument( "-B", "--int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level", help = "(default: 50) the size of the bin (in base pairs) for collecting size distributions at the single-cell level", type = int, default = 50 )
        
        args = parser.parse_args( )

        l_path_file_bam_input = args.l_path_file_bam_input
        l_path_folder_output = args.l_path_folder_output
        n_threads = args.n_threads
        int_num_base_pairs_in_a_batch = args.int_num_base_pairs_in_a_batch
        int_num_samples_analyzed_concurrently = args.int_num_samples_analyzed_concurrently
        float_memory_in_GiB = args.float_memory_in_GiB
        verbose = args.verbose
        int_min_mapq = args.int_min_mapq
        str_seq_r1 = args.str_seq_r1
        str_seq_tso = args.str_seq_tso
        float_error_rate = args.float_error_rate
        int_length_cb = args.int_length_cb
        int_length_umi = args.int_length_umi
        path_file_valid_cb = args.path_file_valid_cb
        int_max_num_cell_expected = args.int_max_num_cell_expected
        int_len_sliding_window_internal_polyT = args.int_len_sliding_window_internal_polyT
        int_len_window_internal_polyT = args.int_len_window_internal_polyT
        float_min_T_fraction = args.float_min_T_fraction
        int_min_n_overlap_kmer_for_clustering_umi = args.int_min_n_overlap_kmer_for_clustering_umi
        int_len_kmer_for_clustering_umi = args.int_len_kmer_for_clustering_umi
        float_min_proportion_read_to_select_kmer_representative_for_clustering_umi = args.float_min_proportion_read_to_select_kmer_representative_for_clustering_umi
        int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level = args.int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level

    """
    Start of the pipeline
    """
    logger.info(str_description)
    logger.info(
        "Ouro-Tools LongExtractBarcodeFromBAM, a pipeline for preprocessing BAM file for extracting barcode information from user-aligned BAM file using the FASTQ file pre-processed by 'LongFilterNSplit' "
    )
    logger.info(f"Started.")

    """ handle special cases and invalid inputs """
    if l_path_file_bam_input is None : # when inputs are not given
        logger.error(
            "Required argument(s) is missing. to view help message, type -h or --help"
        )
        return -1

    """ process required input directories """

    """ process input directory  """
    l_path_file_bam_input = list(
        os.path.abspath(path_file_bam_input)
        for path_file_bam_input in l_path_file_bam_input
    )
    if l_path_folder_output is not None:
        """# when a valid list of output folders were given # ensure directories of the output folder ends with '/' characters"""
        l_path_folder_output = list(
            os.path.abspath(path_folder) + "/" for path_folder in l_path_folder_output
        )
    else:
        """# compose a list of default 'path_folder_output' values for the given list of input BAM files"""
        l_path_file_bam_input_reversed = deepcopy(
            l_path_file_bam_input[::-1]
        )  # reverse the input file paths so that pop operation yield the element located at the front
        l_path_folder_output = []
        for str_mode_ouro_count in l_str_mode_ouro_count:
            path_file_bam = l_path_file_bam_input_reversed.pop()
            path_folder_output = (
                f"{path_file_bam.rsplit( '/', 1 )[ 0 ]}LongExtractBarcodeFromBAM_output/"
            )
            l_path_folder_output.append(path_folder_output)

    """ 
    Fixed Settings
    """
    # internal settings
    int_highest_mapq = 60
    # define interger representation of the CIGAR operations used in BAM files
    int_cigarop_S = 4
    int_cigarop_H = 5
    # output file setting
    l_col_read_analysis = [  ]
    # calculate padding 
    int_length_cb_umi = int_length_cb + int_length_umi 
    int_length_cb_umi_padding = int( np.ceil( int_length_cb_umi * float_error_rate ) )
    int_length_cb_umi_including_padding = int_length_cb_umi + int_length_cb_umi_padding

    """ 
    Load shared data
    """
    # retrieve set of valid cell barcodes
    set_valid_cb = set( pd.read_csv( path_file_valid_cb, header = None, sep = '\t' ).iloc[ :, 0 ].values )

    """
    Exit early when no samples is anlayzed
    """
    # if no samples will be analyzed, return
    if len(l_path_folder_output) == 0:
        logger.info(f"no output folders were given, exiting")
        return

    """
    Initiate pipelines for off-loading works
    """
    pipelines = bk.Offload_Works(
        None
    )  # no limit for the number of works that can be submitted.

    int_num_samples_analyzed_concurrently = min(
        len(l_path_folder_output), int_num_samples_analyzed_concurrently
    )  # if the number of samples are smaller than 'int_num_samples_analyzed_concurrently', adjust 'int_num_samples_analyzed_concurrently' so that it matches the number of samples

    n_threads = int(
        np.ceil(n_threads / int_num_samples_analyzed_concurrently)
    )  # divide the number of processes that can be used by each pipeline by the number of pipelines that will be run concurrently.

    """
    Pipeline specific functions
    """    
    """
    Ouro-Tools
    BAM tags description
    # 2023-08-13 17:59:06 

    'CB', 'Z' : corrected cell barcode sequence (cell barcode correction at sample level, using all reads from a BAM)
    'UB', 'Z' : corrected UMI sequence using UMI clustering (UMI clustering was performed using all reads with the same each poly(A)-tail attached site, considering the direction of the transcript)
    'UR', 'Z' : uncorrected UMI sequence before UMI clustering (a subsequence of the raw sequence from the 'CU' tag).
    'XR', 'i' : the number of errors for identification of R1 adaptor (marks the 3' end of the original cDNA molecule for 10x GEX 3' products, where cell barcode and UMI sequences can be found). -1 indicates that the adaptor was not identified.
    'XT', 'i' : the number of errors for identification of TSO adaptor (marks the 5' end of the original cDNA molecule). -1 indicates that the adaptor was not identified.
    'CU', 'Z' : the uncorrected raw CB-UMI sequence before processing
    'IA', 'i' : the length of detected internal poly(A) priming region in the genomic alignment. 
    'LE', 'i' : the total length of genomic regions that are actually covered by the read, excluding spliced introns (the sum of exons).
    """
    def _detect_poly_t_length(
        seq_after_softclipping,
        int_len_window_internal_polyT=30,
        int_len_sliding_window_internal_polyT=10,
        float_min_T_fraction=0.8,
    ):
        """ # 2023-08-08 23:22:52 
        detect the length of poly T tract
        """
        ba = bitarray(len(seq_after_softclipping))
        ba.setall(0)

        for index, base in enumerate(seq_after_softclipping):
            ba[index] = base == "T"

        int_len_internal_polyT = 0
        if (
            ba[:int_len_sliding_window_internal_polyT].count()
            / int_len_sliding_window_internal_polyT
            >= float_min_T_fraction
        ):
            int_len_internal_polyT = int_len_sliding_window_internal_polyT
            for index in range(
                1, int_len_window_internal_polyT - int_len_sliding_window_internal_polyT + 1
            ):
                if (
                    ba[index : index + int_len_sliding_window_internal_polyT].count()
                    / int_len_sliding_window_internal_polyT
                    < float_min_T_fraction
                ):
                    break
                int_len_internal_polyT += 1
        return int_len_internal_polyT
    
    l_name_type_dist = [
        'aligned_to_genome', # 0

        'aligned_to_genome__R1__TSO', # 1
        'aligned_to_genome__no_R1__TSO', # 2
        'aligned_to_genome__R1__no_TSO', # 3
        'aligned_to_genome__no_R1__no_TSO', # 4

        'aligned_to_genome__R1__no_valid_CB', # 5
        'aligned_to_genome__R1__valid_CB', # 6
        'aligned_to_genome__R1__valid_CB__internal_polyA', # 7
        'aligned_to_genome__R1__valid_CB__no_internal_polyA', # 8

        'aligned_to_genome__R1__valid_CB__UMI_deduplicated', # 9

        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__1', # 10
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__2to3', # 11
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__4to7', # 12
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__8to15', # 13
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__16to31', # 14
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__32to63', # 15
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__64to127', # 16
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__128to255', # 17
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__256to511', # 18
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__512to1023', # 19
        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__above1024', # 20
    ] # list of distribution types
    def _initialize_dict_arr_dist( ) :
        """ # 2023-08-13 21:32:17 
        initialize 'dict_arr_dist'
        different from, LongFilterNSplit, length of molecule is calculated as the total length of the genomic regions actually covered by the aligned read (the total length of the exons covered by the read.)
        """
        return dict( (e, None) for e in l_name_type_dist )
    
    def run_pipeline():
        """# 2023-10-03 23:36:11 
        analyze a pipeline for a given list of samples
        """
        # retrieve id of the pipeline
        str_uuid_pipeline = bk.UUID()
        logger.info(
            f"[Pipeline Start] Forked Pipeline (id={str_uuid_pipeline}) Started."
        )

        """
        Initiate workers for off-loading works
        """
        workers = bk.Offload_Works(
            None
        )  # no limit for the number of works that can be submitted.

        """
        Run pipeline for each sample
        """
        for path_file_bam_input, path_folder_output in zip( l_path_file_bam_input, l_path_folder_output ) :  # retrieve an output folder for the current sample
            """
            define a function to release a lock
            """
            def release_lock():
                """# 2023-01-14 20:36:17
                release the lock file
                """
                path_file_lock = (
                    f"{path_folder_output}ourotools.lock"
                )

                # check the existence of output files for the output folder of each input file of the current sample
                flag_all_output_files_exist = True  # initialize the flag
                
                if not os.path.exists(
                    f"{path_folder_output}pipeline_completed.txt"
                ):
                    flag_all_output_files_exist = False

                # check the existence of the lock file
                if (
                    os.path.exists(path_file_lock) and flag_all_output_files_exist
                ):  # if all output files exist and the lock file exists
                    # check whether the lock file has been created by the current pipeline
                    with open(path_file_lock, "rt") as file_lock:
                        str_uuid_pipeline_lock = file_lock.read() # retrieve uuid of lock
                        flag_lock_acquired = str_uuid_pipeline_lock == str_uuid_pipeline
                    if (
                        flag_lock_acquired
                    ):  # if the lock file has been created by the current pipeline, delete the lock file
                        os.remove(path_file_lock)
                        # lock has been released
                        if verbose:
                            logger.warning(
                                f"[{path_folder_output}] The forked pipeline (id={str_uuid_pipeline}) released the lock"
                            )
                    else :
                        # lock has been released
                        if verbose:
                            logger.warning(
                                f"[{path_folder_output}] The lock belongs to the forked pipeline (id={str_uuid_pipeline_lock}), and the lock was not released."
                            )
                else:
                    if verbose:
                        logger.warning(
                            f"[{path_folder_output}] The forked pipeline (id={str_uuid_pipeline}) attempted to release the lock, but some output files are missing, and the lock will not be released, yet."
                        )

            """
            Run pipeline for each sample
            """
            """
            create a lock
            """
            os.makedirs(path_folder_output, exist_ok=True)
            path_file_lock = f"{path_folder_output}ourotools.lock"

            # check the existence of the lock file
            if os.path.exists(path_file_lock):
                logger.warning( f"[Output folder unavailable] the output folder {path_folder_output} contains a lock file, which appears to be processed by a different process. Therefore, the output folder will be skipped." )
                continue
            flag_lock_acquired = False  # initialize 'flag_lock_acquired'
            try:
                # create the lock file
                with open(path_file_lock, "wt") as newfile_lock:
                    newfile_lock.write(str_uuid_pipeline)
                # check whether the lock file has been created correctly (check for collision).
                with open(path_file_lock, "rt") as file_lock:
                    flag_lock_acquired = file_lock.read() == str_uuid_pipeline
            except Exception as e:
                logger.critical(
                    e, exc_info=True
                )  # if an exception occurs, print the error message
            if not flag_lock_acquired:
                logger.warning(
                    f"[Output folder unavailable] an attempt to acquire a lock for the output folder {path_folder_output} failed, which appears to be processed by a different process. Therefore, the output folder will be skipped."
                )
                continue
            # lock has been acquired

            """
            Run pipeline for each input file
            """
            # define folders and directories
            path_file_bam_input = os.path.abspath(path_file_bam_input)
            if path_folder_output is None:  # set default 'path_folder_output'
                path_folder_output = (
                    f"{path_file_bam.rsplit( '/', 1 )[ 0 ]}LongExtractBarcodeFromBAM_output/"
                )
            path_folder_output = os.path.abspath(path_folder_output)
            path_folder_output += "/"
            path_folder_temp = f"{path_folder_output}temp/"
            path_folder_graph = f"{path_folder_output}graph/"

            """ if the output folder already exists """
            if os.path.exists(path_folder_output):
                """check whether the pipeline has been completed"""
                if os.path.exists( f"{path_folder_output}pipeline_completed.txt" ) :  # intermediate files should not exists, while all output files should exist
                    logger.info(
                        f"[Output folder Already Exists] the output folder {path_folder_output} contains valid output files. Therefore, the output folder will be skipped."
                    )
                    release_lock( ) # release the lock
                    continue  # skip if the pipeline has been completed for the output folder
                else:
                    """if required output files does not exist or the an intermediate file exists, remove the entire output folder, and rerun the pipeline"""
                    if (
                        len(glob.glob(f"{path_folder_output}*/")) > 0
                    ):  # detect a folder inside the output folder and report the presence of the existing folders.
                        logger.info(
                            f"[Output folder Already Exists] the output folder {path_folder_output} does not contain valid output files. The output folder will be cleaned and the pipeline will start anew at the folder."
                        )
                    # delete the folders
                    for path_folder in glob.glob(f"{path_folder_output}*/"):
                        shutil.rmtree(path_folder, ignore_errors = True)
                    # delete the files, excluding the lock file that has been acquired by the current pipeline
                    for path_file in glob.glob(f"{path_folder_output}*"):
                        if (
                            path_file_lock != path_file
                        ):  # does not delete the lock file
                            os.remove(path_file)

            """ create directories """
            for path_folder in [
                path_folder_output,
                path_folder_temp,
                path_folder_graph,
            ]:
                os.makedirs(path_folder, exist_ok=True)

            """
            Report program arguments
            """
            # record arguments used for the program (metadata)
            dict_program_setting = {
                "version": _version_,  # record version
                # external
                "flag_usage_from_command_line_interface" : flag_usage_from_command_line_interface,
                "path_file_bam_input" : path_file_bam_input,
                "path_folder_output" : path_folder_output,
                "n_threads" : n_threads,
                "int_num_samples_analyzed_concurrently" : int_num_samples_analyzed_concurrently,
                "int_num_base_pairs_in_a_batch" : int_num_base_pairs_in_a_batch,
                "int_min_mapq" : int_min_mapq,
                "float_memory_in_GiB" : float_memory_in_GiB,
                # internal
                "path_folder_temp": path_folder_temp,
                "path_folder_graph": path_folder_graph,
            }
            logger.info(
                f"[Setting] program will be run with the following setting for the input file {path_file_bam_input} : {str( dict_program_setting )}"
            )

            """ export program setting """
            path_file_json_setting_program = (
                f"{path_folder_output}program_setting.json"
            )
            if os.path.exists(path_file_json_setting_program):
                with open(path_file_json_setting_program, "r") as file:
                    j = json.load(file)
                if j != dict_program_setting:
                    logger.info(
                        f"[Warning] the current program setting is different from the previous program setting recorded in the pipeline folder. The previous setting will be used."
                    )
                    with open(path_file_json_setting_program, "r") as file:
                        dict_program_setting = json.load(
                            file
                        )  # override current program setting with previous program setting
            with open(path_file_json_setting_program, "w") as newfile:
                json.dump(dict_program_setting, newfile)
                
            """
            Define a generator for partitioning input file
            """
            def gen_batch( ):
                """# 2023-07-30 18:37:49 
                create batch from the input BAM file
                """
                with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
                    gen_r = samfile.fetch( ) # open the generator
                    ''' retrieve the first valid read '''
                    while True :
                        ''' retrieve a read '''
                        try :
                            r = next( gen_r ) # retrieve the first read
                        except StopIteration : # if the bam file is emtpy, end the generator
                            return
                        
                        """ filter read """
                        if r.mapq < int_min_mapq : # filter out reads with low mapq
                            continue
                        if r.seq is None : # consider only the primary alignment
                            continue
                        len_seq = len( r.seq ) # retrieve the length of the sequence
                        cigartuples, flags = r.cigartuples, r.flag # retrieve attributes
                        if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                            continue 
                        if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                            continue
                        ''' once the first valid read is found, continue to the next step '''
                        break

                    # initialize the batch
                    ns_batch = { 'int_num_base_pairs_encountered_for_a_batch' : len( r.seq ), 'start__reference_name' : r.reference_name, 'start__reference_start' : r.reference_start, } # initialize the dictionary containing information about the batch using the first valid read # counts of base pairs in a batch
                    
                    while True :
                        """ retrieve a read """ 
                        try :
                            r = next( gen_r )
                        except StopIteration : # once all reads were analyzed, exit the loop
                            yield ns_batch # yield the last batch
                            break
                        
                        """ filter read """
                        if r.mapq < int_min_mapq : # filter out reads with low mapq
                            continue
                        if r.seq is None : # consider only the primary alignment
                            continue
                        len_seq = len( r.seq ) # retrieve the length of the sequence
                        cigartuples, flags = r.cigartuples, r.flag # retrieve attributes
                        if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                            continue 
                        if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                            continue
                            
                        """ when contig changes """
                        if r.reference_name != ns_batch[ 'start__reference_name' ] :
                            yield ns_batch # yield the last batch for the last contig
                            # initialize the next batch
                            ns_batch = { 'int_num_base_pairs_encountered_for_a_batch' : 0 } # initialize the counter
                            ns_batch[ 'start__reference_name' ] = r.reference_name
                            ns_batch[ 'start__reference_start' ] = r.reference_start
                            
                        ns_batch[ 'int_num_base_pairs_encountered_for_a_batch' ] += len_seq # increase the base pair count
                        if int_num_base_pairs_in_a_batch <= ns_batch[ 'int_num_base_pairs_encountered_for_a_batch' ] : # once the batch is full, yield the batch and consume remaining reads starting at the reference start position, so that the reads of the same reference start position are processed together. # pipe overloading might happens, causing dead lock. in this case, 'int_num_base_pairs_in_a_batch' can be lowered.
                            # update batch information
                            ns_batch[ 'end__reference_start' ] = r.reference_start
                            while True :
                                """ retrieve a read """ 
                                try :
                                    r = next( gen_r )
                                except StopIteration : # once all reads were analyzed, exit the loop
                                    break
                                    
                                """ filter read """
                                if r.mapq < int_min_mapq : # filter out reads with low mapq
                                    continue
                                if r.seq is None : # consider only the primary alignment
                                    continue
                                len_seq = len( r.seq ) # retrieve the length of the sequence
                                cigartuples, flags = r.cigartuples, r.flag # retrieve attributes
                                if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                                    continue 
                                if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                                    continue
                                    
                                """ check boundary condition """
                                if ns_batch[ 'end__reference_start' ] != r.reference_start : # when the 'reference_start' position changes, finish the batch
                                    break
                                
                                ns_batch[ 'int_num_base_pairs_encountered_for_a_batch' ] += len_seq # increase the counter
                            yield ns_batch # yield the batch
                            ns_batch = { 'int_num_base_pairs_encountered_for_a_batch' : len_seq } # initialize the counter
                            ns_batch[ 'start__reference_name' ] = r.reference_name
                            ns_batch[ 'start__reference_start' ] = r.reference_start
                            
            def process_batch(pipe_receiver, pipe_sender):
                """ # 2023-08-09 00:26:28 
                """
                """
                initialize the worker 
                # 2023-08-01 12:19:06 
                """
                str_uuid = bk.UUID()  # retrieve id
                if verbose:
                    logger.info(f"[Started] start working (worker_id={str_uuid})")
                    
                """ open output files """
                path_file_bam_preprocessed = f"{path_folder_temp}{str_uuid}.preprocessed.bam"
                with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
                    newsamfile = pysam.AlignmentFile( path_file_bam_preprocessed, 'wb', template = samfile ) # open the new samfile, based on the input BAM file
                
                while True:
                    ins = pipe_receiver.recv()
                    if ins is None:
                        break
                    ns_batch = ins  # parse input
                    
                    """
                    define batch-specific function
                    """
                    
                    """
                    open and process the input BAM file
                    """
                    int_total_num_records_processed = 0
                    l_cb_umi = [ ] # collect cb_umi sequences
                    start__reference_name, start__reference_start = ns_batch[ 'start__reference_name' ], ns_batch[ 'start__reference_start' ]
                    end__reference_start = ns_batch[ 'end__reference_start' ] if 'end__reference_start' in ns_batch else None
                    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
                        for r in samfile.fetch( start__reference_name, start__reference_start, end__reference_start + 1 ) if end__reference_start is not None else samfile.fetch( start__reference_name, start__reference_start ) : # include the end position by adding +1 # if 'end__reference_start' is None, retrieve all reads remaining for the contig
                            ''' if the batch has not been started, skip the read '''
                            reference_start = r.reference_start
                            if reference_start < start__reference_start : 
                                continue
                            
                            """ filter read """
                            if r.mapq < int_min_mapq : # filter out reads with low mapq
                                continue
                            seq = r.seq
                            if seq is None : # consider only the primary alignment
                                continue
                            len_seq = len( seq ) # retrieve the length of the sequence
                            cigartuples, flags = r.cigartuples, r.flag # retrieve attributes
                            if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                                continue 
                            if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                                continue
                                
                            ''' if the batch has been completed, exit the loop '''
                            if end__reference_start is not None and reference_start > end__reference_start : 
                                break
                            
                            ''' process read '''
                            
                            """
                            (Assumes the aligned FASTQ files are already pre-processed by ouro-tools and poly A tail is located in the downstream of the read.)
                            
                            not reverse complemented:
                                - poly A and cell barcodes (reverse complemented) located at the right
                            
                            reverse complemented:
                                - poly T and cell barcodes located at the left
                            """
                            int_total_num_records_processed += 1 # update 'int_total_num_records_processed'
                            # check whether the read was reverse complemented
                            flag_is_reverse_complemented = _check_binary_flags( flags, 4 ) 
                            
                            # estimate the size of the molecule (total length of genomic regions covered by the molecule, excluding the soft-clipped regions)
                            l_seg, int_total_length_covering_genome = SAM.Retrive_List_of_Mapped_Segments( cigartuples, reference_start )                          

                            # retrieve soft-clipped sequences
                            flag_left_softclipped = int_cigarop_S == cigartuples[ 0 ][ 0 ]
                            flag_right_softclipped = int_cigarop_S == cigartuples[ -1 ][ 0 ]
                            if not ( flag_left_softclipped and flag_right_softclipped ) : # skip reads that does not contain soft-clipped reads at both ends (adaptors not detected at least one end)
                                continue 
                            int_length_softclipped_left = cigartuples[ 0 ][ 1 ]
                            int_length_softclipped_right = cigartuples[ -1 ][ 1 ]
                            seq_sc_left = seq[ : int_length_softclipped_left ]
                            seq_sc_right = SEQ.Reverse_Complement( seq[ - int_length_softclipped_right : ] )
                            
                            # search for R1 and TSO adaptor sequences
                            seq_sc_with_r1, seq_sc_with_tso, seq_r1_is_located_left, int_length_softclipped_with_r1 = ( seq_sc_left, seq_sc_right, seq, int_length_softclipped_left ) if flag_is_reverse_complemented else ( seq_sc_right, seq_sc_left, SEQ.Reverse_Complement( seq ), int_length_softclipped_right )
                            res_r1 = STR.Search_Subsequence( seq_sc_with_r1, str_seq_r1, float_error_rate )
                            res_tso = STR.Search_Subsequence( seq_sc_with_tso, str_seq_tso, float_error_rate )
                            
                            # initialize the tags that will be added to the SAM record
                            l_tags = [ ( 'XR', res_r1[ 'num_errors' ], 'i' ), ( 'XT', res_tso[ 'num_errors' ], 'i' ), ( 'LE', int_total_length_covering_genome, 'i' ) ] # add the number of errors from R1 and TSO adaptor search results as tags
                            
                            ''' Retrieve Cell Barcode and Check for Internal PolyA Priming (looking for reference-derived polyT next to Cell Barcode in the aligned reads) '''
                            # retrieve cell barcode and UMI
                            int_start_cb_umi = res_r1[ 'index_end_subsequence' ]
                            if int_start_cb_umi != -1 : # if R1 adaptor sequence was identified
                                seq_cb_umi = seq_r1_is_located_left[ int_start_cb_umi : int_start_cb_umi + int_length_cb_umi_including_padding ] # retrieve cb-umi sequence # including sequences that are 'accidently' aligned to the genome
                                # Check for Internal PolyA Priming 
                                seq_after_softclipping = seq_r1_is_located_left[ int_length_softclipped_with_r1 : int_length_softclipped_with_r1 + int_len_window_internal_polyT ]
                                int_length_internal_polyT = _detect_poly_t_length( seq_after_softclipping, int_len_window_internal_polyT, int_len_sliding_window_internal_polyT, float_min_T_fraction )
                                int_count_T_in_a_window = seq_after_softclipping.count( 'T' )
                                # add tags
                                l_tags += [ ("CU", seq_cb_umi, 'Z'), ('IA', int_length_internal_polyT, 'i') ] # add uncorrected cb and umi sequence as a tag # add the identified poly T length as a tag
                                # collect data 
                                l_cb_umi.append( seq_cb_umi ) # collect 'seq_cb_umi'

                            ''' write the SAM record ''' 
                            r.set_tags( l_tags ) # set tags
                            newsamfile.write( r ) # write the record to the output BAM file
                            
                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : int_total_num_records_processed, # record the actual number of records processed for the batch
                        'l_cb_umi' : l_cb_umi,
                        'ns_batch' : ns_batch, # for debugging
                    } )  # report the number of processed records

                """ close output files """
                newsamfile.close( )
                # index the resulting BAM file
                pysam.index( path_file_bam_preprocessed )
                
                """ report the worker has completed all works """
                pipe_sender.send( 'completed' )  
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")

            ns = { 'int_num_read_currently_processed' : 0, 'int_num_records_with_cb_umi' : 0, 'l_cb_umi' : [ ], 'l_l' : [ ] }  # define a namespace # initialize total number of reads processed by the algorithm
            name_file_input = path_file_bam_input.rsplit( '/', 1 )[ 1 ] if '/' in path_file_bam_input else path_file_bam_input # retrieve name of the input file

            def post_process_batch(res):
                # update data using the received result
                ns["int_num_read_currently_processed"] += res[ 'int_total_num_records_for_a_batch' ]
                ns["int_num_records_with_cb_umi"] += len( res["l_cb_umi"] ) # update ns["int_num_records_with_cb_umi"]
                ns["l_cb_umi"] += res["l_cb_umi"]
                logger.info( f"[{path_file_bam_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed. CB/UMI sequence identification rate is {np.round(ns['int_num_records_with_cb_umi'] / ns['int_num_read_currently_processed'], 2 ) if ns['int_num_read_currently_processed'] > 0 else np.nan}" )  # report
                ns[ 'l_l' ].append( [ res[ 'ns_batch' ][ 'int_num_base_pairs_encountered_for_a_batch' ], res[ 'ns_batch' ][ 'start__reference_name' ], res[ 'ns_batch' ][ 'start__reference_start' ], np.nan if 'end__reference_start' not in res[ 'ns_batch' ] else res[ 'ns_batch' ][ 'end__reference_start' ], res[ 'int_total_num_records_for_a_batch' ] ] ) # for debugging
            
            """
            Analyze an input BAM file
            """
            if verbose:
                logger.info( f"[{path_file_bam_input}] the analysis pipeline will be run with {n_threads} number of threads" )
            bk.Multiprocessing_Batch_Generator_and_Workers(
                gen_batch=gen_batch(),
                process_batch=process_batch,
                post_process_batch=post_process_batch,
                int_num_threads=n_threads
                + 2,  # one thread for generating batch, another thread for post-processing of the batch
                flag_wait_for_a_response_from_worker_after_sending_termination_signal = True, # wait until all worker exists before resuming works in the main process
            )
            
            ''' write temporary objects for debugging '''
            df_bookmark = pd.DataFrame( ns[ 'l_l' ], columns = [ 'int_num_base_pairs_encountered_for_a_batch', 'start__reference_name', 'start__reference_start', 'end__reference_start', 'int_total_num_records_for_a_batch' ] )
            df_bookmark.to_csv( f"{path_folder_output}df_bookmark.tsv.gz", index = False, sep = '\t' ) # for debugging
            # bk.PICKLE_Write( f"{path_folder_output}l_cb_umi.pickle", ns["l_cb_umi"] ) # ❤️
            
            """ combine results into a single output BAM file """
            path_file_bam_preprocessed = f"{path_folder_temp}preprocessed.bam"
            l_path_file = glob.glob( f"{path_folder_temp}*.preprocessed.bam" ) # retrieve a list of BAM files to combine
            pysam.merge( '--threads', str( min( n_threads, 10 ) ), '-c', '-p', path_file_bam_preprocessed, * l_path_file ) # merge output BAM files
            for path_file in l_path_file : # delete the temporary files
                os.remove( path_file )
            pysam.index( path_file_bam_preprocessed ) # index the input BAM file

            """ 
            Correct and Assign Cell Barcodes to Each Read
            """
            # internal settings
            n_droplets = 100000 # number of droplets generated in 10X Chromium instruments (number of barcoded beads) - with sufficient margin
            n_minimum_count_cb_before_correction = 3 # threshold for filtering cell barcodes before correction
            ''' retrieve list of potentially valid barcodes '''
            l_cb_umi = ns["l_cb_umi"]
            # retrieve all cb sequence by its normal length, and count each unique cb
            s_cb_count = pd.Series( bk.COUNTER( list( e[ : int_length_cb ] for e in l_cb_umi ) ) )
            # filtering uncorrected cb
            s_cb_count = s_cb_count[ s_cb_count >= n_minimum_count_cb_before_correction  ]
            # drop invalid barcodes
            s_cb_count = bk.Series_Subset( s_cb_count, set_valid_cb )
            # sort valid cb by their counts and only retrieve 'int_max_num_cell_expected' number of cell barcodes
            s_cb_count = s_cb_count.sort_values( ascending = False ).iloc[ : int_max_num_cell_expected ]

            # retrieve a list of valid barcodes
            l_cb_valid = s_cb_count.index.values
            set_cb_valid = set( l_cb_valid ) # retrieve a set of valid cell barcodes

            ''' 1) use pre-computed error-correcting dictionary to identify cell barcodes with a single error (substitution/insertion/deletion) '''
            # build a list of possible errors for each base (except for the first base)
            dict_base_to_l_error = dict( )
            for str_base in 'ATGC' :
                l = [ '' ] # simulate 'deletion' of the the current base
                for str_base_error in 'ATGC' :
                    l.append( str_base + str_base_error ) # simulate 'insertion' of a single base after the current base
                    if str_base != str_base_error : # simulate 'substitution' of the current base
                        l.append( str_base_error )
                dict_base_to_l_error[ str_base ] = l
            # build a list of possible errors for the first base
            dict_first_base_to_l_error = dict( )
            for str_base in 'ATGC' :
                l = [ '' ] # simulate 'deletion' of the the current base
                for str_base_error in 'ATGC' :
                    l.extend( [ str_base + str_base_error, str_base_error + str_base ] ) # simulate 'insertion' of a single base BEFORE AND after the current base (consider a single base insertion before the base because the current read is the first base of the raw cell barcode identified from the read)
                    if str_base != str_base_error : # simulate 'substitution' of the current base
                        l.append( str_base_error )
                dict_first_base_to_l_error[ str_base ] = l
            # retrieve mapping between cb with error to error-free cb
            dict_cb_with_error_to_cb = dict( ( cb_valid, cb_valid ) for cb_valid in l_cb_valid ) # initialize 'dict_cb_with_error_to_cb' using self-mapping. if introducing an error to a valid barcode makes the sequence same as a different valid barcode, the error correction should be avoided.
            for cb in l_cb_valid :
                for pos in range( int_length_cb ) :
                    str_base = cb[ pos ]
                    for error in ( dict_first_base_to_l_error[ str_base ] if pos == 0 else dict_base_to_l_error[ str_base ] ) : # if the current position is the first position, use possible errors of the first position
                        cb_with_error = STR.Replace_a_character_at_an_index( cb, pos, error )
                        if cb_with_error in dict_cb_with_error_to_cb :
                            dict_cb_with_error_to_cb[ cb_with_error ] = None # record collision
                        else :
                            dict_cb_with_error_to_cb[ cb_with_error ] = cb # record error-free cb for each cb with an introduced error
            dict_cb_with_error_to_cb = dict( ( kmer, dict_cb_with_error_to_cb[ kmer ] ) for kmer in dict_cb_with_error_to_cb if dict_cb_with_error_to_cb[ kmer ] is not None )

            ''' 2) Using varying number of kmer to identify cell barcodes with many number of errors '''
            dict_len_kmer_to_kmer_from_cb_to_cb = dict( )
            dict_len_kmer_to_kmer_from_cb_to_cb[ 'l_len_kmer' ] = list( range( int( np.floor( int_length_cb / 2 ) ), int_length_cb, 1 ) ) # length of kmer for cb identification
            for int_length_kmer_for_cb_ident in dict_len_kmer_to_kmer_from_cb_to_cb[ 'l_len_kmer' ] :
                dict_kmer_from_cb_to_cb = dict( )
                for cb in l_cb_valid :
                    for str_kmer in SEQ.Generate_Kmer( cb, int_length_kmer_for_cb_ident ) :
                        if str_kmer in dict_kmer_from_cb_to_cb :
                            dict_kmer_from_cb_to_cb[ str_kmer ] = None # record the occurrence of collision
                        else :
                            dict_kmer_from_cb_to_cb[ str_kmer ] = cb # record the cb from which the kmer was derived
                dict_kmer_from_cb_to_cb = dict( ( kmer, dict_kmer_from_cb_to_cb[ kmer ] ) for kmer in dict_kmer_from_cb_to_cb if dict_kmer_from_cb_to_cb[ kmer ] is not None ) # remove the kmer that are containing collisions
                dict_len_kmer_to_kmer_from_cb_to_cb[ int_length_kmer_for_cb_ident ] = dict_kmer_from_cb_to_cb
            dict_len_kmer_to_kmer_from_cb_to_cb[ 'l_len_kmer' ] = sorted( dict_len_kmer_to_kmer_from_cb_to_cb[ 'l_len_kmer' ] )[ : : -1 ] # sort the list from largest kmer length to the smallest kmer length (for identifing cell barcode from which the current cb would likely to be derived from)
            
            # bk.PICKLE_Write( f"{path_folder_output}dict_cb_with_error_to_cb.pickle", dict_cb_with_error_to_cb ) # ❤️
            # bk.PICKLE_Write( f"{path_folder_output}dict_len_kmer_to_kmer_from_cb_to_cb.pickle", dict_len_kmer_to_kmer_from_cb_to_cb ) # ❤️
            # bk.PICKLE_Write( f"{path_folder_output}set_cb_valid.pickle", set_cb_valid ) # ❤️

            ''' define a function for correcting CB sequences retrieved from reads using the different levels of dictionaries '''
            def _correct_cell_barcode( cb_umi_padded : str ) :
                """ # 2023-08-13 17:38:20 
                correct a single cell barcode using various approaches (direct matching, error-correction using a precomputed dictionary, k-mer based identification)
                """
                cb_corrected = np.nan # initialized to 'cb not found'
                ''' 0) if the cell barcode matches a cell barcode in the set of valid cell barcodes, return the cell barcode as-is '''
                cb_from_read_correct_length = cb_umi_padded[ : int_length_cb ] # cell barcode (uncorrecd) with the length of a valid cell barcode
                if cb_from_read_correct_length in set_cb_valid :
                    return cb_from_read_correct_length
                
                ''' 1) use pre-computed error-correcting dictionary to identify cell barcodes with 1 error '''
                for cb_from_read in [ cb_from_read_correct_length, cb_umi_padded[ : int_length_cb - 1 ], cb_umi_padded[ : int_length_cb + 1 ] ] :
                    if cb_from_read in dict_cb_with_error_to_cb :
                        cb_corrected = dict_cb_with_error_to_cb[ cb_from_read ]
                        break
                if not isinstance( cb_corrected, float ) : # if the correct cell barcode was assigned, return the result
                    return cb_corrected
                ''' 2) Using varying number of kmer to identify cell barcodes with many number of errors '''
                cb_from_read = cb_from_read_correct_length # only consider 'int_length_cb' number of bases
                for len_kmer in dict_len_kmer_to_kmer_from_cb_to_cb[ 'l_len_kmer' ] : # from largest kmer length to the smallest kmer length, identify cell barcode from which the current cb would likely to be derived from.
                    dict_kmer_from_cb_to_cb = dict_len_kmer_to_kmer_from_cb_to_cb[ len_kmer ]
                    for seq_kmer in SEQ.Generate_Kmer( cb_from_read, len_kmer ) :
                        if seq_kmer in dict_kmer_from_cb_to_cb :
                            cb_corrected = dict_kmer_from_cb_to_cb[ seq_kmer ]
                            break
                    if not isinstance( cb_corrected, float ) : # if the cell barcode were identified using the given length of kmer, skip the remaining correction process
                        break
                return cb_corrected
            
            ''' define a function for correcting UMI sequences retrieved from reads using the k-mer-based linear clustering methods '''
            def _cluster_umi( l_umi_for_clustering : list ) :
                """ # 2023-08-14 23:49:02 
                cluster UMI sequences of a single cell barcode, and return the lists of corrected UMI sequences and the number of duplicated UMIs for the corrected UMI sequence.
                """
                ''' handle simple cases '''
                if len( l_umi_for_clustering ) <= 1 : # does not perform UMI clustering when the number of given UMI sequences is <= 1
                    if len( l_umi_for_clustering ) == 1 :
                        return l_umi_for_clustering, [ 1 ]
                    else :
                        return [ ], [ ]
                
                ''' perform UMI clustering '''
                # cluster umi with extracted k-mer
                dict_cluster = SEQ.Cluster_with_Kmer( dict_seq_count = bk.COUNTER( l_umi_for_clustering ), int_min_n_overlap_kmer = int_min_n_overlap_kmer_for_clustering_umi, len_kmer = int_len_kmer_for_clustering_umi, float_min_proportion_read_to_select_kmer_representative = float_min_proportion_read_to_select_kmer_representative_for_clustering_umi )

                ''' retrieve a representative UMI sequence for each cluster based on the frequency, and replace the UMI sequences of the cluster with the representative sequence of the cluster '''
                dict_seq_umi_to_str_name_umi_cluster = dict( ) # retrieve mapping
                dict_seq_umi_to_int_num_duplicated_umis = dict( ) # for retrival of umi duplication levels
                for id_c in dict_cluster : # for each cluster
                    c = dict_cluster[ id_c ] # retrieve the cluster
                    # set the seq_umi with the largest count as a name of current umi cluster
                    str_name_umi_cluster, _ = bk.DICTIONARY_Find_Max( c[ 'seq_count' ] )
                    # assign the current umi cluster name to the umi sequences belonging to the current cluster
                    int_num_duplicated_umis = c["n_seq"] # retrieve the total number of molecules in the cluster
                    for seq_umi in c[ 'seq_count' ] :
                        dict_seq_umi_to_str_name_umi_cluster[ seq_umi ] = str_name_umi_cluster
                        dict_seq_umi_to_int_num_duplicated_umis[ seq_umi ] = int_num_duplicated_umis

                l_umi_corrected = list( dict_seq_umi_to_str_name_umi_cluster[ umi_uncorrected ] for umi_uncorrected in l_umi_for_clustering ) # retrieve list of corrected UMI sequences
                l_num_duplicated_umis = list( dict_seq_umi_to_int_num_duplicated_umis[ umi_uncorrected ] for umi_uncorrected in l_umi_for_clustering ) # retrieve list of number of duplicated UMIs
                return l_umi_corrected, l_num_duplicated_umis
            
            """
            Re-analyze pre-processed BAM files
            """
            def process_batch(pipe_receiver, pipe_sender):
                """ # 2023-08-09 00:26:41 
                # 2022-04-24 01:29:59
                """
                """
                initialize the worker 
                # 2023-08-01 12:19:06 
                """
                str_uuid = bk.UUID()  # retrieve id
                if verbose:
                    logger.info(f"[Started] start working (worker_id={str_uuid})")
                    
                """
                Initiate workers for off-loading works for processing each batch
                """
                int_num_workers_for_bucket_processing = 3 # the number of workers for bucket processing
                workers_for_bucket_processing = bk.Offload_Works( int_num_workers_for_bucket_processing )  # 💛 adjustment of the number of workers might be needed.
                    
                """ open output files """
                path_file_bam_barcoded = f"{path_folder_temp}{str_uuid}.barcoded.bam"
                with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
                    newsamfile = pysam.AlignmentFile( path_file_bam_barcoded, 'wb', template = samfile ) # open the new samfile, based on the input BAM file
                
                while True:
                    ins = pipe_receiver.recv()
                    if ins is None:
                        break
                    name_contig = ins  # parse input
                    
                    """
                    define batch-specific function
                    """
                    
                    """
                    open and process the input BAM file
                    """
                    int_max_bucket_deletion_count_before_reinitialize = 10000 # the max number of bucket deletion count before re-initializing the bucket container (python dictionary, when too many keys are deleted, lead to 'memory leak')
                    int_max_num_records_in_a_batch_of_buckets = 20000 # initialize the total number of records in a batch of buckets
                    int_max_num_batches_in_the_result_container_before_flushing = 1 # the max number of batches whose results can be stored in the container before flushing the result to the storage. if this number is too large, the process will consume too much memory 
                    
                    ns = { 'int_total_num_records_processed' : 0, 'int_bucket_deletion_count' : 0 } # create a namespace for buckets # a counter counting the number of bucket deleted from 'dict_poly_a_site_to_l_l'. if the number exceed
                    ns[ 'dict_poly_a_site_to_l_l' ] = dict( ) # a container to collect reads for alignment end position
                    ns[ 'dict_arr_dist' ] = _initialize_dict_arr_dist( ) # initialize 'dict_arr_dist'
                    ns[ 'dict_arr_dist_single_cell_level' ] = dict( ) 

                    reference_name_current = None # initialize the current contig name
                    reference_start_current = None # initialize the current position
                    
                    # set_e = set( ) # ❤️ # for debugging
                    """
                    l_name_type_dist = [
                        'aligned_to_genome', # 0

                        'aligned_to_genome__R1__TSO', # 1
                        'aligned_to_genome__no_R1__TSO', # 2
                        'aligned_to_genome__R1__no_TSO', # 3
                        'aligned_to_genome__no_R1__no_TSO', # 4

                        'aligned_to_genome__R1__no_valid_CB', # 5
                        'aligned_to_genome__R1__valid_CB', # 6
                        'aligned_to_genome__R1__valid_CB__internal_polyA', # 7
                        'aligned_to_genome__R1__valid_CB__no_internal_polyA', # 8

                        'aligned_to_genome__R1__valid_CB__UMI_deduplicated', # 9

                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__1', # 10
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__2to3', # 11
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__4to7', # 12
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__8to15', # 13
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__16to31', # 14
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__32to63', # 15
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__64to127', # 16
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__128to255', # 17
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__256to511', # 18
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__512to1023', # 19
                        'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__above1024', # 20
                    ] # list of distribution types
                    """
                    """
                    define functions for offloading works for multiprocessing
                    """
                    def _process_buckets( l_l_dict_tags_existing ) :
                        """ # 2023-08-14 15:46:16 
                        process reads in a bucket, and return the reads (for offloading)
                        """
                        # initialize
                        l_cb_combined, l_umi_corrected_combined, l_umi_uncorrected_combined = [ ], [ ], [ ] # lists that will contain values of all buckets
                        # initialize 
                        l_l_len = list( [ ] for _ in range( len( l_name_type_dist ) ) ) # each list represents the list of lengths for each type of distribution in the same order
                        dict_arr_len_single_cell_level = dict( )
                        ''' iterate each bucket '''
                        for l_dict_tags_existing in l_l_dict_tags_existing : # analyze all the reads belonging to each bucket as a group.
                            l_seq_cb_umi = list( dict_tags[ 'CU' ] for dict_tags in l_dict_tags_existing ) # retrieve 'l_seq_cb_umi'
                            """ Correct the cell barcode sequence, and collect the read for UMI clustering """
                            l_cb = [ ]
                            for seq_cb_umi in l_seq_cb_umi : # retrieve cb_umi sequence
                                cb_corrected = _correct_cell_barcode( seq_cb_umi ) # correct cell barcode
                                l_cb.append( cb_corrected )

                            ''' retrieve raw UMI sequence (uncertain boundaries, might contain more or less bases than the actual sequenced UMI sequence) '''
                            l_umi_uncorrected = [ ]
                            for seq_cb_umi, cb_assigned in zip( l_seq_cb_umi, l_cb ) :
                                l_umi_uncorrected.append( np.nan if isinstance( cb_assigned, float ) else seq_cb_umi[ int_length_cb : ] )

                            ''' correct UMI sequences '''
                            if len( l_umi_uncorrected ) <= 1 : # if there is less than 2 umi sequences, clustering is not required.
                                ''' handle simple cases '''
                                l_umi_corrected = l_umi_uncorrected
                                l_num_duplicated_umis = np.ones( len( l_umi_corrected ), dtype = int ) # number of duplication UMIs is 1 in this simple case
                            else :
                                l_umi_uncorrected = np.array( l_umi_uncorrected, dtype = object ) # convert to numpy array
                                l_umi_corrected = np.zeros( len( l_umi_uncorrected ), dtype = object ) # initialize empty arrays
                                l_num_duplicated_umis = np.zeros( len( l_umi_uncorrected ), dtype = object ) 
                                # l_umi = l_umi_uncorrected
                                ''' cluster UMI sequences for each cell barcodes '''
                                dict_index = _index_array( l_cb ) # index cell barcode values
                                for cb in dict_index :
                                    l_index = dict_index[ cb ] # retrieve indices of the entries for the current cell barcode
                                    l_umi_corrected[ l_index ], l_num_duplicated_umis[ l_index ] = _cluster_umi( l_umi_uncorrected[ l_index ] ) # perform UMI clustering

                            ''' classify read and collect molecule size (before deduplication) '''
                            dict_cb_umi_to_max_length = dict( ) # mapping for recording max molecule size for each unique cb-umi pair (since the smaller molecules with a UMI are likely fragments of the larger molecule with the same UMIs)
                            for dict_tags_existing, cb, umi_corrected, int_num_duplicated_umis in zip( l_dict_tags_existing, l_cb, l_umi_corrected, l_num_duplicated_umis ) :
                                int_molecule_size = dict_tags_existing[ 'LE' ] # molecule size excluding adaptors (only genomic regions covered by the read are counted)
                                l_l_len[ 0 ].append( int_molecule_size ) # 'aligned_to_genome', # 0

                                # retrieve flags for classifications
                                flag_R1 = dict_tags_existing[ 'XR' ] > -1
                                flag_TSO = dict_tags_existing[ 'XT' ] > -1
                                if flag_R1 and flag_TSO : 
                                    l_l_len[ 1 ].append( int_molecule_size ) # 'aligned_to_genome__R1__TSO', # 1 
                                elif flag_TSO : 
                                    l_l_len[ 2 ].append( int_molecule_size ) # 'aligned_to_genome__no_R1__TSO', # 2 
                                elif flag_R1 : 
                                    l_l_len[ 3 ].append( int_molecule_size ) # 'aligned_to_genome__R1__no_TSO', # 3 
                                else : 
                                    l_l_len[ 4 ].append( int_molecule_size ) # 'aligned_to_genome__no_R1__no_TSO', # 4 
                                
                                if flag_R1 : # only when R1 adaptor was found
                                    flag_valid_CB = isinstance( cb, str ) # retrieve a flag indicating a valid cell barcode has been detected
                                    
                                    if not flag_valid_CB : 
                                        l_l_len[ 5 ].append( int_molecule_size ) # 'aligned_to_genome__R1__no_valid_CB', # 5
                                    else : # if valid CB was detected
                                        l_l_len[ 6 ].append( int_molecule_size ) # 'aligned_to_genome__R1__valid_CB', # 6
                                        int_length_internal_polyA = dict_tags_existing[ 'IA' ] # retrieve the length of the detected internal poly A
                                        if int_length_internal_polyA > 0 : # when internal poly(A) tract has been detected
                                            l_l_len[ 7 ].append( int_molecule_size ) # 'aligned_to_genome__R1__valid_CB__internal_polyA', # 7
                                        else :
                                            l_l_len[ 8 ].append( int_molecule_size ) # 'aligned_to_genome__R1__valid_CB__no_internal_polyA', # 8
                                        l_l_len[ min( 10 + math.floor( math.log2( int_num_duplicated_umis ) ), 20 ) ].append( int_molecule_size ) # 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__1', # 10 ~ # 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__above1024', # 20
                                        seq_cb_umi = cb + umi_corrected # compose cb-umi combination
                                        # collect max molecule size for each unique cb-umi pair
                                        if seq_cb_umi in dict_cb_umi_to_max_length : # if cb-umi pair already exist
                                            if int_molecule_size > dict_cb_umi_to_max_length[ seq_cb_umi ] : # check whether 'int_molecule_size' is larger than the record max molecule size
                                                dict_cb_umi_to_max_length[ seq_cb_umi ] = int_molecule_size # and update the max molecule size if needed
                                        else :
                                            dict_cb_umi_to_max_length[ seq_cb_umi ] = int_molecule_size # cb-umi pair does not exist, add the current molecule size
                                        
                            ''' classify read and collect molecule size (after deduplication, bulk-level) '''
                            l_l_len[ 9 ].extend( list( dict_cb_umi_to_max_length.values( ) ) ) # 'aligned_to_genome__R1__valid_CB__UMI_deduplicated', # 9
                            
                            ''' classify read and collect molecule size (after deduplication, single-cell level)  '''
                            for cb_umi in dict_cb_umi_to_max_length : # for each cb-umi pair
                                int_molecule_size = dict_cb_umi_to_max_length[ cb_umi ] # retrieve the molecule size
                                cb = cb_umi[ : int_length_cb ] # retrieve cell barcode
                                if cb not in dict_arr_len_single_cell_level : # if the list for the current cell barcode does not exist, initialize the list
                                    dict_arr_len_single_cell_level[ cb ] = [ ]
                                dict_arr_len_single_cell_level[ cb ].append( math.floor( int_molecule_size / int_size_bin_in_base_pairs_for_collecting_size_distributions_at_single_cell_level ) ) # record molecule size at single-cell level
                                
                            ''' collect the CB/UMI correction results '''
                            l_cb_combined.extend( l_cb )
                            l_umi_corrected_combined.extend( l_umi_corrected )
                            l_umi_uncorrected_combined.extend( l_umi_uncorrected )
                        dict_arr_len = dict( (n, l) for n, l in zip( l_name_type_dist, l_l_len ) if len( l ) > 0 ) # compose 'dict_arr_len' # include the list only when list is not empty
                        return l_cb_combined, l_umi_corrected_combined, l_umi_uncorrected_combined, dict_arr_len, dict_arr_len_single_cell_level

                    def _post_process_buckets( res, l_l ) :
                        """ # 2023-08-12 16:11:02 
                        perform post-processing of the bucket analysis result
                        (1) write the result of the processed bucket and (2) update the summary metrics using the analysis result of the bucket.
                        """
                        l_cb, l_umi_corrected, l_umi_uncorrected, dict_arr_len, dict_arr_len_single_cell_level = res # parse the result
                        ns[ 'int_total_num_records_processed' ] += len( l_l ) # update 'int_total_num_records_processed'
                        """ write records """
                        for str_cb, str_umi_corrected, str_umi_uncorrected, t_r_and_tags in zip( l_cb, l_umi_corrected, l_umi_uncorrected, l_l ) : # save the sam records to the file
                            # parse records
                            r, dict_tags_existing = t_r_and_tags
                            if isinstance( str_cb, str ) : # if valid cell barcode was assigned, add tags to the reads
                                l_tags = r.get_tags( with_value_type = True ) # initialize l_tags using existing tags (unfortunately, it is not possible to simply add tags to existing tags)
                                l_tags.extend( [ ( 'CB', str_cb, 'Z' ), ( 'UB', str_umi_corrected, 'Z' ), ( 'UR', str_umi_uncorrected, 'Z' ) ] ) # add new tags to the reads
                                r.set_tags( l_tags )
                            newsamfile.write( r ) 
                        ''' update distributions '''
                        ns[ 'dict_arr_dist_single_cell_level' ] = _batch_update_dictionary_of_size_distributions( dict_arr_dist = ns[ 'dict_arr_dist_single_cell_level' ], dict_l_len = dict_arr_len_single_cell_level ) # update single-cell distributions
                        ns[ 'dict_arr_dist' ] = _batch_update_dictionary_of_size_distributions( dict_arr_dist = ns[ 'dict_arr_dist' ], dict_l_len = dict_arr_len ) # update bulk-level distributions

                    def _write_results_from_offloaded_works( flag_wait_all : bool = False ) :
                        """ # 2023-08-27 18:11:43 
                        flag_wait_all : bool = False # if True, wait until all processes completed their works, if False, write results currently contained in the workers object.
                        """
                        for res in ( workers_for_bucket_processing.wait_all if flag_wait_all else workers_for_bucket_processing.collect_results )( flag_return_results = True ).values( ) : # 'flag_wait_all' is True, wait until all processes completed their works. # wait for all submitted works to be completed, and retrieve results for each work
                            _post_process_buckets( res[ 'result' ], res[ 'associated_data' ] ) # save the sam records to the file
                    
                    def _initialize_a_batch_of_buckets( ) :
                        """ # 2023-08-14 16:24:10 
                        """
                        ns[ 'l_l_dict_tags_existing' ] = [ ] # a container of the list of buckets
                        ns[ 'l_l' ] = [ ] # a container of the data associated with the list of buckets
                        ns[ 'int_total_num_records_in_a_batch_of_buckets' ] = 0 # initialize the total number of records in a batch of buckets

                    def _flush_the_current_batch_of_buckets( ) :
                        """ # 2023-08-27 18:11:36 
                        flush the current batch of buckets
                        
                        flag_wait_all : bool = False # if True, wait until all processes completed their works, if False, write results currently contained in the workers object.
                        """
                        if not workers_for_bucket_processing.is_worker_available : # if all workers are working, wait for a while until all workers are idle 
                            _write_results_from_offloaded_works( flag_wait_all = False ) # flush results from offloaded computations, without waiting all works to be completed 
                            _write_results_from_offloaded_works( flag_wait_all = True ) # flush results from offloaded computations, waiting all works to be completed 
                        elif workers_for_bucket_processing.int_num_completed_results >= int_max_num_batches_in_the_result_container_before_flushing : # if the result container became too large, empty the container
                            _write_results_from_offloaded_works( flag_wait_all = False ) # flush results from offloaded computations, without waiting all works to be completed 
                        workers_for_bucket_processing.submit_work( _process_buckets, args = ( ns[ 'l_l_dict_tags_existing' ], ), associated_data = ns[ 'l_l' ] ) # submit the work for offloading (append the list of pysam objects as the data associated with the work) # flush the current batch of the buckets
                        _initialize_a_batch_of_buckets( ) # initialize the next batch of the buckets

                    def _empty_bucket( t_poly_a_site ) :
                        """ # 2023-08-10 21:21:29 
                        empty bucket for the 't_poly_a_site' by clustering UMI of the reads of the bucket and write the reads to the output BAM file
                        """
                        l_l = ns[ 'dict_poly_a_site_to_l_l' ].pop( t_poly_a_site ) # remove the bucket
                        ns[ 'int_bucket_deletion_count' ] += 1 # increase the counter
                        
                        # if t_poly_a_site in set_e : # ❤️
                        #     logger.warn( f"{t_poly_a_site} already processed!" ) # ❤️
                        # set_e.add( t_poly_a_site ) # ❤️
                        
                        ''' if the number of deletion count exceed the deletion count, re-initialize the bucket container '''
                        if ns[ 'int_bucket_deletion_count' ] > int_max_bucket_deletion_count_before_reinitialize :
                            dict_poly_a_site_to_l_l = dict( ) # create a new dictionary for reinitialization
                            for e in ns[ 'dict_poly_a_site_to_l_l' ] : # copy the dictionary
                                dict_poly_a_site_to_l_l[ e ] = ns[ 'dict_poly_a_site_to_l_l' ][ e ]
                            ns[ 'dict_poly_a_site_to_l_l' ] = dict_poly_a_site_to_l_l # reinitialize the container
                            ns[ 'int_bucket_deletion_count' ] = 0 # reset the counter

                        ''' correct CB / cluster UMI sequences and save results to the file '''
                        ns[ 'l_l_dict_tags_existing' ].append( list( dict_tags_existing for r, dict_tags_existing in l_l ) ) # exclude pysam data objects (cannot be pickled) # add bucket to the batch
                        ns[ 'l_l' ].extend( l_l ) # append associated data
                        ns[ 'int_total_num_records_in_a_batch_of_buckets' ] += len( l_l ) # update the number of reads in a batch
                        if ns[ 'int_total_num_records_in_a_batch_of_buckets' ] > int_max_num_records_in_a_batch_of_buckets : # if the number of records in the batch of buckets exceed the limit, flush the current batch
                            _flush_the_current_batch_of_buckets( ) # flush the current batch of the buckets
                            _initialize_a_batch_of_buckets( ) # initialize the next batch of the buckets
                    
                    _initialize_a_batch_of_buckets( ) # initialize the first batch of the buckets
                    with pysam.AlignmentFile( path_file_bam_preprocessed, 'rb' ) as samfile :
                        for r in samfile.fetch( name_contig ) : # analyze all reads (since the pre-processed BAM file only contains valid reads that are already filtered) for the given chromosome
                            seq, cigartuples, flags, reference_name, reference_start, reference_end = r.seq, r.cigartuples, r.flag, r.reference_name, r.reference_start, r.reference_end # retrieve attributes
                            ''' process reads for each 'bucket' (reads with the same poly A tail attachment sites) '''
                            ''' when the contig has changed, empty all buckets '''
                            if reference_name_current != reference_name :  # retrieve a flag for emptying the buckets (when the contig changes)
                                for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                                    _empty_bucket( t_poly_a_site )
                                reference_name_current = reference_name # update the current contig name
                            
                            ''' when the position has changed, detect buckets that should be emptied '''
                            if reference_start_current != reference_start :
                                ''' determine whether to empty bucket or not, based on the current position on the sorted BAM file '''
                                for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                                    ''' regardlesss of whether poly A site is located at the left or the right side of the read, when the current position passes the poly A site, process the bucket '''
                                    flag_is_reverse_complemented, pos = t_poly_a_site # parse 't_poly_a_site'
                                    if pos < reference_start :
                                        _empty_bucket( t_poly_a_site )
                                reference_start_current = reference_start # update 'reference_start_current'
                                
                            ''' process read '''
                            """
                            (Assumes the aligned FASTQ files are already pre-processed by ouro-tools and poly A tail is located in the downstream of the read.)
                            
                            not reverse complemented:
                                - poly A and cell barcodes (reverse complemented) located at the right
                            
                            reverse complemented:
                                - poly T and cell barcodes located at the left
                            """
                            # check whether the read was reverse complemented
                            flag_is_reverse_complemented = _check_binary_flags( flags, 4 ) 
                            
                            dict_tags_existing = dict( r.get_tags( ) ) # retrieve tags
                            if 'CU' in dict_tags_existing : # if cb_umi sequence is present
                                ''' retrieve poly (A) tail attachment site (specific definition: the alignment 'end' position that are closer to the poly (A) tail) '''
                                t_poly_a_site = ( flag_is_reverse_complemented, ( reference_start if flag_is_reverse_complemented else reference_end ) ) # retrieve a tuple indicating the aligned direction and poly A tail attachment position (alignment end position closer to the identified poly A tail)
                                if t_poly_a_site not in ns[ 'dict_poly_a_site_to_l_l' ] : # initialize 'dict_poly_a_site_to_l_l' for 't_poly_a_site'
                                    ns[ 'dict_poly_a_site_to_l_l' ][ t_poly_a_site ] = [ ]
                                ns[ 'dict_poly_a_site_to_l_l' ][ t_poly_a_site ].append( [ r, dict_tags_existing ] )
                            else :
                                ''' write the SAM record (record that does not contain the cell barcode - UMI sequence) ''' 
                                ns[ 'int_total_num_records_processed' ] += 1
                                newsamfile.write( r ) # write the record to the output BAM file
                                # update relevant distributions 
                                ns[ 'dict_arr_dist' ][ 'aligned_to_genome' ] = _update_size_distribution( new_size = dict_tags_existing[ 'LE' ], arr_dist = ns[ 'dict_arr_dist' ][ 'aligned_to_genome' ] ) 
                                name_type_dist = 'aligned_to_genome__no_R1__no_TSO' if dict_tags_existing[ 'XT' ] == -1 else 'aligned_to_genome__no_R1__TSO' # classify read based on the TSO search result
                                ns[ 'dict_arr_dist' ][ name_type_dist ] = _update_size_distribution( new_size = dict_tags_existing[ 'LE' ], arr_dist = ns[ 'dict_arr_dist' ][ name_type_dist ] ) # update the size distribution associated with the read
                            
                        ''' when all reads of the contig were read, empty all buckets and flush the batch, and wait until all computation has been completed '''
                        for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                            _empty_bucket( t_poly_a_site )
                        _flush_the_current_batch_of_buckets( ) # flush the last batch of the buckets
                        _write_results_from_offloaded_works( flag_wait_all = True ) # wait for all works to be completed, and flush results from offloaded computations

                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : ns[ 'int_total_num_records_processed' ], # record the actual number of records processed for the batch
                        'dict_arr_dist' : ns[ 'dict_arr_dist' ], # return results
                        'dict_arr_len_single_cell_level' : ns[ 'dict_arr_dist_single_cell_level' ], 
                    } )  # report the number of processed records

                """ close output files """
                newsamfile.close( )
                # sort the output sam file
                path_file_bam_barcoded_sorted = f"{path_folder_temp}{str_uuid}.barcoded.sorted.bam"
                pysam.sort( "-o", path_file_bam_barcoded_sorted, '-@', str( min( n_threads, 5 ) ), path_file_bam_barcoded )
                os.remove( path_file_bam_barcoded ) # remove the temporary file
                pysam.index( path_file_bam_barcoded_sorted ) # index the resulting BAM file
                
                """ report the worker has completed all works """
                pipe_sender.send( 'completed' )  
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")

            ''' initialize a data structure that will be analyzed in bulk level '''
            ns = { 'int_num_read_currently_processed' : 0 }  # define a namespace for combining results
            ns[ 'dict_arr_dist' ] = _initialize_dict_arr_dist( ) # initialize 'dict_arr_dist' (bulk level, as all reads samples belong to the current BAM file are analyzed.)
            ns[ 'dict_arr_dist_single_cell_level' ] = dict( ) # initialize 'dict_arr_dist_single_cell_level', a summarized (binned) size distribution
            name_type_dist_to_collect_single_cell_level = 'aligned_to_genome__R1__valid_CB__UMI_deduplicated' # the name of the type of the distribution to collect at the single-cell level
            def post_process_batch(res):
                # update data using the received result
                ns["int_num_read_currently_processed"] += res[ 'int_total_num_records_for_a_batch' ]
                logger.info( f"[{path_file_bam_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed." ) # report
                
                # combine distributions (bulk)
                ns[ 'dict_arr_dist' ] = _combine_dictionary_of_size_distributions( dict_arr_dist_existing = ns[ 'dict_arr_dist' ], dict_arr_dist_new = res[ 'dict_arr_dist' ] ) # combine and update the global distributions (bulk-level)
                # combine distributions (single-cell)
                ns[ 'dict_arr_dist_single_cell_level' ] = _combine_dictionary_of_size_distributions( dict_arr_dist_existing = ns[ 'dict_arr_dist_single_cell_level' ], dict_arr_dist_new = res[ 'dict_arr_len_single_cell_level' ] ) # update the container
                    
                logger.info( f"[{path_file_bam_input}] total {np.sum(ns[ 'dict_arr_dist' ][ 'aligned_to_genome' ])} number of reads has been processed. (recalculation using 'dict_arr_dist')" ) # report
                    
            """
            Analyze an input BAM file
            """
            if verbose:
                logger.info( f"[{path_file_bam_input}] the analysis pipeline will be run with {n_threads} number of threads" )
            bk.Multiprocessing_Batch_Generator_and_Workers(
                gen_batch = iter( SAM.Get_contig_names_from_bam_header( path_file_bam_preprocessed ) ), # analyze the pre-processed BAM file for each chromosome
                process_batch=process_batch,
                post_process_batch=post_process_batch,
                int_num_threads=n_threads
                + 2,  # one thread for generating batch, another thread for post-processing of the batch
                flag_wait_for_a_response_from_worker_after_sending_termination_signal = True, # wait until all worker exists before resuming works in the main process
            )
            os.remove( path_file_bam_preprocessed ) # delete the temporary file

            """ 
            post-processing
            """
            def post_processing():  # off-loading a single-core work
                logger.info(
                    f"[{path_file_bam_input}] post-processing started"
                )
                
                # combine results into a single output file (initial read analysis)
                """ combine results into a single output BAM file """
                l_path_file = glob.glob( f"{path_folder_temp}*.barcoded.sorted.bam" ) # retrieve a list of BAM files to combine
                pysam.merge( '--threads', str( min( n_threads, 10 ) ), '-c', '-p', f"{path_folder_output}barcoded.bam", * glob.glob( f"{path_folder_temp}*.barcoded.sorted.bam" ) ) # merge output BAM files
                for path_file in l_path_file : # delete the temporary files
                    os.remove( path_file )
                pysam.index( f"{path_folder_output}barcoded.bam" ) # index the input BAM file

                ''' summarize distributions '''
                dict_arr_dist = ns[ 'dict_arr_dist' ] # retrieve 'dict_arr_dist'
                l_l = [ ]
                for e in dict_arr_dist :
                    arr = dict_arr_dist[ e ]
                    if arr is None :
                        continue
                    arr_bins = np.arange( len( arr ) ) # retrieve bin size of the histograms
                    int_num_reads, int_num_base_pairs = arr.sum( ), ( arr * arr_bins ).sum( )
                    int_avg_length_base_pairs = int_num_base_pairs / int_num_reads
                    float_standard_deviation_length_base_pairs = np.sqrt( np.average((arr_bins - int_avg_length_base_pairs)**2, weights=arr) )
                    l_l.append( [ e, int_num_reads, int_num_base_pairs, int_avg_length_base_pairs, float_standard_deviation_length_base_pairs ] )
                df_summary_of_distributions = pd.DataFrame( l_l, columns = [ 'name_type_distribution', 'int_num_reads', 'int_num_base_pairs', 'int_avg_length_base_pairs', 'float_standard_deviation_length_base_pairs' ] )
                df_summary_of_distributions.to_csv( f"{path_folder_output}df_summary_of_distributions.tsv.gz", sep = '\t', index = False ) # export 'df_summary_of_distributions'
                
                """
                Draw plots of distributions
                """
                # create output folders
                path_folder_graph_noninteractive, path_folder_graph_interactive = f"{path_folder_graph}noninteractive_graph/", f"{path_folder_graph}interactive_graph/"
                for path_folder in [ path_folder_graph_noninteractive, path_folder_graph_interactive ] :
                    os.makedirs( path_folder, exist_ok = True )

                ''' draw simple line plots '''
                # plot settings
                int_max_molecule_size_plot = 6500
                for name_cat_dist in _initialize_dict_arr_dist( ) : # for each category
                    if dict_arr_dist[ name_cat_dist ] is not None :
                        len_max_molecule_size_data = len( dict_arr_dist[ name_cat_dist ] ) # retrieve max molecule size 
                        plt.plot( np.arange( min( int_max_molecule_size_plot, len_max_molecule_size_data ) ), dict_arr_dist[ name_cat_dist ] if len_max_molecule_size_data <= int_max_molecule_size_plot else dict_arr_dist[ name_cat_dist ][ : int_max_molecule_size_plot ] )
                        plt.title( f"{name_cat_dist} ({dict_arr_dist[ name_cat_dist ].sum( )} molecules)" )
                        bk.MPL_SAVE( f"{name_cat_dist}.distribution", folder = path_folder_graph_noninteractive, l_format=['.pdf', '.png'] )
                        
                ''' draw interactive stacked bar graphs '''
                df_bar = _get_df_bar( dict_arr_dist, int_size_bin_in_base_pairs = 50, int_max_size_in_base_pairs = int_max_molecule_size_plot ) # retrieve a dataframe for drawing a bar graph
                for flag_use_proportion in [ True, False ] :
                    _draw_bar_plot( 
                        df_bar, 
                        [ 'aligned_to_genome__R1__TSO', 'aligned_to_genome__no_R1__TSO',  'aligned_to_genome__R1__no_TSO', 'aligned_to_genome__no_R1__no_TSO', ],
                        title = f"R1 and TSO Adaptor Identification of '{name_file_input}'",
                        flag_use_proportion = flag_use_proportion,
                        flag_save_figure = True, path_folder_graph = path_folder_graph_interactive,
                    )
                    _draw_bar_plot( 
                        df_bar, 
                        [ 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__1', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__2to3', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__4to7', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__8to15', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__16to31', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__32to63', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__64to127', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__128to255', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__256to511', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__512to1023', 'aligned_to_genome__R1__valid_CB__UMI_duplication_rate__above1024' ],
                        title = f"UMI Duplication Counts of '{name_file_input}'",
                        flag_use_proportion = flag_use_proportion,
                        flag_save_figure = True, path_folder_graph = path_folder_graph_interactive,
                    )
                    _draw_bar_plot( 
                        df_bar, 
                        ['aligned_to_genome__R1__valid_CB__no_internal_polyA', 'aligned_to_genome__R1__valid_CB__internal_polyA', 'aligned_to_genome__R1__no_valid_CB', 'aligned_to_genome__no_R1__TSO', 'aligned_to_genome__no_R1__no_TSO'],
                        title = f"Internal poly(A) Detection of '{name_file_input}'",
                        flag_use_proportion = flag_use_proportion,
                        flag_save_figure = True, path_folder_graph = path_folder_graph_interactive,
                    )
                
                ''' export pickle files '''
                # write distribution data as pickle files
                bk.PICKLE_Write( f"{path_folder_output}dict_arr_dist.pkl", dict_arr_dist )
                bk.PICKLE_Write( f"{path_folder_output}dict_arr_dist_single_cell_level.pkl", ns[ 'dict_arr_dist_single_cell_level' ] )
                
                # write a flag indicating that the processing has been completed
                with open( f"{path_folder_output}pipeline_completed.txt", 'w' ) as newfile :
                    newfile.write( 'completed' )

                # delete temporary files
                shutil.rmtree( path_folder_temp, ignore_errors = True )
                    
                release_lock()  # release the lock
                logger.info(
                    f"[{path_file_bam_input}] post-processing completed"
                )

            workers.submit_work(post_processing)

            release_lock()  # release the lock

        # wait all the single-core works offloaded to the workers to be completed.
        workers.wait_all()
        logger.info(
            f"[Pipeline Completion] Forked Pipeline (id={str_uuid_pipeline}) Completed."
        )

    for _ in range(
        int_num_samples_analyzed_concurrently
    ):  # run 'int_num_samples_analyzed_concurrently' number of pipelines
        pipelines.submit_work(run_pipeline)

    # wait all pipelines to be completed
    pipelines.wait_all()
    logger.info(f"Completed.")
    return 

def LongCreateReferenceSizeDistribution(
    flag_usage_from_command_line_interface: bool = False,
    l_path_file_distributions: Union[List[str], None] = None, # list of path to the 'dict_arr_dist.pkl' output file of the 'LongExtractBarcodeFromBAM' pipeline for each sample
    l_name_file_distributions : Union[ None, List[ str ] ] = None, # list of the name representing each 'dict_arr_dist.pkl' output file. Should be unique and non-redundant. if None is given, the absolute, real (soft-link resolved) path of the pickle file will be be used as the name representing the file
    path_folder_output: Union[ str, None ] = None, # path to the output folder of the 'LongCreateReferenceSizeDistribution' pipeline
    # smoothening process
    float_sigma_gaussian_filter : float = 8.0, # the standard deviation of the Gaussian filter that will be applied to the base-pair resolution distribution data (histogram with a bin size of 1bp)
    # peak removal
    int_min_total_read_count_for_a_peak : int = 50, # the minimum number of reads in a peak to be considered as a valid peak
    float_min_ratio_read_count_peak_to_baseline : float = 0.05, # the minimum ratio of the number of reads in the peak to the number of reads included in the baseline
    float_min_ratio_peak_height_to_baseline_height : float = 0.5, # the minimum ratio of the peak height to the height of the baseline
    int_size_window_surveying_surrounding_values : int = 4, # the size of the window for estimation of the height of the peak base (since the peak is identified at 50% of its height)
    int_num_iterative_peak_removal : int = 3, # the number of iterative peak removal process for each distribution
    # correction ratio calculation & confidence estimation
    float_max_correction_ratio : float = 10.0, # the maximum correction ratio allowed for estimating the confident molecule size ranges for count matrix normalization, where the correction ratio used for adjusting the distribution of the sample to that of the reference distribution is always below the given threshold, 'float_max_correction_ratio'
    t_distribution_range_of_interest : List[ int ] = [ 1000, 3500 ], # define a range of distribution of interest for searching optimal coefficient for calculating normalization ratios
    # generic    
    float_memory_in_GiB: float = 50,
    verbose: bool = True,
) -> None :
    """# 2023-08-23 21:27:06 
    l_path_file_distributions: Union[List[str], None] = None, # list of path to the 'dict_arr_dist.pkl' output file of the 'LongExtractBarcodeFromBAM' pipeline for each sample
    l_name_file_distributions : Union[ None, List[ str ] ] = None, # list of the name representing each 'dict_arr_dist.pkl' output file. Should be unique and non-redundant. if None is given, the absolute, real (soft-link resolved) path of the pickle file will be be used as the name representing the file
    path_folder_output: Union[ str, None ] = None, # path to the output folder of the 'LongCreateReferenceSizeDistribution' pipeline
    # smoothening process
    float_sigma_gaussian_filter : float = 8.0, # the standard deviation of the Gaussian filter that will be applied to the base-pair resolution distribution data (histogram with a bin size of 1bp)
    # peak removal
    int_min_total_read_count_for_a_peak : int = 50, # the minimum number of reads in a peak to be considered as a valid peak
    float_min_ratio_read_count_peak_to_baseline : float = 0.05, # the minimum ratio of the number of reads in the peak to the number of reads included in the baseline
    float_min_ratio_peak_height_to_baseline_height : float = 0.5, # the minimum ratio of the peak height to the height of the baseline
    int_size_window_surveying_surrounding_values : int = 4, # the size of the window for estimation of the height of the peak base (since the peak is identified at 50% of its height)
    int_num_iterative_peak_removal : int = 3, # the number of iterative peak removal process for each distribution
    # correction ratio calculation & confidence estimation
    float_max_correction_ratio : float = 10.0, # the maximum correction ratio allowed for estimating the confident molecule size ranges for count matrix normalization, where the correction ratio used for adjusting the distribution of the sample to that of the reference distribution is always below the given threshold, 'float_max_correction_ratio'
    t_distribution_range_of_interest : List[ int ] = [ 1000, 3500 ], # define a range of distribution of interest for searching optimal coefficient for calculating normalization ratios
    # generic    
    float_memory_in_GiB: float = 50,
    verbose: bool = True,
    
    returns
    """
    """
    Parse arguments
    """
    if flag_usage_from_command_line_interface:  # parse arguments
        """parse arguments when the function was called from the command-line interface"""
        # {  } # unused arguments
        # command line arguments
        parser = argparse.ArgumentParser(
            description=str_description,
            usage="ourotools LongCreateReferenceSizeDistribution",
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument("LongCreateReferenceSizeDistribution")
        
        arg_grp_general = parser.add_argument_group("General")
        arg_grp_general.add_argument(
            "-i",
            "--l_path_file_distributions",
            help="list of path to the 'dict_arr_dist.pkl' output file of the 'LongExtractBarcodeFromBAM' pipeline for each sample",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-n",
            "--l_name_file_distributions",
            help="list of the name representing each 'dict_arr_dist.pkl' output file. Should be unique and non-redundant. if None is given, the absolute, real (soft-link resolved) path of the pickle file will be be used as the name representing the file",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-o",
            "--path_folder_output",
            help="path to the output folder of the 'LongCreateReferenceSizeDistribution' pipeline",
        )
        arg_grp_general.add_argument(
            "-m",
            "--float_memory_in_GiB",
            help="(default: 50) the maximum memory usage of the pipeline in GiB",
            default=50,
            type=float,
        )
        arg_grp_general.add_argument(
            "-v", 
            "--verbose", 
            help="turn on verbose mode", 
            action="store_true"
        )
        
        arg_grp_gaussian_filter = parser.add_argument_group("Gaussian Filter")
        arg_grp_gaussian_filter.add_argument(
            "-S",
            "--float_sigma_gaussian_filter",
            help="(default: 8.0) the standard deviation of the Gaussian filter that will be applied to the base-pair resolution distribution data (histogram with a bin size of 1bp).",
            default=8.0,
            type=float,
        )
        
        arg_grp_peak_removal = parser.add_argument_group("Peak Removal")
        arg_grp_peak_removal.add_argument(
            "-P",
            "--int_min_total_read_count_for_a_peak",
            help="(default: 50) the minimum number of reads in a peak to be considered as a valid peak.",
            default=50,
            type=int,
        )
        arg_grp_peak_removal.add_argument(
            "-R",
            "--float_min_ratio_read_count_peak_to_baseline",
            help="(default: 0.05) the minimum ratio of the number of reads in the peak to the number of reads included in the baseline.",
            default=0.05,
            type=float,
        )
        arg_grp_peak_removal.add_argument(
            "-H",
            "--float_min_ratio_peak_height_to_baseline_height",
            help="(default: 0.5) the minimum ratio of the peak height to the height of the baseline.",
            default=0.5,
            type=float,
        )
        arg_grp_peak_removal.add_argument(
            "-W",
            "--int_size_window_surveying_surrounding_values",
            help="(default: 4) the size of the window for estimation of the height of the peak base (since the peak is identified at 50% of its height).",
            default=4,
            type=int,
        )
        arg_grp_peak_removal.add_argument(
            "-I",
            "--int_num_iterative_peak_removal",
            help="(default: 3) the number of iterative peak removal process for each distribution.",
            default=3,
            type=int,
        )

        arg_grp_correction = parser.add_argument_group("Distribution Correction") # correction ratio calculation & confidence estimation
        arg_grp_correction.add_argument(
            "-C",
            "--float_max_correction_ratio",
            help="(default: 10.0) the maximum correction ratio allowed for estimating the confident molecule size ranges for count matrix normalization, where the correction ratio used for adjusting the distribution of the sample to that of the reference distribution is always below the given threshold, 'float_max_correction_ratio'.",
            default=10.0,
            type=float,
        )
        arg_grp_correction.add_argument(
            "-t",
            "--t_distribution_range_of_interest",
            help="(default: [ 1000, 3500 ]) define a range of distribution of interest for searching optimal coefficient for calculating normalization ratios.",
            default=[ 1000, 3500 ],
            nargs="*",
        )

        args = parser.parse_args()

        l_path_file_distributions = args.l_path_file_distributions
        l_name_file_distributions = args.l_name_file_distributions
        path_folder_output = args.path_folder_output
        float_memory_in_GiB = args.float_memory_in_GiB
        verbose = args.verbose
        float_sigma_gaussian_filter = args.float_sigma_gaussian_filter
        int_min_total_read_count_for_a_peak = args.int_min_total_read_count_for_a_peak
        float_min_ratio_read_count_peak_to_baseline = args.float_min_ratio_read_count_peak_to_baseline
        float_min_ratio_peak_height_to_baseline_height = args.float_min_ratio_peak_height_to_baseline_height
        int_size_window_surveying_surrounding_values = args.int_size_window_surveying_surrounding_values
        int_num_iterative_peak_removal = args.int_num_iterative_peak_removal
        float_max_correction_ratio = args.float_max_correction_ratio
        t_distribution_range_of_interest = args.t_distribution_range_of_interest
    
    """
    Start of the pipeline
    """
    logger.info(str_description)
    logger.info(
        "Ouro-Tools LongCreateReferenceSizeDistribution, a pipeline for refining and combining size distributions of multiple samples calculated using 'LongExtractBarcodeFromBAM' to create a reference size distribution, and finding optimal correction ratios for correction for count normalization for each sample."
    )
    logger.info(f"Started.")

    """ handle special cases and invalid inputs """
    if ( l_path_file_distributions is None ) or ( path_folder_output is None ) : # check whether the required input paths were given
        logger.error(
            "Required argument(s) is missing. to view help message, type -h or --help"
        )
        return -1

    """ process input directory  """
    l_path_file_distributions = list(
        os.path.abspath(path_file_distributions)
        for path_file_distributions in l_path_file_distributions
    )
    """# when a valid list of output folders were given # ensure directories of the output folder ends with '/' characters"""
    path_folder_output = os.path.abspath(path_folder_output) + "/"

    """ initialize 'l_name_file_distributions' when not given """
    if l_name_file_distributions is None : # if valid 'l_name_file_distributions' is not given, initialize using real, absolute path
        l_name_file_distributions = list( os.path.abspath( os.path.realpath( e ) ) for e in l_path_file_distributions ) # initialize
            
    """ 
    Fixed Settings
    """
    # internal settings
    int_highest_mapq = 60

    """
    Exit early when no samples is anlayzed
    """
    # if no samples will be analyzed, return
    if len(l_path_file_distributions) == 0:
        logger.info(f"no output folders were given, exiting")
        return
        
    """
    Pipeline specific functions
    """

    """
    run pipeline
    """
    
    """
    internal setting
    """
    from scipy.ndimage import gaussian_filter
    from scipy.signal import find_peaks
    from scipy import optimize
    import plotly.express as px
    
    logger.setLevel( logging.INFO ) # reset logging info after importing

    name_type_dist_for_creating_reference = 'aligned_to_genome__R1__valid_CB__UMI_deduplicated' # define name of the type of the distribution for creating the reference
    dict_kw_gaussian_filter = {
        'sigma' : float_sigma_gaussian_filter, 
        'output' : float, 
        'mode' : 'nearest', 
        'truncate' : 4
    }
    dict_kw_find_peaks = {
        'rel_height' : 0.5, 
        'width' : ( 0.5, 15 ), 
        'prominence' : 1, 
        'threshold' : 1
    }

    ''' read distributions '''
    logger.info(f"reading distributions")
    l_arr_dist = list( bk.PICKLE_Read( path_file_distributions )[ name_type_dist_for_creating_reference ] for path_file_distributions in l_path_file_distributions ) # retrieve distributions

    ''' perform peak removal '''
    logger.info(f"performing peak removal")
    def _remove_peaks( 
        arr,
        int_min_total_read_count_for_a_peak : int = int_min_total_read_count_for_a_peak, 
        float_min_ratio_read_count_peak_to_baseline : float = float_min_ratio_read_count_peak_to_baseline, 
        float_min_ratio_peak_height_to_baseline_height : float = float_min_ratio_peak_height_to_baseline_height, 
        int_size_window_surveying_surrounding_values : int = int_size_window_surveying_surrounding_values, 
        dict_kw_find_peaks : dict = dict_kw_find_peaks,
        flag_plot_graph : bool = False,
        figsize : tuple = ( 30, 5 ),
    ) :
        """ # 2023-08-22 22:16:36 
        remove peaks from the given distributions
        """
        ''' search peaks '''
        arr_peaks, properties = find_peaks( arr, ** dict_kw_find_peaks ) 

        ''' search for significant peaks '''
        l_idx_of_significant_peak = [ ]
        for idx_peak in range( len( arr_peaks ) ) : # iterate each peak
            pos_peak_start = math.floor( properties["left_ips"][idx_peak] )
            pos_peak_end = math.ceil( properties["right_ips"][idx_peak] )
            int_height_peak = properties["prominences"][idx_peak]
            int_height_baseline_and_peak = arr[ arr_peaks[ idx_peak ] ]
            int_height_baseline = int_height_baseline_and_peak - int_height_peak # calculate baseline height
            int_total_count_excluding_peak = ( ( arr[ pos_peak_start ] + arr[ pos_peak_end ] ) / 2 ) * ( pos_peak_end - pos_peak_start + 1 )
            int_total_count_including_peak = np.sum( arr[ pos_peak_start : pos_peak_end + 1 ] )
            int_total_count_peak = int_total_count_including_peak - int_total_count_excluding_peak # retrieve total count of peak

            if ( int_total_count_peak >= int_min_total_read_count_for_a_peak ) and ( ( int_total_count_peak / int_total_count_excluding_peak ) >= float_min_ratio_read_count_peak_to_baseline ) and ( ( int_height_baseline <= 0 ) or ( ( int_height_peak / int_height_baseline ) > float_min_ratio_peak_height_to_baseline_height ) ) : # identify significant peak
                l_idx_of_significant_peak.append( idx_peak ) # collect significant peak

        ''' filter out insignificant peaks '''
        # filter and retain only the significant peaks
        arr_peaks = arr_peaks[ l_idx_of_significant_peak ]
        for k in list( properties ) :
            properties[ k ] = properties[ k ][ l_idx_of_significant_peak ]

        ''' remove peaks from the distribution '''
        arr_without_peak = arr.copy( ) # initialize 'arr_without_peak'

        def __log_avg( a ) :
            return np.exp( np.log( a ).mean( ) )
        for idx_peak in range( len( arr_peaks ) ) : # iterate each peak
            pos_peak_start, pos_peak_end = max( 0, math.floor( properties["left_ips"][idx_peak] - int_size_window_surveying_surrounding_values / 2 ) ), min( len( arr ) - 1, math.ceil( properties["right_ips"][idx_peak] + int_size_window_surveying_surrounding_values / 2 ) ) # shift peak start and end positions by 'int_size_window_surveying_surrounding_values'
            val_start, val_end = __log_avg( arr[ pos_peak_start - int_size_window_surveying_surrounding_values : pos_peak_start + 1 ] ), __log_avg( arr[ pos_peak_end : pos_peak_end + int_size_window_surveying_surrounding_values + 1 ] ) # retrieve the log-average of read count values around the peak start and end positions
            slope = ( val_end - val_start ) / ( pos_peak_end - pos_peak_start ) # retrieve the slope of the graph after removing the peak # linear interpolation
            for pos in range( pos_peak_start, pos_peak_end + 1 ) : # remove the peak for the positions covered by the peak
                arr_without_peak[ pos ] = val_start + slope * ( pos - pos_peak_start )

        ''' plot '''
        if flag_plot_graph :
            fig, ax = plt.subplots( 1, 1, figsize = figsize )
            ax.plot( arr, color = 'b', alpha = 0.2 ) # plot original distribution
            ax.plot( arr_without_peak, color = 'b', alpha = 0.5 ) # plot distribution without peak
            arr_dist_smoothened = gaussian_filter( arr_without_peak, sigma = 5, output = float, mode = 'nearest', truncate = 4 ) # smoothen the distribution
            ax.plot( arr_dist_smoothened, '-', lw = 6.5, color = 'g' ) # plot smoothened graph
            # annotate peaks
            ax.plot( arr_peaks, arr[ arr_peaks ], 'x', color = 'C1' ) 
            ax.vlines(x=arr_peaks, ymin=arr[arr_peaks] - properties["prominences"], ymax = arr[arr_peaks], color = "C1", alpha = 0.4)
            ax.hlines(y=properties["width_heights"], xmin=properties["left_ips"], xmax=properties["right_ips"], color = "C1", alpha = 0.4)
        return arr_without_peak # return the result
    def _iterative_peak_removal( arr, int_num_iterative_peak_removal : int = 3 ) :
        """ # 2023-08-22 22:37:33 
        perform iterative peak removal
        """
        arr_without_peak = arr.copy( ) # initialize 'arr_without_peak'
        for _ in range( int_num_iterative_peak_removal ) :
            arr_without_peak = _remove_peaks( arr_without_peak ) # remove peaks
        return arr_without_peak # return the result
    l_arr_dist_peak_removed = list( _iterative_peak_removal( arr, int_num_iterative_peak_removal = int_num_iterative_peak_removal ) for arr in l_arr_dist )

    ''' smoothen the distributions '''
    logger.info(f"smoothening the distributions")
    l_arr_dist_smoothened = list( gaussian_filter( arr, ** dict_kw_gaussian_filter ) for arr in l_arr_dist_peak_removed ) # smoothen the distributions

    ''' calculate the log average of the distributions (reference distribution) '''
    logger.info(f"creating the reference distribution")
    def _log_average_of_distributions( l_arr_dist : List ) :
        """ # 2023-08-22 22:51:51 
        calculate log average of distributions + 1 - 1 (adding + 1 to each distribution and calculating log average of the distributions and subtract 1 from the distribution)
        """
        # copy 'arr_dist'
        l_arr_dist = list( arr.copy( ) for arr in l_arr_dist )
        # add +1
        for arr in l_arr_dist :
            arr += 1   
        # perform log transformation
        l_arr_dist = list( np.log( arr ) for arr in l_arr_dist )
        # calculate average of the combined distributions
        arr_dist_combined = None
        for arr_dist in l_arr_dist :
            arr_dist_combined = _combine_size_distribution( arr_dist_combined, arr_dist )
        arr_dist_combined /= len( l_arr_dist )
        # perform inverse of log transformation, and subtract 1
        arr_dist_combined = np.exp( arr_dist_combined ) - 1
        return arr_dist_combined
    arr_dist_combined = _log_average_of_distributions( l_arr_dist_smoothened ) 

    ''' calculate optimal correction ratios '''
    logger.info(f"calculating the optimal correction ratios")
    def _calculate_ratio( arr_dist_numerator, arr_dist_denominator, value_to_replace_nan = 0 ) :
        len_numerator, len_denominator = len( arr_dist_numerator ), len( arr_dist_denominator )
        # expand the other array whose length is shorter than the other
        if len_numerator > len_denominator :
            arr = np.zeros( len_numerator )
            arr[ : len_denominator ] = arr_dist_denominator
            arr_dist_denominator = arr
        else :
            arr = np.zeros( len_denominator )
            arr[ : len_numerator ] = arr_dist_numerator
            arr_dist_numerator = arr
        res = arr_dist_numerator / arr_dist_denominator
        res[ np.isnan( res ) ] = value_to_replace_nan # replace nan values
        return res

    l_optimal_coefficient = [ ]
    l_arr_ratio_to_ref = [ ]
    slice_distribution_range_of_interest = slice( * t_distribution_range_of_interest ) # define a range of distribution of interest for searching optimal coefficient for calculating normalization ratios
    for arr in l_arr_dist_smoothened :
        def f( ratio ) : # define a function to optimize
            """ # 2023-08-23 14:57:49 
            """
            arr_ratio = _calculate_ratio( arr_dist_combined[ slice_distribution_range_of_interest ], arr[ slice_distribution_range_of_interest ] * ratio ) # use 'slice_distribution_range_of_interest' range for calculating score for searching the optimal cooefficient
            score = np.abs( np.log( arr_ratio[ ( arr_ratio != 0 ) & ( ~ np.isinf( arr_ratio ) ) ] ) ).mean( )
            return score

        grid = (0, 100, 0.1) # define the range for the search
        float_optimal_coefficient = optimize.brute(f, (grid, ))[ 0 ]
        l_optimal_coefficient.append( float_optimal_coefficient )
        l_arr_ratio_to_ref.append( _calculate_ratio( arr_dist_combined, arr * float_optimal_coefficient ) )
        # print( f"'{n}' x {np.round( float_optimal_coefficient, 2 )},\t{np.round( f( 1 ), 2 )} > {np.round( f( float_optimal_coefficient ), 2 )}" ) # print the optimization results

    ''' plot graph (interactive) '''
    # create output folders
    path_folder_graph = f"{path_folder_output}graph/"
    path_folder_graph_noninteractive, path_folder_graph_interactive = f"{path_folder_graph}noninteractive_graph/", f"{path_folder_graph}interactive_graph/"
    for path_folder in [ path_folder_graph, path_folder_graph_noninteractive, path_folder_graph_interactive ] :
        os.makedirs( path_folder, exist_ok = True )
    # display correction ratios
    # compose a dataframe for plotting
    l_l = [ ] # initialize the container
    for arr, name in zip( l_arr_ratio_to_ref, l_name_file_distributions ) :
        for i in range( len( arr ) ) :
            l_l.append( [ name, i, arr[ i ] ] )
    df_ratio = pd.DataFrame( l_l, columns = [ 'name_file_distributions', 'molecule_size_in_base_pairs', 'correction_ratio_to_reference' ] )

    # plot a plotly graph
    fig = px.line( df_ratio, x = 'molecule_size_in_base_pairs', y = 'correction_ratio_to_reference', color = 'name_file_distributions' )
    fig.update_yaxes( range = [ 0, float_max_correction_ratio ] )
    fig.write_html( f'{path_folder_graph_interactive}correction_ratio_to_reference.html' ) # write a html page

    ''' record the size range where correction can be performed confidently  '''
    l_l = [ ]
    for arr, name, path_file in zip( l_arr_ratio_to_ref, l_name_file_distributions, l_path_file_distributions ) :
        arr_pos_invalid = [ 0 ] + list( np.where( arr > float_max_correction_ratio )[ 0 ] ) + [ len( arr ) ] # add start and end of the array as the boundary
        st_range_of_interest, en_range_of_interest = None, None
        for i in range( len( arr_pos_invalid ) - 1 ) :
            st, en = arr_pos_invalid[ i ], arr_pos_invalid[ i + 1 ]
            int_size_overlap = bk.INTERVAL_Overlap( [ st, en ], t_distribution_range_of_interest )
            if int_size_overlap > 0 : # if the overlap exists
                st_range_of_interest, en_range_of_interest = st, en # record the start and end positions of interest.
                break
        l_l.append( [ path_file, name, st_range_of_interest, en_range_of_interest, int_size_overlap, t_distribution_range_of_interest ] )
    df_range_confident = pd.DataFrame( l_l, columns = [ 'path_file', 'name', 'start_range_of_interest', 'end_range_of_interest', 'int_size_overlap', 't_distribution_range_of_interest' ] ) # compose the dataframe
    df_range_confident.to_csv( f'{path_folder_output}df_range_confident.tsv.gz', sep = '\t', index = False ) # save as a file

    ''' plot graph (non-interactive) '''
    for name_col in [ 'start_range_of_interest', 'end_range_of_interest' ] :
        bk.MPL_1D_Sort_Plot( df_range_confident[ name_col ] )
        bk.MPL_SAVE( f"{name_col}", folder = path_folder_graph_noninteractive, l_format=['.pdf', '.png'] )

    ''' export data '''
    logger.info(f"exporting output data")
    # compose a namespace
    dict_output = {
        'setting' : {
            'flag_usage_from_command_line_interface' : flag_usage_from_command_line_interface,
            'l_path_file_distributions' : l_path_file_distributions,
            'l_name_file_distributions' : l_name_file_distributions,
            'float_sigma_gaussian_filter' : float_sigma_gaussian_filter,
            'int_min_total_read_count_for_a_peak' : int_min_total_read_count_for_a_peak,
            'float_min_ratio_read_count_peak_to_baseline' : float_min_ratio_read_count_peak_to_baseline,
            'float_min_ratio_peak_height_to_baseline_height' : float_min_ratio_peak_height_to_baseline_height,
            'int_size_window_surveying_surrounding_values' : int_size_window_surveying_surrounding_values,
            'int_num_iterative_peak_removal' : int_num_iterative_peak_removal,
            'float_max_correction_ratio' : float_max_correction_ratio,
            't_distribution_range_of_interest' : t_distribution_range_of_interest,
            'float_memory_in_GiB' : float_memory_in_GiB,
            'verbose' : verbose,
            'name_type_dist_for_creating_reference' : name_type_dist_for_creating_reference,
            'dict_kw_gaussian_filter' : dict_kw_gaussian_filter,
        },
        'l_arr_dist' : l_arr_dist,
        'l_arr_dist_peak_removed' : l_arr_dist_peak_removed,
        'l_arr_dist_smoothened' : l_arr_dist_smoothened,
        'l_optimal_coefficient' : l_optimal_coefficient,
        'l_arr_ratio_to_ref' : l_arr_ratio_to_ref,
        'arr_dist_combined' : arr_dist_combined,
        'df_range_confident' : df_range_confident,
    }
    bk.PICKLE_Write( f"{path_folder_output}dict_output.pickle", dict_output )
        
    logger.info(f"Completed.")
    return 

def LongExportNormalizedCountMatrix(
    scidx: Union[dict, None] = None,
    flag_usage_from_command_line_interface: bool = False,
    path_folder_ref: Union[str, None] = None,
    path_file_fa_genome: Union[str, None] = None,
    path_file_gtf_genome: Union[str, None] = None,
    path_file_fa_transcriptome: Union[str, None] = None,
    l_path_file_bam_input: Union[list, None] = None,
    l_path_folder_output: [list[str], None] = None,
    path_folder_reference_distribution : Union[ str, None ] = None, # a folder containing the reference distribution, the output of the 'LongCreateReferenceSizeDistribution'
    l_name_distribution : Union[ List[ str ], str, None ] = None, # the name of each sample that was used to build the reference distribution. the distribution of each sample and pre-calculated correction ratios will be retrieved from the data stored in the reference distribution folder using the given names.
    l_str_l_t_distribution_range_of_interest : Union[ List[ str ], str, None ] = None, # define a range of distribution of interest for exporting normalized count matrix
    n_threads: int = 16,
    float_memory_in_GiB: float = 50,
    int_num_sam_records_for_each_chunk: int = 300000,
    str_name_gtf_attr_for_name_transcript: str = "transcript_name",
    int_min_mapq_unique_mapped_for_gex_data: int = 255,
    int_min_mapq_unique_mapped_for_atac_data: int = 60,
    int_n_bases_padding_around_interval: int = 10,
    path_file_tsv_repeatmasker_ucsc: Union[str, None] = None,
    l_repClass_repeatmasker_ucsc: list[str] = [
        "SINE",
        "LINE",
        "LTR",
        "DNA",
        "Retroposon",
    ],
    int_min_length_repeatmasker_ucsc: int = 100,
    flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger: bool = False,
    flag_include_read_aligned_to_opposite_strand: bool = False,
    flag_include_read_aligned_to_intron: bool = False,
    str_name_gtf_attr_for_id_gene: str = "gene_id",
    str_name_gtf_attr_for_name_gene: str = "gene_name",
    str_name_gtf_attr_for_id_transcript: str = "transcript_id",
    path_file_gff_regulatory_element=None,
    str_name_gff_attr_id_regulatory_element="ID",
    flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation: bool = False,
    int_min_length_regulatory_element: int = 50,
    int_bp_padding_regulatory_element_anno: int = 2000,
    float_max_prop_unfiltered_rpmk=1,
    flag_does_not_delete_sequence_and_sequence_qual: bool = False,
    flag_skip_read_analysis_summary_output_bam_file: bool = False,
    flag_skip_read_analysis_summary_output_tsv_file: bool = False,
    flag_turn_off_catching_all_reads_by_binning: bool = False,
    int_bp_for_bins: int = 100,
    flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning: bool = False,
    flag_include_strand_specific_counts: bool = False,
    verbose: bool = True,
    l_str_mode_scarab_count: list[
        Literal[
            "gex5prime-single-end",
            "gex5prime-paired-end",
            "gex3prime-single-end",
            "gex3prime",
            "gex",
            "gex5prime",
            "atac",
            "multiome",
        ]
    ] = ["gex3prime"],
    int_bp_padding_for_defining_promoter_from_transcript_start: int = 2000,
    int_min_mapq_minimap2_tx_assignment : int = 0,
    flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome: bool = False,
    flag_does_not_make_gene_names_unique: bool = False,
    str_name_bam_tag_cb_corrected: str = "CB",
    str_name_bam_tag_cb_uncorrected: str = "CR",
    str_name_bam_tag_umi_corrected: str = "UB",
    str_name_bam_tag_umi_uncorrected: str = "UR",
    flag_skip_exon_and_splice_junc_counting: bool = False,
    path_file_fa_for_cram: Union[str, None] = None,
    int_num_samples_analyzed_concurrently: int = 2,
    flag_does_not_collect_variant_information: bool = False,
    flag_skip_intron_retention_counting: bool = False,
    flag_skip_internal_polyA_filtered_feature_counting: bool = False,
    int_min_length_intron_for_detecting_intron_retention_event: int = 10,
    flag_output_variant_information_with_annotations: bool = False,
    int_min_num_of_reads_for_filtering_genomic_variant: int = 10,
    float_min_prop_of_reads_for_filtering_genomic_variant=0.1,
    path_file_vcf_for_filtering_variant: Union[str, None] = None,
    int_min_count_features_for_filtering_barcodes: int = 50,
    int_length_of_polya_to_append_to_transcript_sequence_during_realignment : int = 50, # during re-alignment analysis for unique transcript assignment, append poly A sequence of given length at the 3' end of transcript sequences, which aids identification of the correct isoform from which the read is likely originated.
    flag_enforce_transcript_end_site_matching_for_long_read_during_realignment : bool = False, # should only be used when (1) all read contains poly A sequences (long-read full-length sequencing results) (2) read is stranded so that its directionality (5'->3') matches that of the original mRNA molecule. For long-read, it is recommanded to turn this setting on. When this mode is active, it also use internal-polyA-tract priming information (the length of the internal poly A tract, recorded as a BAM record tag with the tag name 'str_name_bam_tag_internal_polya_length' for all reads), and does not perform TES matching if the read appear to be primed by internal-polyA-tract. To enable this behavior, 'int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment' should be set to non-negative values, and alignments to transcripts should be filtered based on the softclipping status of the alignment.
    str_name_bam_tag_internal_polya_length : str = 'IA', # the name of the BAM record tag that contains the length of internal poly A tract. The tag should be available for all reads if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is set to True, and TES matching mode is active.
    int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment : int = 10, # Rather than aligning an entire sequence of the read, exclude soft clipped regions and align the portion of read that was aligned to the genome. Since this portion of read should be perfectly match the transcript without softclipping if the read was indeed originated from the transcript, during realignment, alignments with extensive softclipping longer than the given threshold will be filtered out. Additionally, alignment to transcript with insertion and deletion longer than this length will be filtered out, too. To disable this behavior, set this value to negative values (e.g., -1).
    int_max_distance_from_transcript_start_and_end_for_tss_and_tes_matching_during_realignment : int = 100, # the maximum distance (in base pairs, bp) from the transcript start or end coordinates for a read to be classified as a specific transcript. Currently, only TES matching is supported, due to 5' frequent degradation events. This argument will be only effective if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is True.
    str_mappy_aligner_preset_for_realignment : str = 'map-ont', # minimap2 presets for re-alignment analysis: 'sr' for single-end short read data; 'map-pb' for PacBio long read data; 'map-ont' for Oxford Nanopore long read data. Please avoid using the 'splice' preset, since re-alignment to transcripts should not contain 'splicing', or large deletions.
    dict_num_manager_processes_for_each_data_object: dict = {
        'dict_it_promoter' : 0,
        'dict_t_splice_junc_to_info_genome' : 0,
        'dict_it_exon' : 0,
        'dict_it_exon_transcriptome' : 3,
        'dict_it_splice_junc_transcriptome' : 3,
        'dict_it_splice_donor_and_acceptor_genome' : 3,
        'dict_it_rpmk' : 5,
        'dict_it_reg' : 3,
        'dict_fa_transcriptome' : 2,
    },
    l_seqname_to_skip: list = ["MT"],
    flag_no_strand_specificity : bool = False,
) -> dict:
    """# 2023-11-03 17:30:23 
    perform secondary analysis of cell-ranger output (barcoded BAM)

    l_str_mode_scarab_count : list[ Literal[ "gex5prime-single-end", 'gex5prime-paired-end', "gex3prime-single-end", 'gex3prime', 'gex', 'gex5prime', 'atac', 'multiome' ] ] = [ 'gex3prime' ], # list of scarab_count operation mode

    scidx : Union[ dict, None ] = None, # a loaded scarab index object. if given, the object will be used instead loading the index from the disk.
    path_file_fa_for_cram : Union[ str, None ] = None, # path to the fasta file used for CRAM. If the fasta file has not been indexed, it will be automatically indexed
    flag_does_not_collect_variant_information : bool = False, # does not collect the variant information at all. it will improve performance at the cost of the reduced output information.
    flag_output_variant_information_with_annotations : bool = False, # If True, record variants for each individual feature (gene, isoform, genome bin, etc.). If False, variant information will be recorded for only the 'variant' feature type (require 'path_file_vcf_for_filtering_variant' argument to be active).
    int_min_num_of_reads_for_filtering_genomic_variant : int = 10, # for a variant to be exported as a feature, at least this number of reads should be present for the variant
    float_min_prop_of_reads_for_filtering_genomic_variant = 0.1, # for a variant to be exported as a feature, at least this proportion of reads should be present for the variant
    path_file_vcf_for_filtering_variant : Union[ str, None ] = None, # A path to the vcf file for filtering variants. When a valid VCF file is given, variant filtering criteria, 'float_min_prop_of_reads_for_filtering_genomic_variant' and 'int_min_num_of_reads_for_filtering_genomic_variant' will be ignored. Also, a new feature type 'variant' will be added in the count matrix containing coverate of each variant and its reference allele at single-cell level. (warning) due to the internal algorithm for distributing workloads across the workers, count records for 'variant' features can be duplicated (matrix contains more than one records describing counts of a unique pair of cell and feature).
    int_min_count_features_for_filtering_barcodes : int = 50, # the minimum number of features for filtering barcodes.
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    l_seqname_to_skip : list = [ 'MT' ], # the list of names of the chromosomes of the reference genome to skip the analysis. By default, reads aligned to the mitochondrial genomes will be skipped. Because gene boundaries of the mitochondrial genome-encoded genes are often overlapping, an entire mitochondrial genome often assigned as a single chunk, creating a huge bottleneck in the analysis pipeline.
    flag_skip_internal_polyA_filtered_feature_counting: bool = False, # skip exporting counts excluding internal polyA primed reads (if it is set to False, in addition to the default counting mode, internal polyA primed-read filtered reads will be counted for most of the features).
    flag_skip_intron_retention_counting: bool = False, # skip exporting counts of intron retention events
    int_min_length_intron_for_detecting_intron_retention_event: int = 10, # the minimum length of intron to be present in a read in order to detect an intron retention event in the read.
    dict_num_manager_processes_for_each_data_object : dict = {
        'dict_it_promoter' : 0,
        'dict_t_splice_junc_to_info_genome' : 0,
        'dict_it_exon' : 0,
        'dict_it_exon_transcriptome' : 3,
        'dict_it_splice_junc_transcriptome' : 3,
        'dict_it_splice_donor_and_acceptor_genome' : 3,
        'dict_it_rpmk' : 5,
        'dict_it_reg' : 3,
        'dict_fa_transcriptome' : 2,
    }
    flag_no_strand_specificity : bool = False, # flag indicating whether to ignore strand information of the reads in the input BAM files.
    path_folder_reference_distribution : Union[ str, None ] = None, # a folder containing the reference distribution, the output of the 'LongCreateReferenceSizeDistribution'
    l_name_distribution : Union[ List[ str ], str, None ] = None, # the name of each sample that was used to build the reference distribution. the distribution of each sample and pre-calculated correction ratios will be retrieved from the data stored in the reference distribution folder using the given names.
    l_str_l_t_distribution_range_of_interest : Union[ List[ str ], str, None ] = None, # define a range of distribution of interest for exporting normalized count matrix. a list of string for setting the size distrubution ranges of interest for exporting normalized count matrix. if 'raw' is given, no size-based normalization will be performed, and raw counts of all molecules will be exported. example arguments are the followings: 'raw,50-5000,1000-3500' for exporting raw count and size-normalized count matrices for molecules of 50-5000bp and 1000-3500bp (total three output matrices). if only one argument is given, the argument will be applied to all samples.

    ------ for re-alignment -------
    int_min_mapq_minimap2_tx_assignment : int = 0, # (default = 0, meaning no filtering). a value between 0~60. Minimum mapping quality required for assigning reads to a unique isoform using minimap2 realignment.
    int_length_of_polya_to_append_to_transcript_sequence_during_realignment : int = 50, # during re-alignment analysis for unique transcript assignment, append poly A sequence of given length at the 3' end of transcript sequences, which aids identification of the correct isoform from which the read is likely originated.
    flag_enforce_transcript_end_site_matching_for_long_read_during_realignment : bool = False, # hould only be used when (1) all read contains poly A sequences (long-read full-length sequencing results) (2) read is stranded so that its directionality (5'->3') matches that of the original mRNA molecule. For long-read, it is recommanded to turn this setting on. When this mode is active, it also use internal-polyA-tract priming information (the length of the internal poly A tract, recorded as a BAM record tag with the tag name 'str_name_bam_tag_internal_polya_length' for all reads), and does not perform TES matching if the read appear to be primed by internal-polyA-tract. To enable this behavior, 'int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment' should be set to non-negative values, and alignments to transcripts should be filtered based on the softclipping status of the alignment.
    int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment : int = 10, # Rather than aligning an entire sequence of the read, exclude soft clipped regions and align the portion of read that was aligned to the genome. Since this portion of read should be perfectly match the transcript without softclipping if the read was indeed originated from the transcript, during realignment, alignments with extensive softclipping longer than the given threshold will be filtered out. Additionally, alignment to transcript with insertion and deletion longer than this length will be filtered out, too. To disable this behavior, set this value to negative values (e.g., -1).
    str_name_bam_tag_internal_polya_length : str = 'IA', # the name of the BAM record tag that contains the length of internal poly A tract. The tag should be available for all reads if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is set to True, and TES matching mode is active.
    int_max_distance_from_transcript_start_and_end_for_tss_and_tes_matching_during_realignment : int = 100, # the maximum distance (in base pairs, bp) from the transcript start or end coordinates for a read to be classified as a specific transcript. Currently, only TES matching is supported, due to 5' frequent degradation events. This argument will be only effective if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is True.
    str_mappy_aligner_preset_for_realignment : str = 'sr', # minimap2 presets for re-alignment analysis: 'sr' for single-end short read data; 'map-pb' for PacBio long read data; 'map-ont' for Oxford Nanopore long read data. Please avoid using the 'splice' preset, since re-alignment to transcripts should not contain 'splicing', or large deletions.
    
    # the number of manager processes to use for each data object that will be shared across the forked processes. If 0 is given, no manager process will be used. Instead, the object will be directly accessed in the forked process, incurring memory bloating.
    # generally, it is better to use more number of manager processes for data object that are more frequently accessed. If increasing the number of manager processes does not improve performance, considering not using the manager process and accessing the object directly.
    # the expected size of bloated memory per process for each data object is given below.
    #
    #   'object name'                                       'the size of bloated memory per process'
    #   dict_it_exon_transcriptome                          1.617437 GB per process
    #   dict_it_rpmk                                        1.452151 GB per process
    #   dict_it_splice_junc_transcriptome                   1.381314 GB per process
    #   dict_it_splice_donor_and_acceptor_genome            ???????? GB per process (not measured)
    #   dict_fa_transcriptome                               0.460438 GB per process
    #   dict_it_exon                                        0.271540 GB per process
    #   dict_it_reg                                         0.271540 GB per process
    #   dict_t_splice_junc_to_info_genome                   0.188898 GB per process
    #   dict_it_promoter                                    0.141673 GB per process
    #   dict_fa_genome                                      0.082643 GB per process
    #   dict_id_tx_to_id_gene                               0.070837 GB per process
    #   dict_id_tx_to_name_tx                               0.070837 GB per process
    #   dict_it_gene                                        0.059031 GB per process
    #   dict_id_gene_to_l_id_tx                             0.059031 GB per process
    #   dict_index_df_gtf_gene                              0.047224 GB per process
    #   arr_data_df_gtf_gene                                0.047224 GB per process
    #   dict_seqname_to_mask_gtf_reg                        0.035418 GB per process
    #   dict_seqname_to_mask_gtf_intron_near_splice_site    0.035418 GB per process
    #   dict_seqname_to_mask_gtf_rpmk_unfiltered            0.035418 GB per process
    #   dict_seqname_to_mask_gtf_rpmk_filtered              0.035418 GB per process
    #   dict_seqname_to_mask_gtf_exon                       0.035418 GB per process
    #
    # if pre-loaded 'scidx' is given, this argument will be ignored.

    returns
    a loaded scarab-count index ('scidx')
    """
    """
    Parse arguments
    """
    if flag_usage_from_command_line_interface:  # parse arguments
        """parse arguments when the function was called from the command-line interface"""
        # { 'K', 'k' } # unused arguments
        # command line arguments
        parser = argparse.ArgumentParser(
            description=str_description,
            usage="scarab count",
            formatter_class=argparse.RawTextHelpFormatter,
        )
        parser.add_argument("count")
        arg_grp_general = parser.add_argument_group("General")
        arg_grp_general.add_argument(
            "-b",
            "--l_path_file_bam_input",
            help="A barcoded BAM file (or a list of such files, separated by spaces), sorted by read alignment position, from cellranger (GEX, ATAC, or ARC for multiome data) or similar pipelines. To process more than one samples, paths of multiple BAM files can be given. For a multiome sample, two barcoded BAM file paths (the BAM file path for the GEX data comes first, followed by the BAM file path for the ATAC data) should be given for each sample (with an appropriate argument for the '--l_str_mode_scarab_count' argument)",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-o",
            "--l_path_folder_output",
            help="(default: scarab_output/ subdirectory of the folder where an input BAM file for the sample resides. For a multiome sample, scarab_output/ subdirectory of the folder where an input GEX BAM file resides) a directory of an output folder or a list of output folders, separated by spaces. The number of output folders should be same as the number of samples",
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-i",
            "--l_str_mode_scarab_count",
            help="one of [ 'gex', 'atac', 'multiome' ] (default: 'gex') an operating mode of scarab short-read. 'gex' for gene expression profiling, 'atac' for chromatin accessibility profiling, 'multiome' for gene expression profiling and chromatin accessibility profiling for the same single cells. If more than a single scarab count mode is used for different samples (e.g. 'gex' for the first sample and 'atac' for the second sample), multiple arguments separated by a space can be given. If multiple operating modes are given, the number of modes should be same as the number of samples",
            default=["gex"],
            nargs="*",
        )
        arg_grp_general.add_argument(
            "-a",
            "--path_folder_ref",
            help="path_folder_ref (Default: ref/ folder inside the 'path_folder_output' folder) if given, instead of processing and building annotation data from the input files, the reference data saved in the given folder will be used to load annotations required for analyzing reads.",
        )
        arg_grp_general.add_argument(
            "-t",
            "--n_threads",
            help="n_threads (default : 16) scarab-count pipeline is scalable upto ~50 threads (processes). The use of Python Multiprocessing Manager process to share some of data across workers limits the number of processes.",
            default=16,
            type=int,
        )
        arg_grp_general.add_argument(
            "-j",
            "--float_memory_in_GiB",
            help="float_memory_in_GiB",
            default=70,
            type=float,
        )
        arg_grp_general.add_argument(
            "--int_min_count_features_for_filtering_barcodes",
            help="int_min_count_features_for_filtering_barcodes",
            default=50,
            type=int,
        )
        arg_grp_general.add_argument(
            "--int_num_samples_analyzed_concurrently",
            help="the number of samples that can be analyzed concurrently.",
            default=2,
            type=int,
        )
        arg_grp_general.add_argument(
            "--dict_num_manager_processes_for_each_data_object",
            help="(default: { 'dict_it_promoter' : 0, 'dict_t_splice_junc_to_info_genome' : 0, 'dict_it_exon' : 0, 'dict_it_exon_transcriptome' : 2, 'dict_it_splice_junc_transcriptome' : 2, 'dict_it_rpmk' : 3, 'dict_it_reg' : 2, 'dict_fa_transcriptome' : 2, }) the number of manager processes to use for each data object that will be shared across the forked processes. If 0 is given, no manager process will be used. Instead, the object will be directly accessed in the forked process, incurring memory bloating.",
            default="{ 'dict_it_promoter' : 0, 'dict_t_splice_junc_to_info_genome' : 0, 'dict_it_exon' : 0, 'dict_it_exon_transcriptome' : 2, 'dict_it_splice_junc_transcriptome' : 2, 'dict_it_rpmk' : 3, 'dict_it_reg' : 2, 'dict_fa_transcriptome' : 2, }",
            type=str,
        )
        arg_grp_general.add_argument(
            "-B", "--verbose", help="turn on verbose mode", action="store_true"
        )

        arg_grp_annotation_gene = parser.add_argument_group("Annotation - Genes")
        arg_grp_annotation_gene.add_argument(
            "-g",
            "--path_file_fa_genome",
            help="path_file_fa_genome. Using Ensembl version is highly recommended (chromosome name should not contain 'chr')",
        )
        arg_grp_annotation_gene.add_argument(
            "--l_seqname_to_skip",
            help="(default: [ 'MT' ]) the list of names of the chromosomes of the reference genome to skip the analysis. By default, reads aligned to the mitochondrial genomes will be skipped. Because gene boundaries of the mitochondrial genome-encoded genes are often overlapping, an entire mitochondrial genome often assigned as a single chunk, creating a huge bottleneck in the analysis pipeline.",
            default=["MT"],
            nargs="*",
        )
        arg_grp_annotation_gene.add_argument(
            "-G",
            "--path_file_gtf_genome",
            help="path_file_gtf_genome Using Ensembl annotation is highly recommended (chromosome name should not contain 'chr')",
        )
        arg_grp_annotation_gene.add_argument(
            "-T", "--path_file_fa_transcriptome", help="path_file_fa_transcriptome"
        )
        arg_grp_annotation_gene.add_argument(
            "-J",
            "--flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome",
            help="flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome (Default: False) does not drop the integer 'version_info' after the '.' at the end of the 'id_transcript' from the given fasta file containing transcriptome sequences",
            action="store_true",
        )

        # process gene annotations
        arg_grp_annotation_gene.add_argument(
            "-D",
            "--str_name_gtf_attr_for_id_gene",
            help="str_name_gtf_attr_for_id_gene",
            default="gene_id",
        )
        arg_grp_annotation_gene.add_argument(
            "-N",
            "--str_name_gtf_attr_for_name_gene",
            help="str_name_gtf_attr_for_name_gene",
            default="gene_name",
        )
        arg_grp_annotation_gene.add_argument(
            "-Y",
            "--str_name_gtf_attr_for_id_transcript",
            help="str_name_gtf_attr_for_id_transcript",
            default="transcript_id",
        )
        arg_grp_annotation_gene.add_argument(
            "-Z",
            "--str_name_gtf_attr_for_name_transcript",
            help="str_name_gtf_attr_for_name_transcript",
            default="transcript_name",
        )

        # modify gene annotations for more accurate analysis in downstream applications (e.g. scanpy)
        arg_grp_annotation_gene.add_argument(
            "-u",
            "--flag_does_not_make_gene_names_unique",
            help="(Default: False) flag_does_not_make_gene_names_unique",
            action="store_true",
        )
        arg_grp_annotation_gene.add_argument(
            "-c",
            "--flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation",
            help="(Default: False) flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation",
            action="store_true",
        )
        arg_grp_annotation_gene.add_argument(
            "-p",
            "--int_bp_padding_for_defining_promoter_from_transcript_start",
            help="int_bp_padding_for_defining_promoter_from_transcript_start (default : 2000) the number of base pairs from the transcript start site for defining promoter regions of a gene",
            default=2000,
            type=int,
        )
        
        arg_grp_isoform_realignment = parser.add_argument_group("Isoform assignment - Realignment")
        arg_grp_isoform_realignment.add_argument(
            "--int_max_distance_from_transcript_start_and_end_for_tss_and_tes_matching_during_realignment",
            help="(Default: 100) The maximum distance (in base pairs, bp) from the transcript start or end coordinates for a read to be classified as a specific transcript. Currently, only TES matching is supported, due to 5' frequent degradation events. This argument will be only effective if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is True.",
            default=100,
            type=int,
        )    
        arg_grp_isoform_realignment.add_argument(
            "--str_mappy_aligner_preset_for_realignment",
            help="(Default: 'sr') minimap2 presets for re-alignment analysis: 'sr' for single-end short read data; 'map-pb' for PacBio long read data; 'map-ont' for Oxford Nanopore long read data. Please avoid using the 'splice' preset, since re-alignment to transcripts should not contain 'splicing', or large deletions.",
            default = 'sr',
        )    
        arg_grp_isoform_realignment.add_argument(
            "--int_length_of_polya_to_append_to_transcript_sequence_during_realignment",
            help="(Default: 50) During re-alignment analysis for unique transcript assignment, append poly A sequence of given length at the 3' end of transcript sequences, which aids identification of the correct isoform from which the read is likely originated.",
            default=50,
            type=int,
        )
        arg_grp_isoform_realignment.add_argument(
            "--flag_enforce_transcript_end_site_matching_for_long_read_during_realignment",
            help="(Default: False) Should only be used when (1) all read contains poly A sequences (long-read full-length sequencing results) (2) read is stranded so that its directionality (5'->3') matches that of the original mRNA molecule. For long-read, it is recommanded to turn this setting on. When this mode is active, it also use internal-polyA-tract priming information (the length of the internal poly A tract, recorded as a BAM record tag with the tag name 'str_name_bam_tag_internal_polya_length' for all reads), and does not perform TES matching if the read appear to be primed by internal-polyA-tract. To enable this behavior, 'int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment' should be set to non-negative values, and alignments to transcripts should be filtered based on the softclipping status of the alignment.",
            action="store_true",
        )
        arg_grp_isoform_realignment.add_argument(
            "--int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment",
            help="(Default: 10) Rather than aligning an entire sequence of the read, exclude soft clipped regions and align the portion of read that was aligned to the genome. Since this portion of read should be perfectly match the transcript without softclipping if the read was indeed originated from the transcript, during realignment, alignments with extensive softclipping longer than the given threshold will be filtered out. Additionally, alignment to transcript with insertion and deletion longer than this length will be filtered out, too. To disable this behavior, set this value to negative values (e.g., -1).",
            default=10,
            type=int,
        )    
        arg_grp_isoform_realignment.add_argument(
            "--str_name_bam_tag_internal_polya_length",
            help="(Default: 'IA') The name of the BAM record tag that contains the length of internal poly A tract. The tag should be available for all reads if 'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' is set to True, and TES matching mode is active.",
            default = "IA",
        )    

        # regulatory elements
        arg_grp_annotation_reg = parser.add_argument_group(
            "Annotation - Regulatory Elements"
        )
        arg_grp_annotation_reg.add_argument(
            "-F",
            "--path_file_gff_regulatory_element",
            help="path_file_gff_regulatory_element",
        )
        arg_grp_annotation_reg.add_argument(
            "-R",
            "--str_name_gff_attr_id_regulatory_element",
            help="str_name_gff_attr_id_regulatory_element",
            default="ID",
        )
        arg_grp_annotation_reg.add_argument(
            "-M",
            "--int_min_length_regulatory_element",
            help="int_min_length_regulatory_element",
            default=50,
            type=int,
        )
        arg_grp_annotation_reg.add_argument(
            "-E",
            "--int_bp_padding_regulatory_element_anno",
            help="int_bp_padding_regulatory_element_anno",
            default=2000,
            type=int,
        )
        arg_grp_annotation_reg.add_argument(
            "-v",
            "--float_max_prop_unfiltered_rpmk",
            help="float_max_prop_unfiltered_rpmk (default: 1) (i.e. no filtering based on the proportion of reads overlapped with unfiltered repeatmasker elements)",
            default=1,
            type=float,
        )

        # repeatmasker
        arg_grp_annotation_rpmk = parser.add_argument_group(
            "Annotation - Repeats (UCSC)"
        )
        arg_grp_annotation_rpmk.add_argument(
            "-U",
            "--path_file_tsv_repeatmasker_ucsc",
            help="path_file_tsv_repeatmasker_ucsc: a TSV file downloaded from UCSC Table Browser (all fields)",
        )
        arg_grp_annotation_rpmk.add_argument(
            "-C",
            "--l_repClass_repeatmasker_ucsc",
            help="(default: [ 'SINE', 'LINE', 'LTR', 'DNA', 'Retroposon' ]) l_repClass_repeatmasker_ucsc: list of repClass in the given UCSC repeatmasker annotations to be analyzed",
            default=["SINE", "LINE", "LTR", "DNA", "Retroposon"],
            nargs="*",
        )
        arg_grp_annotation_rpmk.add_argument(
            "-l",
            "--int_min_length_repeatmasker_ucsc",
            help="int_min_length_repeatmasker_ucsc: ignore repeatmasker annotations with length smaller than the given length",
            default=100,
            type=int,
        )

        # repeatmasker
        arg_grp_annotation_rpmk = parser.add_argument_group(
            "Annotation - Catch-All (binning)"
        )
        arg_grp_annotation_rpmk.add_argument(
            "-W",
            "--flag_turn_off_catching_all_reads_by_binning",
            help="flag_turn_off_catching_all_reads_by_binning: Does not collect read counts of the reads confidently aligned to the reference genome and count the reads based on the orientation and position of the alignment for each genomic bins",
            action="store_true",
        )
        arg_grp_annotation_rpmk.add_argument(
            "-L",
            "--int_bp_for_bins",
            help="(default: 100) int_bp_for_bins: number of base pairs for each genomic bin for binning. For example, when 'int_bp_for_bins' is 500, genomic bins will be chr1:1-500, chr1:501-1000, ... ... ",
            default=100,
            type=int,
        )
        arg_grp_annotation_rpmk.add_argument(
            "-I",
            "--flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning",
            help="flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning: if True, reads assigned to features (both gene/miscellaneous) will be excluded when counting reads for each genomic bin. When this flag is not set, EVERY reads confidently aligned to the reference genome will be counted for each corresponding genomic bin",
            action="store_true",
        )

        # for exporting count matrix containing genomic variation information
        arg_grp_variant = parser.add_argument_group("Variant Calling / Filtering")
        arg_grp_variant.add_argument(
            "-V",
            "--flag_does_not_collect_variant_information",
            help="(Default: False) set this flag to True to reduce the output file size and speed up the run by not analyzing variant information present in the sequencing data. The variant information collected and stored in BAM file can be retrieved easily later to look for possible SNPs and allele specific expression",
            action="store_true",
        )
        arg_grp_variant.add_argument(
            "--flag_output_variant_information_with_annotations",
            help="(Default: False) If True, record variants for each individual feature (gene, isoform, genome bin, etc.). If False, only the minimal number of features for recording variants (the number of reads covering the variant for each cell) will be included in the output. Setting this flag to True will slightly decrease the performance",
            action="store_true",
        )
        arg_grp_variant.add_argument(
            "-f",
            "--int_min_num_of_reads_for_filtering_genomic_variant",
            help="(int_min_num_of_reads_for_filtering_genomic_variant) (default: 10) minimum number of UMI counts required to include the information about the detected genomic variation in the output count matrix. Setting this value to a negative number (e.g. -1) will disable the behavior of 'scarab' exporting genomic variation information into the count matrix, which can significantly increase the run time and the size of the output count matrix for a typical 'int_min_num_of_reads_for_filtering_genomic_variant' value between 1~10. ",
            default=10,
            type=int,
        )
        arg_grp_variant.add_argument(
            "-d",
            "--float_min_prop_of_reads_for_filtering_genomic_variant",
            help="(default: 0.1) minimum proportion of UMI counts required to include the information about the detected genomic variation in the output count matrix. The variant should be above both thresholds, 'float_min_prop_of_reads_for_filtering_genomic_variant' and 'int_min_num_of_reads_for_filtering_genomic_variant'",
            default=0.1,
            type=float,
        )
        arg_grp_variant.add_argument(
            "-K",
            "--path_file_vcf_for_filtering_variant",
            help="A path to the vcf file for filtering variants. When a valid VCF file is given, variant filtering criteria, 'float_min_prop_of_reads_for_filtering_genomic_variant' and 'int_min_num_of_reads_for_filtering_genomic_variant' will be ignored. Also, a new feature type 'variant' will be added in the count matrix containing coverate of each variant and its reference allele at single-cell level. (warning) due to the internal algorithm for distributing workloads across the workers, count records for 'variant' features can be duplicated (matrix contains more than one records describing counts of a unique pair of cell and feature).",
        )

        # BAM file processing
        arg_grp_bam_processing = parser.add_argument_group(
            "Barcoded BAM File Processing"
        )
        arg_grp_bam_processing.add_argument(
            "-H",
            "--path_file_fa_for_cram",
            help="path to the fasta file used for CRAM. If the fasta file has not been indexed, it will be automatically indexed.",
        )
        arg_grp_bam_processing.add_argument(
            "-n",
            "--int_num_sam_records_for_each_chunk",
            help="int_num_sam_records_for_each_chunk: default = 300000. the minimum number of SAM records for each chunk that will be processed by a single process during multiprocessing. For a smaller BAM file (1~3 GB), decreasing this number below 100000 reads are recommended to make use of multiple threads of the machine.",
            default=300000,
            type=int,
        )
        arg_grp_bam_processing.add_argument(
            "-m",
            "--int_min_mapq_unique_mapped_for_gex_data",
            help="int_min_mapq_unique_mapped_for_gex_data: default = 255. Minimum mapping quality of aligned reads in the input BAM file (Gene Expression data) to be included in the analysis. For gene expression data (BAM files aligned with STAR Aligner), recommended value is 255.",
            default=255,
            type=int,
        )
        arg_grp_bam_processing.add_argument(
            "-r",
            "--int_min_mapq_unique_mapped_for_atac_data",
            help="int_min_mapq_unique_mapped_for_atac_data: default = 60. Minimum mapping quality of aligned reads in the input BAM file (ATAC-seq data) to be included in the analysis. For ATAC data (BAM file aligned with BWA Aligner), recommended value is 60.",
            default=60,
            type=int,
        )
        arg_grp_bam_processing.add_argument(
            "-e",
            "--int_min_mapq_minimap2_tx_assignment",
            help="(default = 0, meaning no filtering). a value between 0~60. Minimum mapping quality required for assigning reads to a unique isoform using minimap2 realignment.",
            default=0,
            type=int,
        )
        arg_grp_bam_processing.add_argument(
            "-k",
            "--flag_skip_exon_and_splice_junc_counting",
            help="flag_skip_exon_and_splice_junc_counting: (Default = False) Skip exon and splice_junciton detection and counting. This will slighly increase the processing time and the smaller number of features in the output count matrix. when 'flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger' is True, the exon and splice-junction counting behaviors will be disabled",
            action="store_true",
        )

        # for instructing barcoded BAM tag names
        arg_grp_bam_processing.add_argument(
            "-x",
            "--str_name_bam_tag_cb_corrected",
            help="str_name_bam_tag_cb_corrected",
            default="CB",
        )
        arg_grp_bam_processing.add_argument(
            "-y",
            "--str_name_bam_tag_cb_uncorrected",
            help="str_name_bam_tag_cb_uncorrected",
            default="CR",
        )
        arg_grp_bam_processing.add_argument(
            "-z",
            "--str_name_bam_tag_umi_corrected",
            help="str_name_bam_tag_umi_corrected",
            default="UB",
        )
        arg_grp_bam_processing.add_argument(
            "-w",
            "--str_name_bam_tag_umi_uncorrected",
            help="str_name_bam_tag_umi_uncorrected",
            default="UR",
        )

        # count matrix output
        arg_grp_bam_processing.add_argument(
            "-s",
            "--flag_include_strand_specific_counts",
            help="flag_include_strand_specific_counts: in addition to the typical gene and isoform counts, include count of reads aligned in antisense direction as a separate features. For other features, including repeats, regulatory elements, and genomic regions, include separate features for sense and antisense reads",
            action="store_true",
        )
        
        arg_grp_bam_processing.add_argument(
            "--flag_skip_intron_retention_counting",
            help="(Default: False) set this flag to True to skip exporting counts of intron retention events",
            action="store_true",
        )
        
        arg_grp_bam_processing.add_argument(
            "--flag_skip_internal_polyA_filtered_feature_counting",
            help="(Default: False) skip exporting counts excluding internal polyA primed reads (if it is set to False, in addition to the default counting mode, internal polyA primed-read filtered reads will be counted for most of the features).",
            action="store_true",
        )
        
        arg_grp_bam_processing.add_argument(
            "--int_min_length_intron_for_detecting_intron_retention_event",
            help="(Default: 10). The minimum length of intron to be present in a read in order to detect an intron retention event in the read.",
            default=10,
            type=int,
        )
        
        arg_grp_bam_processing.add_argument(
            "--flag_no_strand_specificity",
            help="(Default: False) set this flag to True to consider reads in the input BAM file as reads without strand-specificity information",
            action="store_true",
        )

        # for multiprocessing (bookmarks)
        arg_grp_bam_processing.add_argument(
            "-P",
            "--int_n_bases_padding_around_interval",
            help="int_n_bases_padding_around_interval",
            default=10,
            type=int,
        )

        # 10X cellranger count behavior
        arg_grp_bam_processing.add_argument(
            "-X",
            "--flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger",
            help="(Default: False) flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger. Setting this flag to True will disable exon and splice-junction counting behavior, due to a possible inconsistency between reference annotations used in cellranger and scarab",
            action="store_true",
        )
        arg_grp_bam_processing.add_argument(
            "-O",
            "--flag_include_read_aligned_to_opposite_strand",
            help="(Default: False) flag_include_read_aligned_to_opposite_strand",
            action="store_true",
        )
        arg_grp_bam_processing.add_argument(
            "-A",
            "--flag_include_read_aligned_to_intron",
            help="(Default: False) flag_include_read_aligned_to_intron",
            action="store_true",
        )

        # setting for output annotated BAM file
        arg_grp_bam_output = parser.add_argument_group("Annotated BAM Output")
        arg_grp_bam_output.add_argument(
            "-Q",
            "--flag_skip_read_analysis_summary_output_bam_file",
            help="(Default: False) set this flag to True in order to skip writing the analysis results of individual reads in the BAM file format",
            action="store_true",
        )
        arg_grp_bam_output.add_argument(
            "-q",
            "--flag_skip_read_analysis_summary_output_tsv_file",
            help="(Default: False) set this flag to True in order to skip writing the analysis results of individual reads in the TSV file format",
            action="store_true",
        )
        arg_grp_bam_output.add_argument(
            "-S",
            "--flag_does_not_delete_sequence_and_sequence_qual",
            help="(Default: False) set this flag to True in order to disable a behavior that removes seq and qual records from SAM records to reduce the output file size",
            action="store_true",
        )
        
        # setting for exporting size-based normalized count matrix
        arg_grp_size_norm = parser.add_argument_group("Count Normalization Using the Reference Distribution")
        arg_grp_size_norm.add_argument(
            "--path_folder_reference_distribution",
            help="path_folder_reference_distribution. a folder containing the reference distribution, the output of the 'LongCreateReferenceSizeDistribution'",
        )
        arg_grp_size_norm.add_argument(
            "--l_name_distribution",
            help="the name of each sample that was used to build the reference distribution. the distribution of each sample and pre-calculated correction ratios will be retrieved from the data stored in the reference distribution folder using the given names.",
            default=None,
            nargs="*",
        )
        arg_grp_size_norm.add_argument(
            "--l_str_l_t_distribution_range_of_interest",
            help="a list of string for setting the size distrubution ranges of interest for exporting normalized count matrix. if 'raw' is given, no size-based normalization will be performed, and raw counts of all molecules will be exported. example arguments are the followings: 'raw,50-5000,1000-3500' for exporting raw count and size-normalized count matrices for molecules of 50-5000bp and 1000-3500bp (total three output matrices). if only one argument is given, the argument will be applied to all samples.",
            default=["raw"],
            nargs="*",
        )
        args = parser.parse_args()

        flag_skip_exon_and_splice_junc_counting = (
            args.flag_skip_exon_and_splice_junc_counting
        )
        flag_include_strand_specific_counts = args.flag_include_strand_specific_counts
        flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome = (
            args.flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome
        )
        int_min_mapq_minimap2_tx_assignment = args.int_min_mapq_minimap2_tx_assignment
        int_bp_padding_for_defining_promoter_from_transcript_start = (
            args.int_bp_padding_for_defining_promoter_from_transcript_start
        )
        int_num_sam_records_for_each_chunk = args.int_num_sam_records_for_each_chunk
        int_min_mapq_unique_mapped_for_gex_data = (
            args.int_min_mapq_unique_mapped_for_gex_data
        )
        int_min_mapq_unique_mapped_for_atac_data = (
            args.int_min_mapq_unique_mapped_for_atac_data
        )
        l_str_mode_scarab_count = args.l_str_mode_scarab_count
        verbose = args.verbose
        l_path_file_bam_input = args.l_path_file_bam_input
        l_path_folder_output = args.l_path_folder_output
        path_file_fa_genome = args.path_file_fa_genome
        path_file_gtf_genome = args.path_file_gtf_genome
        path_file_fa_transcriptome = args.path_file_fa_transcriptome
        path_folder_ref = args.path_folder_ref
        n_threads = args.n_threads
        float_memory_in_GiB = args.float_memory_in_GiB
        int_n_bases_padding_around_interval = args.int_n_bases_padding_around_interval
        int_min_num_of_reads_for_filtering_genomic_variant = (
            args.int_min_num_of_reads_for_filtering_genomic_variant
        )
        path_file_tsv_repeatmasker_ucsc = args.path_file_tsv_repeatmasker_ucsc
        l_repClass_repeatmasker_ucsc = args.l_repClass_repeatmasker_ucsc
        int_min_length_repeatmasker_ucsc = args.int_min_length_repeatmasker_ucsc
        flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger = (
            args.flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger
        )
        flag_include_read_aligned_to_opposite_strand = (
            args.flag_include_read_aligned_to_opposite_strand
        )
        flag_include_read_aligned_to_intron = args.flag_include_read_aligned_to_intron
        str_name_gtf_attr_for_id_gene = args.str_name_gtf_attr_for_id_gene
        str_name_gtf_attr_for_name_gene = args.str_name_gtf_attr_for_name_gene
        str_name_gtf_attr_for_id_transcript = args.str_name_gtf_attr_for_id_transcript
        str_name_gtf_attr_for_name_transcript = (
            args.str_name_gtf_attr_for_name_transcript
        )
        path_file_gff_regulatory_element = args.path_file_gff_regulatory_element
        str_name_gff_attr_id_regulatory_element = (
            args.str_name_gff_attr_id_regulatory_element
        )
        int_min_length_regulatory_element = args.int_min_length_regulatory_element
        int_bp_padding_regulatory_element_anno = (
            args.int_bp_padding_regulatory_element_anno
        )
        float_max_prop_unfiltered_rpmk = args.float_max_prop_unfiltered_rpmk
        flag_does_not_make_gene_names_unique = args.flag_does_not_make_gene_names_unique
        flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation = (
            args.flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation
        )
        str_name_bam_tag_cb_corrected = args.str_name_bam_tag_cb_corrected
        str_name_bam_tag_cb_uncorrected = args.str_name_bam_tag_cb_uncorrected
        str_name_bam_tag_umi_corrected = args.str_name_bam_tag_umi_corrected
        str_name_bam_tag_umi_uncorrected = args.str_name_bam_tag_umi_uncorrected
        flag_skip_read_analysis_summary_output_bam_file = (
            args.flag_skip_read_analysis_summary_output_bam_file
        )
        flag_skip_read_analysis_summary_output_tsv_file = (
            args.flag_skip_read_analysis_summary_output_tsv_file
        )
        flag_does_not_delete_sequence_and_sequence_qual = (
            args.flag_does_not_delete_sequence_and_sequence_qual
        )
        flag_does_not_collect_variant_information = (
            args.flag_does_not_collect_variant_information
        )
        float_min_prop_of_reads_for_filtering_genomic_variant = (
            args.float_min_prop_of_reads_for_filtering_genomic_variant
        )
        flag_turn_off_catching_all_reads_by_binning = (
            args.flag_turn_off_catching_all_reads_by_binning
        )
        int_bp_for_bins = args.int_bp_for_bins
        flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning = (
            args.flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning
        )
        path_file_fa_for_cram = args.path_file_fa_for_cram
        int_min_count_features_for_filtering_barcodes = (
            args.int_min_count_features_for_filtering_barcodes
        )
        path_file_vcf_for_filtering_variant = args.path_file_vcf_for_filtering_variant
        flag_output_variant_information_with_annotations = (
            args.flag_output_variant_information_with_annotations
        )
        dict_num_manager_processes_for_each_data_object = ast.literal_eval(
            args.dict_num_manager_processes_for_each_data_object
        )  # parse dictionary in a string
        int_num_samples_analyzed_concurrently = (
            args.int_num_samples_analyzed_concurrently
        )
        l_seqname_to_skip = args.l_seqname_to_skip
        flag_skip_intron_retention_counting = args.flag_skip_intron_retention_counting
        flag_no_strand_specificity = args.flag_no_strand_specificity
        int_min_length_intron_for_detecting_intron_retention_event = args.int_min_length_intron_for_detecting_intron_retention_event
        
        path_folder_reference_distribution = args.path_folder_reference_distribution
        l_name_distribution = args.l_name_distribution
        l_str_l_t_distribution_range_of_interest = args.l_str_l_t_distribution_range_of_interest
        
        int_length_of_polya_to_append_to_transcript_sequence_during_realignment = args.int_length_of_polya_to_append_to_transcript_sequence_during_realignment
        flag_enforce_transcript_end_site_matching_for_long_read_during_realignment = args.flag_enforce_transcript_end_site_matching_for_long_read_during_realignment
        int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment = args.int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment
        str_name_bam_tag_internal_polya_length = args.str_name_bam_tag_internal_polya_length
        str_mappy_aligner_preset_for_realignment = args.str_mappy_aligner_preset_for_realignment
        int_length_of_polya_to_append_to_transcript_sequence_during_realignment = args.int_length_of_polya_to_append_to_transcript_sequence_during_realignment
        flag_skip_internal_polyA_filtered_feature_counting = args.flag_skip_internal_polyA_filtered_feature_counting

    """
    Start of the Scarab-Count Program
    """
    logger.info(str_description)
    logger.info(
        "Scarab Short Read: secondary analysis of 10X barcoded BAM file for isoform and misc. annotation analysis"
    )
    logger.info(f"[Scarab-Count Start] Program Started.")

    """ handle special cases and invalid inputs """
    if l_path_file_bam_input is None or (
        path_folder_ref is None
        and (
            path_file_gtf_genome is None
            or path_file_fa_transcriptome is None
            or path_file_fa_genome is None
        )
    ):
        logger.error(
            "Required argument(s) is missing. to view help message, type -h or --help"
        )
        return -1
    set_valid_str_mode_scarab_count = {
        "gex5prime-single-end",
        "gex5prime-paired-end",
        "gex3prime-single-end",
        "gex3prime",
        "gex",
        "gex5prime",
        "atac",
        "multiome",
    }  # define a valid set of 'str_mode_scarab_count'
    if (
        sum(
            str_mode_scarab_count not in set_valid_str_mode_scarab_count
            for str_mode_scarab_count in l_str_mode_scarab_count
        )
        > 0
    ):
        logger.error(
            f"invalid 'str_mode_scarab_count' mode detected. Only one of {set_valid_str_mode_scarab_count} can be used. exiting"
        )
        return -1
    # replace short-hand to the full-length scarab_count mode
    dict_str_mode_scarab_count_short_hand = {
        "gex3prime": "gex3prime-single-end",
        "gex": "gex3prime-single-end",
        "gex5prime": "gex5prime-single-end",
    }  # define mapping
    l_str_mode_scarab_count = list(
        dict_str_mode_scarab_count_short_hand[e]
        if e in dict_str_mode_scarab_count_short_hand
        else e
        for e in l_str_mode_scarab_count
    )  # replace short-hand to the full-length scarab_count mode
    if (
        len(l_str_mode_scarab_count) == 1
    ):  # if only a single operating mode for scarab count has been given, apply (broadcast) the operating mode to all samples
        l_str_mode_scarab_count = list(l_str_mode_scarab_count) * len(
            l_path_folder_output
        )
    if len(l_path_folder_output) != len(l_str_mode_scarab_count):
        logger.error(
            "the number of samples are not equal to the number of given Scarab operating modes (the 'l_str_mode_scarab_count' argument), exiting"
        )
        return -1
    # check the number of input BAM files
    int_num_of_expected_input_bam_files = sum(
        2 if str_mode_scarab_count == "multiome" else 1
        for str_mode_scarab_count in l_str_mode_scarab_count
    )  # calculate the number of expected input BAM files based on the given 'l_str_mode_scarab_count' argument.
    if len(l_path_file_bam_input) != int_num_of_expected_input_bam_files:
        logger.error(
            f"the number of given input BAM files are {len( l_path_file_bam_input )}, which is different from the expected number of input BAM files (expected {int_num_of_expected_input_bam_files} files) from the given 'l_str_mode_scarab_count' argument, exiting"
        )
        return -1

    """ process required input directories """
    if path_file_fa_genome is not None:
        path_file_fa_genome = os.path.abspath(path_file_fa_genome)
    if path_file_gtf_genome is not None:
        path_file_gtf_genome = os.path.abspath(path_file_gtf_genome)
    if path_file_fa_transcriptome is not None:
        path_file_fa_transcriptome = os.path.abspath(path_file_fa_transcriptome)
    if path_file_fa_for_cram is not None:
        path_file_fa_for_cram = os.path.abspath(path_file_fa_for_cram)

    """ process input directoy  """
    l_path_file_bam_input = list(
        os.path.abspath(path_file_bam_input)
        for path_file_bam_input in l_path_file_bam_input
    )
    if l_path_folder_output is not None:
        """# when a valid list of output folders were given # ensure directories of the output folder ends with '/' characters"""
        l_path_folder_output = list(
            os.path.abspath(path_folder) + "/" for path_folder in l_path_folder_output
        )
    else:
        """# compose a list of default 'path_folder_output' values for the given list of input BAM files"""
        l_path_file_bam_input_reversed = deepcopy(
            l_path_file_bam_input[::-1]
        )  # reverse the input BAM file paths so that pop operation yield the element located at the front
        l_path_folder_output = []
        for str_mode_scarab_count in l_str_mode_scarab_count:
            path_file_bam_input = l_path_file_bam_input_reversed.pop()
            if (
                str_mode_scarab_count == "multiome"
            ):  # consume one more input BAM file path (a file path for ATAC-seq data) for a 'multiome' sample
                _ = l_path_file_bam_input_reversed.pop()
            path_folder_output = (
                f"{path_file_bam_input.rsplit( '/', 1 )[ 0 ]}scarab_output/"
            )
            l_path_folder_output.append(path_folder_output)

    """ set default reference annotation folder (a subdirectory inside the output folder of the first BAM file input) """
    if path_folder_ref is None:
        path_folder_ref = l_path_folder_output[0].rsplit("/", 1)[0] + "ref/"
    path_folder_ref = os.path.abspath(path_folder_ref)
    path_folder_ref += "/"

    # parse optional arguments receiving directories
    if path_file_tsv_repeatmasker_ucsc is not None:
        path_file_tsv_repeatmasker_ucsc = os.path.abspath(
            path_file_tsv_repeatmasker_ucsc
        )
    if path_file_gff_regulatory_element is not None:
        path_file_gff_regulatory_element = os.path.abspath(
            path_file_gff_regulatory_element
        )
        
    ''' pre-process input arguments '''
    l_l_t_distribution_range_of_interest = None # initialize 'l_l_t_distribution_range_of_interest'
    def _t_distribution_range_of_interest_to_str( t_distribution_range_of_interest ) :
        ''' # 2023-08-29 23:58:21 
        function for converting 't_distribution_range_of_interest' to a string
        '''
        return 'raw' if t_distribution_range_of_interest is None else '-'.join( list( map( str, t_distribution_range_of_interest ) ) )
    if l_str_l_t_distribution_range_of_interest is not None : # if 'l_l_t_distribution_range_of_interest' has been given
        # pre-process 'l_str_l_t_distribution_range_of_interest' - fill missing values
        if isinstance( l_str_l_t_distribution_range_of_interest, str ) : # if a string was given, wrap the string in a list
            l_str_l_t_distribution_range_of_interest = [ l_str_l_t_distribution_range_of_interest ]
        if len( l_str_l_t_distribution_range_of_interest ) == 1 and len( l_path_folder_output ) > 1 : # if 'l_str_l_t_distribution_range_of_interest' only contain a single entry and more than one sample has been given, use the string for all the samples
            l_str_l_t_distribution_range_of_interest *= len( l_path_folder_output )

        # parse 'l_str_l_t_distribution_range_of_interest'
        def _process_str_distribution_range_of_interest( e ) :
            e = e.strip( )
            return None if e.lower( ) in { 'raw', 'raw_count', 'rawcount' } else tuple( map( int, e.split( '-' ) ) )
        l_l_t_distribution_range_of_interest = list( list( set( _process_str_distribution_range_of_interest( e ) for e in set( e.strip( ).split( ',' ) ) ) ) for e in l_str_l_t_distribution_range_of_interest ) # parse each 'str_l_t_distribution_range_of_interest'
    
    """ 
    Fixed Settings
    """
    # settings for compatibility
    int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error = (
        75  # size of typically small exon that causes splicing site detection error
    )

    # internal settings
    (
        flag_use_gene_assignment_from_10x_cellranger,
        flag_use_isoform_assignment_from_10x_cellranger,
        flag_use_intronic_read_assignment_from_10x_cellranger,
    ) = (False, False, False)
    if flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger:
        (
            flag_use_gene_assignment_from_10x_cellranger,
            flag_use_isoform_assignment_from_10x_cellranger,
            flag_use_intronic_read_assignment_from_10x_cellranger,
        ) = (True, True, True)
    flag_filtering_alignment_to_transcript_during_realignment_based_on_structural_difference = int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment >= 0 # retrieve a flag indicating whether to filter alignments to transcripts based on the soft-clipping status

    # process arguments
    set_seqname_to_skip = set(l_seqname_to_skip)

    """
    read reference distribution
    """
    dict_name_sample_to_arr_ratio_to_ref = None # initialize 'dict_name_sample_to_arr_ratio_to_ref'
    if path_folder_reference_distribution is not None : # if 'path_folder_reference_distribution' has been given
        path_folder_reference_distribution = os.path.abspath( os.path.realpath( path_folder_reference_distribution ) ) + '/' # preprocess the path
        dict_output = bk.PICKLE_Read( f"{path_folder_reference_distribution}dict_output.pickle" ) # read reference distribution data
        dict_name_sample_to_arr_ratio_to_ref = dict( ( name_sample, arr_ratio_to_ref ) for name_sample, arr_ratio_to_ref in zip( dict_output[ 'setting' ][ 'l_name_file_distributions' ], dict_output[ 'l_arr_ratio_to_ref' ] ) ) # map name_sample to correction ratio data
        del dict_output
    # check whether reference distribution is used for normalization (check all the required arguments were given)
    flag_size_distribution_based_normalization_is_applied = not ( dict_name_sample_to_arr_ratio_to_ref is None or l_l_t_distribution_range_of_interest is None or l_name_distribution is None or len( l_name_distribution ) == 0 ) # retrieve a flag indicating whether size-distribution-based normalization will be applied
    if not flag_size_distribution_based_normalization_is_applied : # if size distribution normalization is not applied, put some dummy values into the list
        l_l_t_distribution_range_of_interest = [ [ None ] ] * len( l_path_folder_output ) # [ None ] representing exporting count matrix without size-distribution-based normalization
        l_name_distribution = [ None ] * len( l_path_folder_output )
    
    """
    Preprocess gene, isoform, and miscellaneous annotations
    """
    # load annotations
    if scidx is None:
        scidx = Preprocess_and_Load_Annotations(
            path_folder_ref,
            path_file_gtf_genome,
            str_name_gtf_attr_for_id_gene,
            str_name_gtf_attr_for_name_gene,
            str_name_gtf_attr_for_id_transcript,
            str_name_gtf_attr_for_name_transcript,
            flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation,
            flag_does_not_make_gene_names_unique,
            int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
            path_file_tsv_repeatmasker_ucsc,
            l_repClass_repeatmasker_ucsc,
            int_min_length_repeatmasker_ucsc,
            path_file_gff_regulatory_element,
            str_name_gff_attr_id_regulatory_element,
            int_min_length_regulatory_element,
            int_bp_padding_regulatory_element_anno,
            path_file_fa_transcriptome,
            path_file_fa_genome,
            int_bp_padding_for_defining_promoter_from_transcript_start,
            flag_does_not_remove_the_version_information_from_id_transcript_in_the_file_fa_transcriptome,
            dict_num_manager_processes_for_each_data_object=dict_num_manager_processes_for_each_data_object,
        )
    # return scidx
    """
    Preprocess variant annotations
    """
    # retrieve a flag indicating how variant filtering
    flag_filter_variant_using_predefined_set = (
        path_file_vcf_for_filtering_variant is not None
    )
    if (
        flag_filter_variant_using_predefined_set
    ):  # read predefined set of variants for filtering variants
        df_vcf = pd.read_csv(
            path_file_vcf_for_filtering_variant,
            comment="#",
            skip_blank_lines=True,
            sep="\t",
            header=None,
            usecols=[0, 1, 3, 4],
        )  # load vcf file, loading only required columns
        df_vcf.columns = ["CHROM", "POS", "REF", "ALT"]
        df_vcf["CHROM"] = (
            df_vcf.CHROM.astype(str)
            .astype(object)
            .apply(__chromosome_name_remove_chr__)
        )  # process chromosome name
        set_var_name_valid = set(
            f"{chrom}:{pos}:{ref}>{alt}" for chrom, pos, ref, alt in df_vcf.values
        )  # retrieve set of variant names for filtering
        """ build an interval tree of genomic positions of a predefined set of variants """
        dict_it_pos_variant_predefined_set = (
            dict()
        )  # a dictionary of interval tree for identifying reads overlapping with a predefined set of variants
        for chrom, pos, ref, alt in df_vcf.values:  # 1-based
            if (
                chrom not in dict_it_pos_variant_predefined_set
            ):  # initialize the interval tree
                dict_it_pos_variant_predefined_set[chrom] = intervaltree.IntervalTree()
            dict_it_pos_variant_predefined_set[chrom][pos - 1 : pos] = (
                chrom,
                pos,
            )  # update interval tree for searching variant positions

    """
    Exit early when no samples is anlayzed
    """
    # if no samples will be analyzed, return
    if len(l_path_folder_output) == 0:
        logger.info(f"[Scarab-Count Start] Program Completed.")
        return scidx  # return the loaded index object

    """
    Initiate pipelines for off-loading works
    """
    pipelines = bk.Offload_Works(
        None
    )  # no limit for the number of works that can be submitted.

    int_num_samples_analyzed_concurrently = min(
        len(l_path_folder_output), int_num_samples_analyzed_concurrently
    )  # if the number of samples are smaller than 'int_num_samples_analyzed_concurrently', adjust 'int_num_samples_analyzed_concurrently' so that it matches the number of samples

    n_threads = int(
        np.ceil(n_threads / int_num_samples_analyzed_concurrently)
    )  # divide the number of processes that can be used by each pipeline by the number of pipelines that will be run concurrently.

    def run_pipeline():
        """# 2023-07-27 12:09:01 
        analyze a pipeline for a given list of samples
        """
        # retrieve id of the pipeline
        str_uuid_pipeline = bk.UUID()
        logger.info(
            f"[Pipeline Start] Forked Pipeline (id={str_uuid_pipeline}) Started."
        )

        """
        Initiate workers for off-loading works
        """
        workers = bk.Offload_Works(
            None
        )  # no limit for the number of works that can be submitted.

        """
        Run scarab count for each sample
        """
        def _reverse_list( l ) :
            ''' # 2023-08-29 22:49:13 
            copy and reverse the list
            '''
            return deepcopy( l[::-1] )
        l_path_file_bam_input_reversed = _reverse_list( l_path_file_bam_input ) # reverse the input BAM file paths so that pop operation yield the element located at the front
        l_l_t_distribution_range_of_interest_reversed = _reverse_list( l_l_t_distribution_range_of_interest )
        for str_mode_scarab_count_for_the_current_sample, path_folder_output_for_the_current_sample, name_distribution in zip( l_str_mode_scarab_count, l_path_folder_output, l_name_distribution ) :  # retrieve scarab count operating mode and an output folder for the current sample
            """settings for each scarab count operating mode"""
            if (
                str_mode_scarab_count_for_the_current_sample == "multiome"
            ):  # Multiome run mode
                l_path_file_bam_input_for_the_current_sample = [
                    l_path_file_bam_input_reversed.pop(),
                    l_path_file_bam_input_reversed.pop(),
                ]  # retrieve two samples for the multiome sample
                l_str_mode_scarab_count_for_the_current_sample = [
                    "gex3prime-single-end",
                    "atac",
                ]
                l_path_folder_output_for_the_current_sample = [
                    f"{path_folder_output_for_the_current_sample}gex/",
                    f"{path_folder_output_for_the_current_sample}atac/",
                ]
                l_int_min_mapq_unique_mapped_for_the_current_sample = [
                    int_min_mapq_unique_mapped_for_gex_data,
                    int_min_mapq_unique_mapped_for_atac_data,
                ]
                l_l_t_distribution_range_of_interest_for_the_current_sample = [
                    l_l_t_distribution_range_of_interest_reversed.pop(),
                    l_l_t_distribution_range_of_interest_reversed.pop(),
                ]
                raise NotImplementedError( 'atac distribution not implemented.' )
                l_arr_ratio_to_ref_for_the_current_sample = [
                    None, None
                ]
            elif (
                str_mode_scarab_count_for_the_current_sample[:4] == "atac"
            ):  # ATAC run mode
                l_path_file_bam_input_for_the_current_sample = [
                    l_path_file_bam_input_reversed.pop()
                ]  # retrieve one sample for the atac sample
                l_str_mode_scarab_count_for_the_current_sample = [
                    str_mode_scarab_count_for_the_current_sample
                ]
                l_path_folder_output_for_the_current_sample = [
                    path_folder_output_for_the_current_sample
                ]
                l_int_min_mapq_unique_mapped_for_the_current_sample = [
                    int_min_mapq_unique_mapped_for_atac_data
                ]
                l_l_t_distribution_range_of_interest_for_the_current_sample = [
                    l_l_t_distribution_range_of_interest_reversed.pop(),
                ]
                l_arr_ratio_to_ref_for_the_current_sample = [ dict_name_sample_to_arr_ratio_to_ref[ name_distribution ], ] if flag_size_distribution_based_normalization_is_applied else [ None, ]
            elif (
                str_mode_scarab_count_for_the_current_sample[:3] == "gex"
            ):  # GEX run mode
                l_path_file_bam_input_for_the_current_sample = [
                    l_path_file_bam_input_reversed.pop()
                ]  # retrieve one sample for the gex sample
                l_str_mode_scarab_count_for_the_current_sample = [
                    str_mode_scarab_count_for_the_current_sample
                ]
                l_path_folder_output_for_the_current_sample = [
                    path_folder_output_for_the_current_sample
                ]
                l_int_min_mapq_unique_mapped_for_the_current_sample = [
                    int_min_mapq_unique_mapped_for_gex_data
                ]
                l_l_t_distribution_range_of_interest_for_the_current_sample = [
                    l_l_t_distribution_range_of_interest_reversed.pop(),
                ]
                l_arr_ratio_to_ref_for_the_current_sample = [ dict_name_sample_to_arr_ratio_to_ref[ name_distribution ], ] if flag_size_distribution_based_normalization_is_applied else [ None, ]

            """
            define a function to release a lock
            """

            def release_lock():
                """# 2023-01-14 20:36:17
                release the lock file
                """
                path_file_lock = (
                    f"{path_folder_output_for_the_current_sample}ourotools.lock"
                )

                # check the existence of output files for the output folder of each BAM file of the current sample
                flag_all_output_files_exist = True  # initialize the flag
                for path_folder_output in l_path_folder_output_for_the_current_sample:
                    if not os.path.exists(
                        f"{path_folder_output}count_matrix.export_completed.txt"
                    ):
                        flag_all_output_files_exist = False
                        break
                    if not flag_skip_read_analysis_summary_output_bam_file:
                        if not os.path.exists(
                            f"{path_folder_output}bam_output.export_completed.txt"
                        ):
                            flag_all_output_files_exist = False
                            break

                # check the existence of the output folder for the current sample
                if str_mode_scarab_count_for_the_current_sample == "multiome":
                    if not os.path.exists(
                        f"{path_folder_output_for_the_current_sample}count_matrix.export_completed.txt"
                    ):
                        flag_all_output_files_exist = False

                # check the existence of the lock file
                if (
                    os.path.exists(path_file_lock) and flag_all_output_files_exist
                ):  # if all output files exist and the lock file exists
                    # check whether the lock file has been created by the current pipeline
                    with open(path_file_lock, "rt") as file_lock:
                        flag_lock_acquired = file_lock.read() == str_uuid_pipeline
                    if (
                        flag_lock_acquired
                    ):  # if the lock file has been created by the current pipeline, delete the lock file
                        os.remove(path_file_lock)
                    # lock has been released
                    if verbose:
                        logger.warning(
                            f"[{path_folder_output_for_the_current_sample}] The forked pipeline (id={str_uuid_pipeline}) released the lock"
                        )
                else:
                    if verbose:
                        logger.warning(
                            f"[{path_folder_output_for_the_current_sample}] The forked pipeline (id={str_uuid_pipeline}) attempted to release the lock, but some output files are missing, and the lock will not be released."
                        )

            """
            Run pipeline for each sample
            """
            """
            create a lock
            """
            os.makedirs(path_folder_output_for_the_current_sample, exist_ok=True)
            path_file_lock = (
                f"{path_folder_output_for_the_current_sample}ourotools.lock"
            )
            # check the existence of the lock file
            if os.path.exists(path_file_lock):
                logger.warning(
                    f"[Output folder unavailable] the output folder {path_folder_output_for_the_current_sample} contains a lock file, which appears to be processed by a different process. Therefore, the output folder will be skipped."
                )
                continue
            flag_lock_acquired = False  # initialize 'flag_lock_acquired'
            try:
                # create the lock file
                with open(path_file_lock, "wt") as newfile_lock:
                    newfile_lock.write(str_uuid_pipeline)
                # check whether the lock file has been created correctly (check for collision).
                with open(path_file_lock, "rt") as file_lock:
                    flag_lock_acquired = file_lock.read() == str_uuid_pipeline
            except Exception as e:
                logger.critical(
                    e, exc_info=True
                )  # if an exception occurs, print the error message
            if not flag_lock_acquired:
                logger.warning(
                    f"[Output folder unavailable] an attempt to acquire a lock for the output folder {path_folder_output_for_the_current_sample} failed, which appears to be processed by a different process. Therefore, the output folder will be skipped."
                )
                continue
            # lock has been acquired

            """
            Run pipeline for each BAM file
            """
            for (
                path_file_bam_input,
                str_mode_scarab_count,
                path_folder_output,
                int_min_mapq_unique_mapped,
                arr_ratio_to_ref,
                l_t_distribution_range_of_interest,
            ) in zip(
                l_path_file_bam_input_for_the_current_sample,
                l_str_mode_scarab_count_for_the_current_sample,
                l_path_folder_output_for_the_current_sample,
                l_int_min_mapq_unique_mapped_for_the_current_sample,
                l_arr_ratio_to_ref_for_the_current_sample,
                l_l_t_distribution_range_of_interest_for_the_current_sample,
            ):
                # define folders and directories
                path_file_bam_input = os.path.abspath(path_file_bam_input)
                if path_folder_output is None:  # set default 'path_folder_output'
                    path_folder_output = (
                        f"{path_file_bam_input.rsplit( '/', 1 )[ 0 ]}scarab_output/"
                    )
                path_folder_output = os.path.abspath(path_folder_output)
                path_folder_output += "/"
                path_folder_temp = f"{path_folder_output}temp/"
                path_folder_graph = f"{path_folder_output}graph/"

                """ if the output folder already exists """
                if os.path.exists(path_folder_output):
                    """check whether the pipeline has been completed"""
                    ''' # 2023-09-07 12:41:22 
                    check whether valid count matrices has been exported
                    '''
                    if os.path.exists( f"{path_folder_output}count_matrix.export_completed.txt" ) : # check the presence of count matrix, which should be present if a typical run has been completed.
                        logger.info( f"[Output folder Already Exists] the output folder {path_folder_output} contains a valid count matrix file. Therefore, the output folder will be skipped." )
                        continue  # skip if the pipeline has been completed for the output folder
                    else:
                        """if required output files does not exist or the an intermediate file exists, remove the entire output folder, and rerun Scarab"""
                        if ( len(glob.glob(f"{path_folder_output}*/")) > 0 ):  # detect a folder inside the output folder and report the presence of the existing folders.
                            logger.info( f"[Output folder Already Exists] the output folder {path_folder_output} does not contain a valid count matrix file. The output folder will be cleaned and the pipeline will start anew." )
                        # delete the folders inside the output folder
                        for path_folder in glob.glob(f"{path_folder_output}*/"):
                            shutil.rmtree(path_folder)
                        # delete the files, excluding the lock file
                        for path_file in glob.glob(f"{path_folder_output}*"):
                            if ( path_file_lock != path_file ):  # does not delete the lock file
                                os.remove(path_file)

                """ create directories """
                for path_folder in [
                    path_folder_output,
                    path_folder_temp,
                    path_folder_graph,
                    path_folder_ref,
                ]:
                    os.makedirs(path_folder, exist_ok=True)

                """ index input bam files if index files do not exist """
                flag_corrupted_bam_file = False
                for path_file_bam in [path_file_bam_input]:
                    path_file_index = (
                        f"{path_file_bam}.bai"
                        if path_file_bam.rsplit(".", 1)[1].lower() == "bam"
                        else f"{path_file_bam}.crai"
                    )  # retrieve path of the index of the input BAM (or CRAM) file
                    if not os.path.exists(path_file_index):
                        try:
                            pysam.index(path_file_bam)
                        except:
                            flag_corrupted_bam_file = True
                if flag_corrupted_bam_file:
                    logger.info(
                        f"[Corrupted input BAM file] the input BAM file {path_file_bam} cannot be indexed. Therefore, the input file will be skipped."
                    )
                    continue

                """
                Report program arguments
                """
                # record arguments used for the program (metadata)
                dict_program_setting = {
                    "version": _version_,  # record version
                    # external
                    "path_file_bam_input": path_file_bam_input,
                    "path_file_fa_genome": path_file_fa_genome,
                    "path_file_gtf_genome": path_file_gtf_genome,
                    "path_file_fa_transcriptome": path_file_fa_transcriptome,
                    "path_folder_output": path_folder_output,
                    "path_folder_ref": path_folder_ref,
                    "n_threads": n_threads,
                    "float_memory_in_GiB": float_memory_in_GiB,
                    "int_min_mapq_unique_mapped": int_min_mapq_unique_mapped,
                    "int_n_bases_padding_around_interval": int_n_bases_padding_around_interval,
                    "path_file_tsv_repeatmasker_ucsc": path_file_tsv_repeatmasker_ucsc,
                    "l_repClass_repeatmasker_ucsc": l_repClass_repeatmasker_ucsc,
                    "int_min_length_repeatmasker_ucsc": int_min_length_repeatmasker_ucsc,
                    "flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger": flag_use_gene_isoform_and_intron_assignment_from_10x_cellranger,
                    "flag_include_read_aligned_to_opposite_strand": flag_include_read_aligned_to_opposite_strand,
                    "flag_include_read_aligned_to_intron": flag_include_read_aligned_to_intron,
                    "str_name_gtf_attr_for_id_gene": str_name_gtf_attr_for_id_gene,
                    "str_name_gtf_attr_for_name_gene": str_name_gtf_attr_for_name_gene,
                    "str_name_gtf_attr_for_id_transcript": str_name_gtf_attr_for_id_transcript,
                    "str_name_gtf_attr_for_name_transcript": str_name_gtf_attr_for_name_transcript,
                    "path_file_gff_regulatory_element": path_file_gff_regulatory_element,
                    "str_name_gff_attr_id_regulatory_element": str_name_gff_attr_id_regulatory_element,
                    "int_min_length_regulatory_element": int_min_length_regulatory_element,
                    "int_bp_padding_regulatory_element_anno": int_bp_padding_regulatory_element_anno,
                    "float_max_prop_unfiltered_rpmk": float_max_prop_unfiltered_rpmk,
                    "flag_does_not_make_gene_names_unique": flag_does_not_make_gene_names_unique,
                    "flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation": flag_does_not_merge_overlapping_genes_with_the_same_gene_name_and_strand_orientation,
                    "str_name_bam_tag_cb_corrected": str_name_bam_tag_cb_corrected,
                    "str_name_bam_tag_cb_uncorrected": str_name_bam_tag_cb_uncorrected,
                    "str_name_bam_tag_umi_corrected": str_name_bam_tag_umi_corrected,
                    "str_name_bam_tag_umi_uncorrected": str_name_bam_tag_umi_uncorrected,
                    "flag_does_not_delete_sequence_and_sequence_qual": flag_does_not_delete_sequence_and_sequence_qual,
                    "flag_does_not_collect_variant_information": flag_does_not_collect_variant_information,
                    "int_min_num_of_reads_for_filtering_genomic_variant": int_min_num_of_reads_for_filtering_genomic_variant,
                    "float_min_prop_of_reads_for_filtering_genomic_variant": float_min_prop_of_reads_for_filtering_genomic_variant,
                    "flag_turn_off_catching_all_reads_by_binning": flag_turn_off_catching_all_reads_by_binning,
                    "int_bp_for_bins": int_bp_for_bins,
                    "flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning": flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning,
                    "int_bp_padding_for_defining_promoter_from_transcript_start": int_bp_padding_for_defining_promoter_from_transcript_start,
                    "flag_skip_exon_and_splice_junc_counting": flag_skip_exon_and_splice_junc_counting,
                    'flag_skip_intron_retention_counting' : flag_skip_intron_retention_counting,
                    'flag_no_strand_specificity' : flag_no_strand_specificity,
                    "path_file_fa_for_cram": path_file_fa_for_cram,
                    "flag_output_variant_information_with_annotations": flag_output_variant_information_with_annotations,
                    "path_file_vcf_for_filtering_variant": path_file_vcf_for_filtering_variant,
                    "int_min_count_features_for_filtering_barcodes": int_min_count_features_for_filtering_barcodes,
                    "dict_num_manager_processes_for_each_data_object": dict_num_manager_processes_for_each_data_object,
                    "l_seqname_to_skip": l_seqname_to_skip,
                    "flag_skip_read_analysis_summary_output_bam_file": flag_skip_read_analysis_summary_output_bam_file,
                    "flag_skip_read_analysis_summary_output_tsv_file": flag_skip_read_analysis_summary_output_tsv_file,
                    'path_folder_reference_distribution' : path_folder_reference_distribution,
                    'l_name_distribution' : l_name_distribution,
                    'l_str_l_t_distribution_range_of_interest' : l_str_l_t_distribution_range_of_interest,
                    'int_length_of_polya_to_append_to_transcript_sequence_during_realignment' : int_length_of_polya_to_append_to_transcript_sequence_during_realignment,
                    'flag_enforce_transcript_end_site_matching_for_long_read_during_realignment' : flag_enforce_transcript_end_site_matching_for_long_read_during_realignment,
                    'int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment' : int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment,
                    'str_name_bam_tag_internal_polya_length' : str_name_bam_tag_internal_polya_length,
                    'flag_skip_internal_polyA_filtered_feature_counting' : flag_skip_internal_polyA_filtered_feature_counting,
                    # internal
                    "path_folder_temp": path_folder_temp,
                    "path_folder_graph": path_folder_graph,
                    "path_folder_ref": path_folder_ref,
                    "int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error": int_window_size_from_the_end_of_exon_to_ignore_mutation_calling_from_splice_site_detection_error,
                    "flag_use_gene_assignment_from_10x_cellranger": flag_use_gene_assignment_from_10x_cellranger,
                    "flag_use_isoform_assignment_from_10x_cellranger": flag_use_isoform_assignment_from_10x_cellranger,
                    "flag_use_intronic_read_assignment_from_10x_cellranger": flag_use_intronic_read_assignment_from_10x_cellranger,
                }
                logger.info(
                    f"[Setting] program will be run with the following setting for the input BAM file {path_file_bam_input} : {str( dict_program_setting )}"
                )

                """ export program setting """
                path_file_json_setting_program = (
                    f"{path_folder_output}program_setting.json"
                )
                if os.path.exists(path_file_json_setting_program):
                    with open(path_file_json_setting_program, "r") as file:
                        j = json.load(file)
                    if j != dict_program_setting:
                        logger.info(
                            f"[Warning] the current program setting is different from the previous program setting recorded in the pipeline folder. The previous setting will be used."
                        )
                        with open(path_file_json_setting_program, "r") as file:
                            dict_program_setting = json.load(
                                file
                            )  # override current program setting with previous program setting
                with open(path_file_json_setting_program, "w") as newfile:
                    json.dump(dict_program_setting, newfile)

                # fixed setting
                name_ref = "genome"
                flag_ignore_record_with_flag_secondary_alignment = True
                flag_ignore_record_with_flag_optical_or_pcr_duplicate = True
                flag_ignore_reads_with_ambiguous_gene_assignment = True

                # define interger representation of the CIGAR operations used in BAM files
                int_cigarop_I = 1
                int_cigarop_D = 2
                int_cigarop_N = 3
                int_cigarop_S = 4
                int_cigarop_H = 5
                
                set_int_cigarop_structural_difference = { int_cigarop_I, int_cigarop_D, int_cigarop_N, int_cigarop_S }
                def _mappy_detect_structural_difference( mappy_cigartuples ) :
                    """ # 2023-09-21 23:06:25 
                    from the cigartuples of the mappy aligned segment record, detect structural difference.
                    return False if a structural difference is present, and return True if no structural difference is present, based on the given cigartuples
                    """
                    for length_of_operation, operation in mappy_cigartuples :
                        if operation in set_int_cigarop_structural_difference and length_of_operation > int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment :
                            return False
                    return True

                """ load setting for reference annotations """
                path_file_json_ref_setting = f"{path_folder_ref}ref_setting.json"
                with open(path_file_json_ref_setting, "r") as file:
                    dict_setting_ref = json.load(
                        file
                    )  # override current program setting with previous program setting
                # parse settings update values in the local scope
                dict_seqname_to_len_seq = dict_setting_ref["dict_seqname_to_len_seq"]

                """ initialize base count array """
                dict_base_to_index = dict(
                    (str_base, index) for index, str_base in enumerate("ATGC")
                )
                int_n_base_types = len(
                    dict_base_to_index
                )  # retrieve the number of base types
                dict_seqname_to_index_chunk_to_base_count = dict(
                    (seqname, dict()) for seqname in dict_seqname_to_len_seq
                )

                """ load local settings """
                flag_is_mode_scarab_count_atac = str_mode_scarab_count == "atac"
                flag_is_5prime, flag_is_paired_end = None, None  # initialize the fiags
                # retrieve a flag indiciating whether to export filter excluding internal polyA primed reads
                flag_export_internal_polyA_filtered_count_as_a_separate_feature = not ( flag_is_mode_scarab_count_atac or flag_skip_internal_polyA_filtered_feature_counting )
                if not flag_is_mode_scarab_count_atac:
                    flag_is_5prime = "gex5prime" in str_mode_scarab_count
                    flag_is_paired_end = "paired-end" in str_mode_scarab_count

                """ define internal parameters according to the scarab count modes """
                flag_use_gene_assignment_from_10x_cellranger_for_the_current_bam_file = flag_use_gene_assignment_from_10x_cellranger  # retrieve 'flag_use_gene_assignment_from_10x_cellranger_for_the_current_bam_file'
                # for a header line for annotated count matrix
                l_col_df_count = ["barcode", "feature", "id_feature", "read_count"] # define 'l_col_df_count'
                if flag_is_mode_scarab_count_atac:
                    """data columns"""
                    # a list of names of columns for defining unique molecules assigned to a single feature for ATAC reads (following 10X cellranger ATAC > v1.2)
                    l_col_for_identifying_unique_molecules = [
                        "CB",
                        "refstart",
                        "refend",
                    ]  # first element should be 'CB'
                    # a list of names of columns for collecting information about the current alignment (full record)
                    # 'str_l_seg' is required for variant calling
                    l_col = [
                        "qname",
                        "mapq",
                        "refname",
                        "flag",
                        "refstart",
                        "refend",
                        "str_l_seg",
                        "CB",
                        "CR",
                        "int_flag_classification",
                        "id_rpmk",
                        "id_gene",
                        "id_reg",
                        "id_promoter",
                        "l_name_variant",
                        'int_total_aligned_length',
                    ]
                    # a list of names of columns for collecting information about the current alignment just for counting (only essential information to reduce memory footprint)
                    l_col_for_counting = [
                        "qname",
                        "mapq",
                        "flag",
                        "str_l_seg",
                        "CB",
                        "refstart",
                        "refend",
                        "int_flag_classification",
                        "l_name_variant",
                        'int_total_aligned_length',
                    ]
                    # a list of names of columns that will be added as BAM tags to the output BAM file # should be a contiguous subset of 'l_col'
                    l_name_col_newanno = l_col[ - 7 : ]

                    flag_use_gene_assignment_from_10x_cellranger_for_the_current_bam_file = (
                        False  # no gene assignment for ATAC data
                    )
                else:
                    # a list of names of columns for defining unique molecules assigned to a single feature for GEX reads (following 10X cellranger  > v1.2)
                    l_col_for_identifying_unique_molecules = [
                        "CB",
                        "UB",
                    ]  # first element should be 'CB'
                    # a list of names of columns for collecting information about the current alignment (full record)
                    # 'str_l_seg' is required for variant calling
                    l_col = [
                        "qname",
                        "mapq",
                        "refname",
                        "flag",
                        "refstart",
                        "refend",
                        "str_l_seg",
                        "CB",
                        "CR",
                        "UB",
                        "UR",
                        "TX",
                        "AN",
                        "GX",
                        "GN",
                        "MM",
                        "RE",
                        "xf",
                        "int_flag_classification",
                        "id_rpmk",
                        "int_max_num_base_pairs_overlap_with_rpmk",
                        "id_gene",
                        "int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id",
                        "int_base_gene_exon_count",
                        "int_base_filtered_rpmk_count",
                        "id_reg",
                        "int_base_unfiltered_rpmk_count",
                        "int_base_reg_count",
                        "id_tx_assigned_by_minimap2",
                        "l_name_variant",
                        'int_total_aligned_length',
                        'flag_internal_polya_tract_primed',
                    ]
                    # a list of names of columns for collecting information about the current alignment just for counting (only essential information to reduce memory footprint)
                    l_col_for_counting = [
                        "qname",
                        "mapq",
                        "flag",
                        "str_l_seg",
                        "CB",
                        "UB",
                        "TX",
                        "RE",
                        "int_flag_classification",
                        "id_tx_assigned_by_minimap2",
                        "l_name_variant",
                        'int_total_aligned_length',
                        'flag_internal_polya_tract_primed',
                    ]
                    # a list of names of columns that will be added as BAM tags to the output BAM file # should be a contiguous subset of 'l_col'
                    l_name_col_newanno = l_col[ -14 : ]

                # shared settings across scarab modes
                dict_name_col_newanno_to_sam_tag_name = {
                    "int_flag_classification": ("XC", "i"),
                    "id_rpmk": ("XR", "Z"),
                    "int_max_num_base_pairs_overlap_with_rpmk": ("YR", "i"),
                    "id_gene": ("XG", "Z"),
                    "int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id": (
                        "YG",
                        "i",
                    ),
                    "id_promoter" : ( 'XP', "Z" ),
                    "int_base_gene_exon_count": ("YX", "i"),
                    "int_base_filtered_rpmk_count": ("YF", "i"),
                    "id_reg": ("XE", "Z"),
                    "int_base_unfiltered_rpmk_count": ("YU", "i"),
                    "int_base_reg_count": ("YE", "i"),
                    "id_tx_assigned_by_minimap2": ("XT", "Z"),
                    "l_name_variant": ("XV", "Z"),
                    'int_total_aligned_length' : ( 'YA', 'i' ),
                    'flag_internal_polya_tract_primed' : ( 'ZI', 'i' ),
                }
                l_index_col_for_counting = list(
                    l_col.index(col) for col in l_col_for_counting
                )  # retrieve a list of column indices used for counting

                """
                int_flag_classification : binary flags

                gene        | 0x1 : a flag indicating an overlap(s) with gene body 
                gene        | 0x2 : a flag indicating that the read can be assigned to more than one genes (gene assignment was ambiguous)
                gene        | 0x4 : a flag indicating completely intronic reads (GEX mode specific)
                gene        | 0x8 : a flag indicating mostly exonic reads (>90% overlaps with exons) (GEX mode specific)
                promoter    | 0x10 : a flag indicating an overlap(s) with promoter regions (ATAC mode specific)
                promoter    | 0x20 : a flag indicating that the read can be assigned to more than one promoter regions (promoter assignment was ambiguous) (ATAC mode specific)
                repeats     | 0x40 : a flag indicating overlap with filtered repeatmasker annotations
                repeats     | 0x80 : a flag indicating that the read can be assigned to more than one filtered repeatmasker annotations (ambiguous assignment to filtered repeatmasker annotations)
                repeats     | 0x100 : a flag indicating complete overlap with filtered repeatmasek
                regulatory  | 0x200 : a flag indicating overlap with regulatory element(s)
                regulatory  | 0x400 : a flag indicating overlap with both repeatmasker annotations and regulatory element(s)
                regulatory  | 0x800 : a flag indicating overlap with regulatory element annotations and not overlapping with unfiltered repeatmasker annotations
                regulatory  | 0x1000 : a flag indicating that the read can be assigned to more than one regulatory elements (ambiguous assignment to regulatory annotations)
                regulatory  | 0x2000 : a flag indicating complete overlap with regulatory element annotations
                """
                int_flag_class_overlap_with_gene_body = 1 << 0
                int_flag_class_ambiguous_assignment_to_gene = 1 << 1
                int_flag_class_completely_intronic = 1 << 2
                int_flag_class_mostly_exonic = 1 << 3
                int_flag_class_overlap_with_promoter = 1 << 4
                int_flag_class_ambiguous_assignment_to_promoter = 1 << 5
                int_flag_class_overlap_with_filtered_rpmk_anno = 1 << 6
                int_flag_class_ambiguous_assignment_to_filtered_rpmk_anno = 1 << 7
                int_flag_class_complete_overlap_with_filtered_rpmk_anno = 1 << 8
                int_flag_class_overlap_with_reg = 1 << 9
                int_flag_class_overlap_with_reg_and_rpmk = 1 << 10
                int_flag_class_overlap_with_reg_not_overlap_with_unfiltered_rpmk_anno = (
                    1 << 11
                )
                int_flag_class_ambiguous_assignment_to_reg = 1 << 12
                int_flag_class_complete_overlap_with_reg_not_overlap_with_unfiltered_rpmk_anno = (
                    1 << 13
                )

                """
                Settings for BAM/CRAM files
                """
                flag_is_cram = (
                    path_file_bam_input.rsplit(".", 1)[0].lower() == "cram"
                )  # retrieve a flag indicating that the file is a CRAM file
                if flag_is_cram:
                    if path_file_fa_for_cram is None:  # show warning
                        if verbose:
                            logger.warning(
                                "a CRAM file has been given but reference genome path has not been given. A default location stored in CRAM file will be used."
                            )
                    else:
                        if not os.path.exists(
                            f"{path_file_fa_for_cram}.fai"
                        ):  # if the reference fasta file has not been indexed
                            bk.OS_Run(
                                ["samtools", "faidx", path_file_fa_for_cram]
                            )  # index the reference fasta file
                    if verbose:
                        logger.warning("path_file_fa_for_cram")
                str_mode_sam = "rc" if flag_is_cram else "rb"

                """ define internal functions and parameters """
                int_max_n_removed_elements = 10000  # maximum number of removed elements that can exist in a dictionary storing analyzed reads (avoid a 'memory leakeage')
                seq_polya_sequence_to_append_to_tx = 'A' * int_length_of_polya_to_append_to_transcript_sequence_during_realignment # retrieve poly A sequence to append
                int_max_distance_from_transcript_end_for_tes_matching_during_realignment = int_max_distance_from_transcript_start_and_end_for_tss_and_tes_matching_during_realignment + int_length_of_polya_to_append_to_transcript_sequence_during_realignment # retrieve 'int_max_distance_from_transcript_end_for_tes_matching_during_realignment', by adding the length of poly A tract sequence added to the end of transcript
                
                def __Get_Genomic_Region__(int_pos, int_bp_for_bins=int_bp_for_bins):
                    """get start and end coordinates (0-based) of the genomic region of the current position (0-based)"""
                    index_region = int(int_pos / int_bp_for_bins)
                    start = index_region * int_bp_for_bins  # 0-based
                    end = start + int_bp_for_bins  # 0-based
                    return start, end  # return 0-based coordinates

                def __get_data_object(name_data):
                    """get data object of the given 'name_data'"""
                    if name_data in scidx:
                        return scidx[name_data]
                    else:

                        def get_proxy_object():
                            return scidx[f"l_managed_{name_data}"][
                                _get_random_integer(
                                    dict_num_manager_processes_for_each_data_object[
                                        name_data
                                    ]
                                )
                            ]  # randomly select and return the proxy object to distribute the workload

                        return (
                            get_proxy_object  # return a function to select proxy object
                        )

                def __data_object_subset(data_object, l_key: list):
                    """# 2023-01-09 18:39:51
                    perform subset operation on a dictionary-like data object

                    l_key : list # list of key names to subset
                    """
                    if hasattr(data_object, "__call__"):  # for proxy object
                        return data_object().subset(l_key)
                    else:
                        dict_subset = dict()
                        for key in l_key:
                            if key in data_object:
                                dict_subset[key] = data_object[key]
                        return dict_subset  # return a subset dictionary of a given data object

                def __data_object_search_query(data_object, seqname, query):
                    """# 2023-01-09 18:39:47
                    perform subset operation on a dictionary-intervaltree-like data object
                    """
                    return (
                        data_object().search_query(seqname, query)
                        if hasattr(data_object, "__call__")
                        else (
                            data_object[seqname][query]
                            if seqname in data_object
                            else []
                        )
                    )

                def __data_object_search_queries(data_object, seqname, queries):
                    """# 2023-01-09 18:39:47
                    perform subset operation on a dictionary-like data object
                    """
                    return (
                        data_object().search_queries(seqname, queries)
                        if hasattr(data_object, "__call__")
                        else (
                            list(data_object[seqname][query] for query in queries)
                            if seqname in data_object
                            else []
                        )
                    )

                """
                Define a function for processing a part of a BAM file
                """
                def process_batch(pipe_receiver, pipe_sender):
                    """
                    # 2022-04-24 01:29:59
                    Requires loading several data objects (using copy-on-write method)

                    receives a bookmark file (either file directory of a tsv file or a dataframe)
                    """
                    str_uuid = bk.UUID()  # retrieve id

                    """ retrieve the sam header lines and open a new bam file using the input bam file as a template """
                    with pysam.AlignmentFile(
                        path_file_bam_input,
                        str_mode_sam,
                        reference_filename=path_file_fa_for_cram,
                    ) as samfile:
                        """retrieve sam header with chromosom names without 'chr' prefix"""
                        samfile_header = samfile.header
                        dict_header = samfile_header.to_dict()
                        """ remove 'chr' prefix from the header if it already exists  """
                        for e in dict_header["SQ"]:
                            e["SN"] = __chromosome_name_remove_chr__(e["SN"])
                        samfile_header = pysam.AlignmentHeader.from_dict(dict_header)
                        """ open a new sam file (output BAM file) """
                        if not flag_skip_read_analysis_summary_output_bam_file:
                            newsamfile = pysam.AlignmentFile(
                                f"{path_folder_temp}{str_uuid}.analysis.{name_ref}.bam",
                                "wb",
                                header=samfile_header,
                            )

                    """ open tsv files """
                    ''' open count matrix output files for each distribution of interest '''
                    dict_t_distribution_range_of_interest_to_newfile_df_count = dict( ( t_distribution_range_of_interest, gzip.open( f"{path_folder_temp}{str_uuid}.count.size_distribution__{_t_distribution_range_of_interest_to_str( t_distribution_range_of_interest )}.tsv.gz", "wb", ) ) for t_distribution_range_of_interest in l_t_distribution_range_of_interest )
                    # open other tsv output files
                    newfile_df_analysis_statistics = gzip.open(
                        f"{path_folder_temp}{str_uuid}.analysis_statistics.tsv.gz", "wb"
                    )  # for analyzing performance issues
                    if not flag_skip_read_analysis_summary_output_tsv_file:
                        newfile = gzip.open( f"{path_folder_temp}{str_uuid}.analysis.{name_ref}.tsv.gz", "wb", )

                    if verbose:
                        logger.info(f"[Started] ({str_uuid})")

                    """ get data objects (either real data object with direct access or a proxy object for accessing data object in the manager process) """
                    data_dict_it_exon_transcriptome = __get_data_object(
                        "dict_it_exon_transcriptome"
                    )
                    data_dict_it_rpmk = __get_data_object("dict_it_rpmk")
                    data_dict_it_splice_junc_transcriptome = __get_data_object(
                        "dict_it_splice_junc_transcriptome"
                    )
                    data_dict_it_exon = __get_data_object("dict_it_exon")
                    data_dict_it_splice_donor_and_acceptor = __get_data_object("dict_it_splice_donor_and_acceptor_genome")
                    data_dict_it_reg = __get_data_object("dict_it_reg")
                    data_dict_it_promoter = __get_data_object("dict_it_promoter")
                    data_dict_fa_transcriptome = __get_data_object(
                        "dict_fa_transcriptome"
                    )
                    data_dict_t_splice_junc_to_info_genome = __get_data_object(
                        "dict_t_splice_junc_to_info_genome"
                    )
                    
                    """
                    Initiate workers for off-loading works for processing each batch
                    """
                    int_num_workers_for_bucket_processing = 3 # the number of workers for bucket processing
                    workers_for_bucket_processing = bk.Offload_Works( int_num_workers_for_bucket_processing )  # 💛 adjustment of the number of workers might be needed.
  
                    int_max_bucket_deletion_count_before_reinitialize = 10000 # the max number of bucket deletion count before re-initializing the bucket container (python dictionary, when too many keys are deleted, lead to 'memory leak')
                    int_max_num_records_in_a_batch_of_buckets = 200000 # initialize the total number of records in a batch of buckets
                    int_max_num_batches_in_the_result_container_before_flushing = 2 # the max number of batches whose results can be stored in the container before flushing the result to the storage. if this number is too large, the process will consume too much memory 
                    
                    ns = { 'int_bucket_deletion_count' : 0 } # create a namespace for buckets # a counter counting the number of bucket deleted (for recreating the dictionary container to avoid 'memory leakage')
                                      
                    """
                    functions for initializing reads
                    """
                    def __Initialize_Reads__():
                        """initialize Reads object"""
                        return {"data": dict(), "int_n_removed_elements": 0}

                    def __Initialize_gene_and_isoform_data__(
                        reads, refname_gene, refstart_gene, refend_gene, id_gene
                    ):
                        """# 2022-05-21 16:35:02
                        initialized data for gene and isoform annotation

                        refstart_gene, refend_gene : 0-based coordinates
                        """
                        dict_data = dict()
                        dict_data["annotation_type"] = "gene_and_isoform"
                        dict_data["id_anno"] = id_gene
                        dict_data["refname_anno"] = refname_gene
                        dict_data["refstart_anno"] = refstart_gene
                        dict_data["refend_anno"] = refend_gene
                        dict_data["wall_time"] = 0  # initialize wall_time
                        # initialize an array for counting reads for each cell-barcode
                        dict_data["l_read"] = []

                        """ when using isoform assignment from minimap2 alignment """
                        if not flag_use_isoform_assignment_from_10x_cellranger:
                            # retrieve list of id_tx for the current id_gene
                            dict_data["l_id_tx"] = (
                                scidx["dict_id_gene_to_l_id_tx"][id_gene]
                                if id_gene in scidx["dict_id_gene_to_l_id_tx"]
                                else []
                            )  # handle when no id_tx is available for the id_gene
                            # retrieve fasta files of transcripts, and save as a fasta file
                            # since creating a dictionary using the one-liner expression (for some reason) cause occasional deadlocks (why??) the multi-line version should be used
                            dict_fasta_for_current_gene = __data_object_subset(
                                data_dict_fa_transcriptome, dict_data["l_id_tx"]
                            )
                            # append poly A sequence to the transcript
                            if len( seq_polya_sequence_to_append_to_tx ) > 0 :
                                dict_fasta_for_current_gene = dict( ( k, dict_fasta_for_current_gene[ k ] + seq_polya_sequence_to_append_to_tx ) for k in dict_fasta_for_current_gene ) # append poly A sequence to the transcript (update 'dict_fasta_for_current_gene')
                            dict_data["dict_fasta_tx"] = dict_fasta_for_current_gene
                                
                            dict_data[
                                "path_file_fasta_tx"
                            ] = f"{path_folder_temp}{bk.UUID( )}.tx.fa.gz"
                            bk.FASTA_Write(
                                dict_data["path_file_fasta_tx"],
                                dict_fasta=dict_data["dict_fasta_tx"],
                            )  # save as a fasta file
                            dict_data["am_tx"] = mappy.Aligner(
                                fn_idx_in=dict_data["path_file_fasta_tx"], preset=str_mappy_aligner_preset_for_realignment # initialize the aligner using the given preset
                            )  # load minimap aligner using the minimap2 index with the single-end short read alignment mode

                        # load the initialized data
                        reads["data"][id_gene] = dict_data

                    def __Initialize_misc_anno_data__(
                        reads,
                        refname_anno,
                        refstart_anno,
                        refend_anno,
                        id_anno,
                        annotation_type="miscellaneous",
                        str_mode_scarab_count="gex",
                    ):
                        """# 2022-04-16 21:17:19
                        initialized data for miscellaneous annotation

                        refstart_anno, refend_anno : 0-based coordinates

                        'str_mode_scarab_count' : 'atac' or 'gex'
                        """
                        dict_data = dict()
                        dict_data["id_anno"] = id_anno
                        dict_data["mode_scarab_count"] = str_mode_scarab_count
                        dict_data["annotation_type"] = annotation_type
                        dict_data["refname_anno"] = refname_anno
                        dict_data["refstart_anno"] = refstart_anno
                        dict_data["refend_anno"] = refend_anno
                        dict_data["wall_time"] = 0  # initialize wall_time
                        # initialize an array for counting reads for each cell-barcode
                        dict_data["l_read"] = []

                        # load the initialized data
                        reads["data"][id_anno] = dict_data
                        
                    """
                    define functions for offloading works for multiprocessing
                    """
                    def _initialize_dict_output( ) :
                        ''' # 2023-08-30 11:49:23 
                        initialize dict_output for processing buckets (a container of output from a child process that will be collected by this process)
                        '''
                        dict_output = {
                            'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' : dict( ( e, BytesIO( ) ) for e in l_t_distribution_range_of_interest ),
                            'bio_newfile_df_analysis_statistics' : BytesIO( ),
                        } # initialize the dictionary containing the outputs as binary stream objects                        
                        return dict_output
                                        
                    def _apply_size_distribution_correction_and_export_count_matrix( 
                        df, 
                        t_distribution_range_of_interest : Union[ None, tuple ] = None, 
                        l_col_for_dropping_duplicates : Union[ List[ str ], None ] = l_col_for_identifying_unique_molecules, 
                        l_col_for_groupby_operation : Union[ List[ str ], None ] = [ 'CB' ], 
                        dict_rename_columns : Union[ dict, None ] = {"CB": "barcode"}, 
                        func_set_id_feature = None, 
                        func_set_feature = None, 
                        inplace : bool = False,
                        str_suffix : str = '',
                    ) : 
                        ''' # 2023-08-30 20:09:49 
                        apply size distribution correction in a dataframe with [ "read_count", "int_total_aligned_length" ] columns, and compose a count matrix
                        
                        str_suffix : str = '', # suffix to add to 'name_feature' and 'id_feature'
                        '''
                        flag_size_distribution_correction_is_applied = t_distribution_range_of_interest is not None
                    
                        if not inplace : # copy dataframe
                            df = deepcopy( df )
                            
                        ''' filter molecules using 't_distribution_range_of_interest' '''
                        if flag_size_distribution_correction_is_applied :
                            arr_len = df[ 'int_total_aligned_length' ].values # retrieve length values
                            df = df[ ( arr_len >= t_distribution_range_of_interest[ 0 ] ) & ( arr_len <= t_distribution_range_of_interest[ 1 ] ) ] # filter out molecules outside the given distribution range
                            df.sort_values( "int_total_aligned_length", inplace = True, ascending = False ) # put longer read to the front (for 'drop duplicated molecule' step)
                        
                        ''' drop duplicated molecule '''
                        df.drop_duplicates( subset = l_col_for_dropping_duplicates, keep="first", inplace=True, ) # remove duplicated molecules, keeping the molecules in the front (longer reads)
                        
                        ''' apply size-distribution correction '''
                        if flag_size_distribution_correction_is_applied :
                            arr_read_count = df[ "read_count" ].values.astype( float ) # retrieve "read_count" values # convert to the float datatype
                            for i, int_len in enumerate( df[ "int_total_aligned_length" ].values ) : # apply size-distribution correction
                                arr_read_count[ i ] *= arr_ratio_to_ref[ int_len ]
                            df[ "read_count" ] = arr_read_count
                        
                        ''' dropping columns that are not required for groupby operations '''
                        set_col_to_retain = set( l_col_for_groupby_operation + [ 'read_count' ] ) # set of columns to retain
                        df.drop( columns = list( col for col in df.columns.values if col not in set_col_to_retain ), inplace = True )
                        
                        ''' perform groupby operation   '''
                        df = getattr( df.groupby( l_col_for_groupby_operation ), 'sum' if flag_size_distribution_correction_is_applied else 'count' )( ) # calculate sum of weights for 'normalized counts', and calculate counts for 'raw counts'
                        df.reset_index(drop=False, inplace=True) # reset index after groupby operation
                            
                        ''' compose count matrix '''
                        if dict_rename_columns is not None : # rename columns
                            df.rename( columns = dict_rename_columns, inplace=True )
                        df[ 'id_feature' ] = func_set_id_feature if isinstance( func_set_id_feature, str ) else func_set_id_feature( df ) # set 'id_feature'
                        df[ 'feature' ] = func_set_feature if isinstance( func_set_feature, str ) else func_set_feature( df ) # set 'feature'
                        
                        ''' add suffix '''
                        if isinstance( str_suffix, str ) and len( str_suffix ) > 0 : # if valid 'str_suffix' has been given, add the suffix
                            df[ 'id_feature' ] = df[ 'id_feature' ] + str_suffix
                            df[ 'feature' ] = df[ 'feature' ] + str_suffix
                        df = df[ l_col_df_count ]  # reorder columns
                        return df # return the resulting dataframe containing count data
                    
                    l_col_for_composing_df_count = [ ] if flag_is_mode_scarab_count_atac else [ 'flag_internal_polya_tract_primed' ] # define a additional list of columns for composing 'df_count'
                    def _process_gene_and_isoform_data( id_anno : str, dict_data : dict ):
                        """ # 2023-08-28 16:25:33 
                        Flush data for gene and isoform
                        assumes uniquely aligned reads
                        requires the following columns:
                        l_col_for_counting = [ 'qname', 'mapq', 'flag', 'str_l_seg', 'CB', 'UB', 'TX', 'RE', 'int_flag_classification', 'id_tx_assigned_by_minimap2', 'l_name_variant', 'int_total_aligned_length' ] + l_col_for_composing_df_count
                        
                        id_anno # name of the annotation
                        dict_data # dictionary data of the annotation
                        """
                        float_time_start = time.time()  # record the start time
                        ''' initialize the output '''
                        dict_output = _initialize_dict_output( )

                        """ retrieve 'name_gene' """
                        (
                            seqname_anno,
                            source_anno,
                            feature_anno,
                            start_anno,
                            end_anno,
                            score_anno,
                            strand_anno,
                            frame_anno,
                            attribute_anno,
                            name_anno,
                        ) = scidx["arr_data_df_gtf_gene"][
                            scidx["dict_index_df_gtf_gene"][id_anno][0]
                        ]  # retrieve name_gene # 1-based
                        if isinstance(
                            name_anno, float
                        ):  # if 'name_gene' is not available, use id_gene as 'name_gene'
                            name_anno = id_anno

                        """ when using isoform assignment from minimap2 alignment """
                        # remove temporary files
                        if not flag_use_isoform_assignment_from_10x_cellranger:
                            os.remove(dict_data["path_file_fasta_tx"])

                        int_num_record = len(
                            dict_data["l_read"]
                        )  # retrieve the number of records for the current annotation
                        df_read = pd.DataFrame(
                            dict_data.pop("l_read"), columns=l_col_for_counting
                        )  # compose a dataframe for the current annotation

                        """ filter aligned reads """
                        df_read = bk.PD_Binary_Flag_Select(
                            df_read, "flag", 10, flag_select=False
                        )  # remove PCR or optical duplicate
                        df_read = bk.PD_Binary_Flag_Select(
                            df_read, "flag", 8, flag_select=False
                        )  # remove secondary alignments

                        # remove reads aligned to introns
                        if not flag_include_read_aligned_to_intron and len(df_read) > 0:
                            df_read = (
                                bk.PD_Select(df_read, RE="N", deselect=True)
                                if flag_use_intronic_read_assignment_from_10x_cellranger
                                else bk.PD_Binary_Flag_Select(
                                    df_read,
                                    "int_flag_classification",
                                    2,
                                    flag_select=False,
                                )
                            )  # remove intronic read based on read alignment classification from 10X cellrange or default classification (that of this program)

                        if len(df_read) > 0:
                            """assign id_tx to reads uniquely aligned to transcripts"""
                            if flag_use_isoform_assignment_from_10x_cellranger:
                                # assign id_tx to reads uniquely aligned to a single transcript (use isoform assignment of 10x cellranger)
                                l_tx_assigned = []
                                for e in df_read.TX.values:
                                    if isinstance(e, float) or (
                                        isinstance(e, str) and len(e) == 0
                                    ):
                                        l_tx_assigned.append(np.nan)
                                    elif e.count(";") > 0:
                                        l_tx_assigned.append(np.nan)
                                    else:
                                        l_tx_assigned.append(e.split(",", 1)[0])
                                df_read["id_tx_assigned"] = l_tx_assigned
                            else:
                                df_read[
                                    "id_tx_assigned"
                                ] = df_read.id_tx_assigned_by_minimap2.replace(
                                    "", np.nan
                                )  # use np.nan to indicate a read has not been aligned uniquely to a single transcript

                            """ count at gene-level """
                            df = df_read[
                                ["qname", "mapq", "flag", "CB", "UB", "id_tx_assigned", "int_total_aligned_length"] + l_col_for_composing_df_count
                            ]
                            # remove records without cell barcodes
                            df.loc[:, l_col_for_identifying_unique_molecules] = df[l_col_for_identifying_unique_molecules].replace(
                                "", np.nan
                            )
                            df.dropna(subset=l_col_for_identifying_unique_molecules, inplace=True)

                            # drop duplicates based on read name, and prioritize reads uniquely aligned to transcript sequences
                            df.sort_values(["id_tx_assigned"], inplace=True)
                            df.sort_values(
                                ["mapq"], inplace=True, ascending=False
                            )  # sort reads in the order of an increasing mapq scores (so that alignments with lower mapping scores are removed)
                            df.drop_duplicates(
                                subset=["qname"], keep="first", inplace=True
                            )
                            df.drop(columns=["qname", "mapq"], inplace=True)

                            # initialize counting read for each cell
                            df["read_count"] = 1

                            def __Count_and_Write_gene_data__(
                                df, feature, id_feature
                            ):
                                """count and write gene annotation count data to the given file handle 'newfile_df_count'"""
                                if len(df) == 0:  # detect empty dataframe
                                    return -1
                                df_count = df[ l_col_for_identifying_unique_molecules + [ "read_count", "int_total_aligned_length" ] + l_col_for_composing_df_count ] # subset the data
                                ''' create normalized count matrix '''
                                for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest', compose output
                                    for _suffix, _df_count in zip( [ '', '|excluding_internal_polyA_primed' ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ '' ], [ df_count, df_count[ df_count.flag_internal_polya_tract_primed == False ] ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ df_count ] ) : # exclude internal poly A primed reads during counting
                                        if len( _df_count ) == 0 : # ignore empty dataframe
                                            continue
                                        _apply_size_distribution_correction_and_export_count_matrix( 
                                            _df_count, 
                                            t_distribution_range_of_interest = t_distribution_range_of_interest, 
                                            l_col_for_dropping_duplicates = l_col_for_identifying_unique_molecules, 
                                            l_col_for_groupby_operation = [ 'CB' ], 
                                            dict_rename_columns = {"CB": "barcode"}, 
                                            func_set_id_feature = id_feature, 
                                            func_set_feature = feature,
                                            str_suffix = _suffix,
                                        ).to_csv( dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ][ t_distribution_range_of_interest ], sep="\t", header=None, index=False )

                            """ count strand specific reads """
                            flag_valid_strand_info = strand_anno in [
                                "-",
                                "+",
                            ]  # retrieve a flag indicating whether the gene's strand information is valid
                            if flag_valid_strand_info:
                                flag_reverse_complemented_is_sense = (
                                    strand_anno == "-"
                                )  # retrieve a flag for selecting sense reads
                                if (
                                    flag_include_strand_specific_counts
                                ):  # count reads aligned to each strand separately if 'flag_include_strand_specific_counts' is True # only counts strand specific counts if the strand information of the gene annotation is valid
                                    df_sense = bk.PD_Binary_Flag_Select(
                                        df,
                                        "flag",
                                        4,
                                        flag_select=flag_reverse_complemented_is_sense,
                                    )
                                    __Count_and_Write_gene_data__(
                                        df_sense,
                                        f"{name_anno}|strand=sense",
                                        f"{id_anno}|strand=sense",
                                    )  # export count data for the reads aligned to sense direction
                                    __Count_and_Write_gene_data__(
                                        bk.PD_Binary_Flag_Select(
                                            df,
                                            "flag",
                                            4,
                                            flag_select=not flag_reverse_complemented_is_sense,
                                        ),
                                        f"{name_anno}|strand=antisense",
                                        f"{id_anno}|strand=antisense",
                                    )  # export count data for the reads aligned to anti-sense direction

                                elif not flag_include_read_aligned_to_opposite_strand:
                                    df_sense = bk.PD_Binary_Flag_Select(
                                        df,
                                        "flag",
                                        4,
                                        flag_select=flag_reverse_complemented_is_sense,
                                    )  # include only the sense reads
                            else:  # if valid strand information is not available, consider all reads as aligned to 'sense' direction
                                df_sense = df
                            if (
                                not flag_include_read_aligned_to_opposite_strand
                            ):  # if only sense reads are being counted, use 'df_sense' as 'df'
                                df = df_sense

                            __Count_and_Write_gene_data__(
                                df, name_anno, id_anno
                            )  # count reads # exclude reads aligned to antisense direction if 'flag_include_read_aligned_to_opposite_strand' is False
                            
                            """ write counts of reads not uniquely assigned to a single transcript """
                            __Count_and_Write_gene_data__(
                                df[ pd.isnull( df.id_tx_assigned ) ], f"{name_anno}|not_uniquely_assigned_to_tx", f"{id_anno}|not_uniquely_assigned_to_tx"
                            )  # count reads # exclude reads aligned to antisense direction if 'flag_include_read_aligned_to_opposite_strand' is False

                            """ count at transcript-level (currently counts reads aligned to both sense and antisense strands) """
                            # count read for each cell
                            df.dropna( subset=["id_tx_assigned"], inplace=True )  # retrieve reads uniquely aligned to a single transcript
                            
                            df_tx_count = df[
                                ["CB", "UB", "id_tx_assigned", "read_count", "int_total_aligned_length"] + l_col_for_composing_df_count
                            ]
                            if len(df_tx_count) > 0:  # if a valid isoform count exists
                                ''' functions for mapping identifiers '''
                                def _tx__func_set_id_feature( df ) :
                                    return df["id_tx_assigned"].values
                                mapping = MAP.Map( scidx["dict_id_tx_to_name_tx"] ).a2b_if_mapping_available_else_Map_a2a # retrieve mapping
                                def _tx__func_set_feature( df ) :
                                    return name_anno + "|tx_name=" + df["id_tx_assigned"].apply( mapping ) + "|tx_id=" + df["id_tx_assigned"]
                                ''' create normalized count matrix '''
                                for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest', compose output
                                    for _suffix, _df_count in zip( [ '', '|excluding_internal_polyA_primed' ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ '' ], [ df_tx_count, df_tx_count[ df_tx_count.flag_internal_polya_tract_primed == False ] ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ df_tx_count ] ) : # exclude internal poly A primed reads during counting
                                        if len( _df_count ) == 0 : # ignore empty dataframe
                                            continue
                                        _apply_size_distribution_correction_and_export_count_matrix( 
                                            _df_count, 
                                            t_distribution_range_of_interest = t_distribution_range_of_interest, 
                                            l_col_for_dropping_duplicates = l_col_for_identifying_unique_molecules + [ 'id_tx_assigned' ], 
                                            l_col_for_groupby_operation = [ "CB", "id_tx_assigned" ], 
                                            dict_rename_columns = {"CB": "barcode"}, 
                                            func_set_id_feature = _tx__func_set_id_feature, 
                                            func_set_feature = _tx__func_set_feature,
                                            str_suffix = _suffix,
                                        ).to_csv( dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ][ t_distribution_range_of_interest ], sep="\t", header=None, index=False )

                            """ count gene/transcript counts for each detected variant """
                            if (
                                not flag_does_not_collect_variant_information
                                and flag_output_variant_information_with_annotations
                            ):
                                """compose a dataframe for counting molecules with variants"""
                                df = df_read[
                                    [
                                        "qname",
                                        "str_l_seg",
                                        "CB",
                                        "UB",
                                        "id_tx_assigned",
                                        "l_name_variant",
                                        "int_total_aligned_length",
                                    ] + l_col_for_composing_df_count
                                ]
                                # remove records without cell barcodes
                                df.loc[:, l_col_for_identifying_unique_molecules] = df[l_col_for_identifying_unique_molecules].replace(
                                    "", np.nan
                                )
                                df.dropna(subset=l_col_for_identifying_unique_molecules, inplace=True)
                                df.sort_values(
                                    ["l_name_variant", "id_tx_assigned"],
                                    ascending=False,
                                    inplace=True,
                                )  # prioritize reads with l_name_variant and id_tx_assigned
                                df.drop_duplicates(
                                    subset=["qname"], keep="first", inplace=True
                                )  # retrieve best record for each unique read
                                df.drop(columns=["qname"], inplace=True)

                                """ count and filter mutations """
                                str_id_var_concatenated = ";".join(
                                    df.l_name_variant[
                                        df.l_name_variant.apply(len) > 0
                                    ].values
                                )
                                # retrieve a list of all 'id_var'
                                l_id_var_all = str_id_var_concatenated.split(";")
                                int_min_num_of_reads = (
                                    1
                                    if flag_filter_variant_using_predefined_set
                                    else max(
                                        int(
                                            np.ceil(
                                                len(l_id_var_all)
                                                * float_min_prop_of_reads_for_filtering_genomic_variant
                                            )
                                        ),
                                        int_min_num_of_reads_for_filtering_genomic_variant,
                                    )
                                )  # retrieve a threshold (maximum of the two thresholds). if 'flag_filter_variant_using_predefined_set' is True, use all available (already filtered) variants
                                # retrieve list of id_mut (in the order of decreasing counts)
                                l_id_var = (
                                    bk.LIST_COUNT(
                                        l_id_var_all,
                                        duplicate_filter=int_min_num_of_reads,
                                    ).index.values
                                    if len(str_id_var_concatenated) > 0
                                    else []
                                )
                                del l_id_var_all
                                if (
                                    len(l_id_var) > 0
                                ):  # if a variant is detected for the current feature
                                    """build an interval tree of detected variants"""
                                    it_var = (
                                        intervaltree.IntervalTree()
                                    )  # interval tree for searching variants. values contains mutation, reference
                                    for id_var in l_id_var:
                                        refname, str_refpos, anno_mut = id_var.split(
                                            ":", 2
                                        )
                                        refpos = (
                                            int(str_refpos) - 1
                                        )  # 1-based > 0-based coord
                                        ref, alt = anno_mut.split(">")
                                        len_ref, len_alt = len(ref), len(alt)
                                        it_var[refpos : refpos + len_ref] = (
                                            id_var,
                                            refname + ":" + str_refpos,
                                        )  # update interval tree for searching variant using intervals
                                    """ retrieve a list of unique molecules for each variant information """
                                    dict_unique_molecule_to_max_molecule_size = (
                                        dict()
                                    )
                                    for (
                                        str_l_seg,
                                        CB,
                                        UB,
                                        id_tx_assigned,
                                        l_name_variant,
                                        int_total_aligned_length,
                                    ) in df[
                                        [
                                            "str_l_seg",
                                            "CB",
                                            "UB",
                                            "id_tx_assigned",
                                            "l_name_variant",
                                            "int_total_aligned_length"
                                        ]
                                    ].values:
                                        set_name_variant = (
                                            set(l_name_variant.split(";"))
                                            if len(l_name_variant) > 0
                                            else set()
                                        )  # retrieve set of variants detected in the current reads
                                        """ collect intervals of variants overlapping with the current aligment """
                                        set_overlapped_interval_var = set()
                                        for e in str_l_seg.split(","):
                                            refstart, refend = list(
                                                map(int, e.split("-"))
                                            )
                                            set_overlapped_interval_var |= it_var[
                                                refstart - 1 : refend
                                            ]
                                        for (
                                            start,
                                            end,
                                            values,
                                        ) in (
                                            set_overlapped_interval_var
                                        ):  # for each interval of a variant
                                            id_var, id_allele_ref = values
                                            flag_read_contains_id_allele_ref = (
                                                id_var not in set_name_variant
                                            )
                                            id_allele = (
                                                id_allele_ref
                                                if flag_read_contains_id_allele_ref
                                                else id_var
                                            )  # retrieve allele of the current read based on the detected list of 'name_variants'

                                            """ count at gene-level """
                                            t = (
                                                CB,
                                                UB,
                                                id_allele,
                                                id_anno,
                                                flag_read_contains_id_allele_ref,
                                                id_allele_ref,
                                            )
                                            if t not in dict_unique_molecule_to_max_molecule_size :
                                                dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length
                                            elif int_total_aligned_length > dict_unique_molecule_to_max_molecule_size[ t ] : 
                                                dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length

                                            """ count at transcript-level """
                                            if isinstance( id_tx_assigned, str ):  # if valid id_tx_assigned exists
                                                t = (
                                                    CB,
                                                    UB,
                                                    id_allele,
                                                    id_tx_assigned,
                                                    flag_read_contains_id_allele_ref,
                                                    id_allele_ref,
                                                )
                                                if t not in dict_unique_molecule_to_max_molecule_size :
                                                    dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length
                                                elif int_total_aligned_length > dict_unique_molecule_to_max_molecule_size[ t ] : 
                                                    dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length

                                    if ( len( dict_unique_molecule_to_max_molecule_size ) > 0 ):
                                        """compose a count matrix of molecules containing genomic variants"""
                                        ''' prepare df_var_count '''
                                        df_var_count = pd.Series(
                                            dict_unique_molecule_to_max_molecule_size
                                        ).reset_index(drop=False)
                                        df_var_count.columns = [
                                            "barcode",
                                            "str_umi",
                                            "id_var",
                                            "id_feature",
                                            "flag_read_contains_id_allele_ref",
                                            "id_allele_ref",
                                            "int_total_aligned_length",
                                        ]
                                        df_var_count[ 'read_count' ] = 1 # initialize 'read_count' column
                                        df_var_count.sort_values( "flag_read_contains_id_allele_ref", inplace=True, )  # put read that does not contains id_allele_ref at the front of the dataframe # for some cases, molecules with same barcode and umi has different alleles, due to PCR chimera formation and sequencing errors. assuming variants are more rare, variants are prioritized over reference alleles to increase the sensitivity. 
                                        
                                        ''' function for annotating the count matrix '''
                                        """ retrieve id_feature > feature mapping (for both gene and transcripts) """
                                        l_id_tx = list( id_tx for id_tx in df.id_tx_assigned.dropna().unique() )  # retrieve list of 'id_tx' in the dataframe for counting variants
                                        dict_id_feature_to_feature = dict( ( id_tx, name_anno + "|tx_name=" + ( scidx["dict_id_tx_to_name_tx"][ id_tx ] if id_tx in scidx["dict_id_tx_to_name_tx"] else id_tx ) + "|tx_id=" + id_tx, ) for id_tx in l_id_tx )  # id_feature > feature mapping for transcripts
                                        dict_id_feature_to_feature[ id_anno ] = name_anno  # id_feature > feature mapping for the current gene
                                        mapping = MAP.Map(dict_id_feature_to_feature).a2b
                                        def _g_tx_var__func_set_id_feature( df ) :
                                            """ compose feature ('id_feature') """
                                            return df.id_feature.apply( mapping ) + "|var_name=" + df.id_var # add 'id_var' as 'var_name' to the feature (name of feature)
                                        def _g_tx_var__func_set_feature( df ) :
                                            return df[ 'id_feature' ]

                                        ''' create normalized count matrix '''
                                        for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest', compose output
                                            _apply_size_distribution_correction_and_export_count_matrix( 
                                                df_var_count, 
                                                t_distribution_range_of_interest = t_distribution_range_of_interest, 
                                                l_col_for_dropping_duplicates = [ "barcode", "str_umi", "id_feature", "id_allele_ref", ], # dropping duplicates for each reference allele for each transcript or the current gene
                                                l_col_for_groupby_operation = ["barcode", "id_var", "id_feature"], # groupby operations for each transcript or the current gene
                                                dict_rename_columns = None, 
                                                func_set_id_feature = _g_tx_var__func_set_id_feature, 
                                                func_set_feature = _g_tx_var__func_set_feature
                                            ).to_csv( dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ][ t_distribution_range_of_interest ], sep="\t", header=None, index=False )
                                        del df_var_count
                                    else:
                                        """debug"""
                                        logger.info( f"{id_anno =}, {name_anno =}, {l_id_var =}" )
                                        df.to_csv(
                                            f"{path_folder_temp}for_debug.gene.{bk.UUID( )}.tsv.gz",
                                            sep="\t",
                                            index=False,
                                        )
                                    del (
                                        it_var,
                                        dict_unique_molecule_to_max_molecule_size,
                                    )
                                del str_id_var_concatenated, l_id_var
                            # release memory
                            del df_read, df, df_tx_count
                            if (
                                not flag_include_read_aligned_to_opposite_strand
                                or flag_include_strand_specific_counts
                            ):
                                del df_sense
                        # update the wall time
                        dict_data["wall_time"] += time.time() - float_time_start
                        # write analysis statistics
                        dict_output[ 'bio_newfile_df_analysis_statistics' ].write(
                            (
                                "\t".join(
                                    list(
                                        map(
                                            str,
                                            [
                                                str_uuid,
                                                id_anno,
                                                dict_data["wall_time"],
                                                int_num_record,
                                            ],
                                        )
                                    )
                                )
                                + "\n"
                            ).encode()
                        )  # 0-based -> 1-based
                        # release memory # remove data of the reads aligned to the annotation (from the memory)
                        del dict_data
                        return dict_output
                    
                    def _process_misc_anno_data( id_anno : str, dict_data : dict ):
                        """# 2023-01-07 23:12:08
                        flush data for miscellaneous annotation
                        assumes uniquely aligned reads
                        required the following columns defined in the 'l_col_for_counting' (which changes depending on the scarab count modes)
                        
                        id_anno # name of the annotation
                        dict_data # dictionary data of the annotation
                        """
                        float_time_start = time.time()  # record the start time
                        ''' initialize the output '''
                        dict_output = _initialize_dict_output( )

                        int_num_record = len(
                            dict_data["l_read"]
                        )  # retrieve the number of records for the current annotation
                        df_read = pd.DataFrame(
                            dict_data.pop("l_read"), columns=l_col_for_counting
                        )  # compose a dataframe for the current annotation

                        """ filter aligned reads """
                        df_read = bk.PD_Binary_Flag_Select(
                            df_read, "flag", 10, flag_select=False
                        )  # remove PCR or optical duplicate
                        df_read = bk.PD_Binary_Flag_Select(
                            df_read, "flag", 8, flag_select=False
                        )  # remove secondary alignments

                        """ retrieve annotation-type specific settings """
                        flag_variant = dict_data["annotation_type"] == "variant"

                        if len(df_read) > 0:
                            """
                            Export counts of reads aligned to the annotation
                            """
                            if (
                                not flag_variant
                            ):  # for 'variant' annotation type, exporting counts of all reads with all variants will be skipped
                                """retain unique reads"""
                                df = df_read[
                                    ["qname", "flag", 'int_total_aligned_length']
                                    + l_col_for_identifying_unique_molecules + l_col_for_composing_df_count
                                ]
                                # remove records without cell barcodes
                                df.loc[:, l_col_for_identifying_unique_molecules] = df[
                                    l_col_for_identifying_unique_molecules
                                ].replace("", np.nan)
                                df.dropna(
                                    subset=l_col_for_identifying_unique_molecules,
                                    inplace=True,
                                )

                                # drop duplicates simply based on read name
                                df.drop_duplicates(
                                    subset=["qname"], keep="first", inplace=True
                                )
                                df.drop(columns=["qname"], inplace=True)

                                # initialize counting read for each cell
                                df["read_count"] = 1

                                """ if the annotation_type is 'gene' retrieve 'name_gene' from 'id_gene', and use 'name_gene' for flushing data """
                                name_anno = np.nan
                                if dict_data["annotation_type"] == "gene":
                                    (
                                        seqname_anno,
                                        source_anno,
                                        feature_anno,
                                        start_anno,
                                        end_anno,
                                        score_anno,
                                        strand_anno,
                                        frame_anno,
                                        attribute_anno,
                                        name_anno,
                                    ) = scidx["arr_data_df_gtf_gene"][
                                        scidx["dict_index_df_gtf_gene"][id_anno][0]
                                    ]  # retrieve name_gene # 1-based
                                if isinstance(
                                    name_anno, float
                                ):  # if 'name_anno' is not available, use id_anno as 'name_anno'
                                    name_anno = id_anno

                                def __Count_and_Write_misc_anno_data__( df, id_anno, name_anno, ):
                                    """count and write misc anno data to the given file handle 'newfile_df_count'"""
                                    if len(df) == 0:  # detect empty dataframe
                                        return -1
                                    df_count = df[ l_col_for_identifying_unique_molecules + [ "read_count", "int_total_aligned_length" ] + l_col_for_composing_df_count ] # subset the data
                                    ''' create normalized count matrix '''
                                    for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest', compose output
                                        for _suffix, _df_count in zip( [ '', '|excluding_internal_polyA_primed' ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ '' ], [ df_count, df_count[ df_count.flag_internal_polya_tract_primed == False ] ] if flag_export_internal_polyA_filtered_count_as_a_separate_feature else [ df_count ] ) : # exclude internal poly A primed reads during counting
                                            if len( _df_count ) == 0 : # ignore empty dataframe
                                                continue
                                            _apply_size_distribution_correction_and_export_count_matrix( 
                                                _df_count, 
                                                t_distribution_range_of_interest = t_distribution_range_of_interest, 
                                                l_col_for_dropping_duplicates = l_col_for_identifying_unique_molecules, 
                                                l_col_for_groupby_operation = [ 'CB' ], 
                                                dict_rename_columns = {"CB": "barcode"}, 
                                                func_set_id_feature = id_anno, 
                                                func_set_feature = name_anno,
                                                str_suffix = _suffix,
                                            ).to_csv( dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ][ t_distribution_range_of_interest ], sep="\t", header=None, index=False )

                                __Count_and_Write_misc_anno_data__(
                                    df,
                                    id_anno,
                                    name_anno,
                                )  # count reads
                                """ count strand specific reads """
                                if (
                                    flag_include_strand_specific_counts
                                ):  # count reads aligned to each strand separately if 'flag_include_strand_specific_counts' is True # for variant annotation type, strand-specific count will not be included
                                    __Count_and_Write_misc_anno_data__(
                                        bk.PD_Binary_Flag_Select(
                                            df, "flag", 4, flag_select=False
                                        ),
                                        f"{id_anno}|strand=+",
                                        f"{id_anno}|strand=+",
                                    )  # select non-reverse complemented reads for + strand reads
                                    __Count_and_Write_misc_anno_data__(
                                        bk.PD_Binary_Flag_Select(
                                            df, "flag", 4, flag_select=True
                                        ),
                                        f"{id_anno}|strand=-",
                                        f"{id_anno}|strand=-",
                                    )  # select reverse complemented reads for - strand reads

                            """ count counts for each detected variant for the current miscellaneous annotation """
                            if not flag_does_not_collect_variant_information and (
                                flag_output_variant_information_with_annotations
                                or flag_variant
                            ):  # export variant information for 'variant' annotation type
                                """compose a dataframe for counting molecules with variants"""
                                df = df_read[
                                    ["qname", "str_l_seg", "l_name_variant", "int_total_aligned_length"]
                                    + l_col_for_identifying_unique_molecules + l_col_for_composing_df_count
                                ]
                                # remove records without cell barcodes
                                df.loc[:, l_col_for_identifying_unique_molecules] = df[
                                    l_col_for_identifying_unique_molecules
                                ].replace("", np.nan)
                                df.dropna(
                                    subset=l_col_for_identifying_unique_molecules,
                                    inplace=True,
                                )
                                df.sort_values(
                                    "l_name_variant", ascending=False, inplace=True
                                )  # prioritize reads with l_name_variant and id_tx_assigned
                                df.drop_duplicates(
                                    subset=["qname"], keep="first", inplace=True
                                )  # retrieve best record for each unique read
                                df.drop(columns=["qname"], inplace=True)
                                """ count and filter mutations """
                                str_id_var_concatenated = ";".join(
                                    df.l_name_variant[
                                        df.l_name_variant.apply(len) > 0
                                    ].values
                                )
                                l_id_var_all = str_id_var_concatenated.split(
                                    ";"
                                )  # retrieve a list of all 'id_var'
                                int_min_num_of_reads = (
                                    1
                                    if flag_filter_variant_using_predefined_set
                                    else max(
                                        int(
                                            np.ceil(
                                                len(l_id_var_all)
                                                * float_min_prop_of_reads_for_filtering_genomic_variant
                                            )
                                        ),
                                        int_min_num_of_reads_for_filtering_genomic_variant,
                                    )
                                )  # retrieve a threshold (maximum of the two thresholds). if 'flag_filter_variant_using_predefined_set' is True, use all available (already filtered) variants
                                # retrieve list of id_mut (in the order of decreasing counts)
                                l_id_var = (
                                    bk.LIST_COUNT(
                                        l_id_var_all,
                                        duplicate_filter=int_min_num_of_reads,
                                    ).index.values
                                    if len(str_id_var_concatenated) > 0
                                    else []
                                )
                                if flag_variant:
                                    pos_name = id_anno.split("variant|pos=", 1)[
                                        1
                                    ]  # retrieve the name of the genomic position of the current 'variant' annotation type feature
                                    l_id_var = list(
                                        e
                                        for e in l_id_var
                                        if pos_name == e.rsplit(":", 1)[0]
                                    )  # exclude variants that are outside of the genomic position of the current 'variant' annotation type feature
                                del l_id_var_all
                                if (
                                    len(l_id_var) > 0
                                ):  # if a variant is detected for the current feature
                                    """build an interval tree of detected variants"""
                                    it_var = (
                                        intervaltree.IntervalTree()
                                    )  # interval tree for searching variants. values contains mutation, reference
                                    for id_var in l_id_var:
                                        refname, str_refpos, anno_mut = id_var.split(
                                            ":", 2
                                        )
                                        refpos = (
                                            int(str_refpos) - 1
                                        )  # 1-based > 0-based coord
                                        ref, alt = anno_mut.split(">")
                                        len_ref, len_alt = len(ref), len(alt)
                                        it_var[refpos : refpos + len_ref] = (
                                            id_var,
                                            refname + ":" + str_refpos,
                                        )  # update interval tree for searching variant using intervals
                                    """ retrieve a list of unique molecules for each variant information """
                                    dict_unique_molecule_to_max_molecule_size = (
                                        dict()
                                    )
                                    for arr in df[
                                        ["str_l_seg", "l_name_variant", "int_total_aligned_length" ]
                                        + l_col_for_identifying_unique_molecules
                                    ].values:
                                        (
                                            str_l_seg,
                                            l_name_variant,
                                            int_total_aligned_length,
                                            CB,
                                            t_unique_identifier_for_a_cell,
                                        ) = (arr[0], arr[1], arr[2], arr[3], tuple(arr[4:]))
                                        set_name_variant = (
                                            set(l_name_variant.split(";"))
                                            if len(l_name_variant) > 0
                                            else set()
                                        )  # retrieve set of variants detected in the current reads
                                        """ collect intervals of variants overlapping with the current aligment """
                                        set_overlapped_interval_var = set()
                                        for e in str_l_seg.split(","):
                                            refstart, refend = list(
                                                map(int, e.split("-"))
                                            )
                                            set_overlapped_interval_var |= it_var[
                                                refstart - 1 : refend
                                            ]
                                        """ iterate each overlapped interval of a variant """
                                        for (
                                            start,
                                            end,
                                            values,
                                        ) in set_overlapped_interval_var:
                                            id_var, id_allele_ref = values
                                            flag_read_contains_id_allele_ref = (
                                                id_var not in set_name_variant
                                            )
                                            id_allele = (
                                                id_allele_ref
                                                if flag_read_contains_id_allele_ref
                                                else id_var
                                            )  # retrieve allele of the current read based on the detected list of 'name_variants'

                                            """ count at annotation-level """
                                            t = (
                                                CB,
                                                t_unique_identifier_for_a_cell,
                                                id_allele,
                                                flag_read_contains_id_allele_ref,
                                                id_allele_ref,
                                            )  # this tuple represent a unique molecule with a specific variant.
                                            if t not in dict_unique_molecule_to_max_molecule_size :
                                                dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length
                                            elif int_total_aligned_length > dict_unique_molecule_to_max_molecule_size[ t ] : 
                                                dict_unique_molecule_to_max_molecule_size[ t ] = int_total_aligned_length
                                                
                                    if len( dict_unique_molecule_to_max_molecule_size ) > 0 :
                                        """compose a count matrix of molecules containing genomic variants"""
                                        ''' prepare df_var_count '''
                                        df_var_count = pd.Series(
                                            dict_unique_molecule_to_max_molecule_size
                                        ).reset_index(drop=False)
                                        df_var_count.columns = [
                                            "barcode",
                                            "t_unique_identifier_for_a_cell",
                                            "id_var",
                                            "flag_read_contains_id_allele_ref",
                                            "id_allele_ref",
                                            "int_total_aligned_length",
                                        ] 
                                        df_var_count[ 'read_count' ] = 1 # initialize 'read_count' column
                                        df_var_count.sort_values( "flag_read_contains_id_allele_ref", inplace=True, )  # put read that does not contains id_allele_ref at the front of the dataframe # for some cases, molecules with same barcode and umi has different alleles, due to PCR chimera formation and sequencing errors. assuming variants are more rare, variants are prioritized over reference alleles to increase the sensitivity. 

                                        ''' function for annotating the count matrix '''
                                        def _misc_var__func_set_id_feature( df ) :
                                            return id_anno + "|var_name=" + df.id_var # add 'id_var' as 'var_name' to the feature (name of feature)
                                        def _misc_var__func_set_feature( df ) :
                                            return df[ 'id_feature' ]

                                        ''' create normalized count matrix '''
                                        for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest', compose output
                                            _apply_size_distribution_correction_and_export_count_matrix( 
                                                df_var_count, 
                                                t_distribution_range_of_interest = t_distribution_range_of_interest, 
                                                l_col_for_dropping_duplicates = [ "barcode", "t_unique_identifier_for_a_cell", "id_allele_ref", ], # dropping duplicates for each reference allele
                                                l_col_for_groupby_operation = [ "barcode", "id_var" ], # groupby for each variant detected
                                                dict_rename_columns = None, 
                                                func_set_id_feature = _misc_var__func_set_id_feature, 
                                                func_set_feature = _misc_var__func_set_feature
                                            ).to_csv( dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ][ t_distribution_range_of_interest ], sep="\t", header=None, index=False )
                                        del df_var_count
                                    else:
                                        """debug"""
                                        logger.info(f"{id_anno =}, {l_id_var = }")
                                        df.to_csv(
                                            f"{path_folder_temp}for_debug.misc.{bk.UUID( )}.tsv.gz",
                                            sep="\t",
                                            index=False,
                                        )
                                    del (
                                        it_var,
                                        dict_unique_molecule_to_max_molecule_size,
                                    )
                                del str_id_var_concatenated, l_id_var
                            # release memory
                            del df_read, df
                        # update the wall time
                        dict_data["wall_time"] += time.time() - float_time_start
                        # write analysis statistics
                        dict_output[ "bio_newfile_df_analysis_statistics" ].write(
                            (
                                "\t".join(
                                    list(
                                        map(
                                            str,
                                            [
                                                str_uuid,
                                                id_anno,
                                                dict_data["wall_time"],
                                                int_num_record,
                                            ],
                                        )
                                    )
                                )
                                + "\n"
                            ).encode()
                        )  # 0-based -> 1-based
                        # release memory # remove data of the reads aligned to the annotation (from the memory)
                        del dict_data
                        return dict_output # return the result

                    def _process_buckets( l_id_anno : List[ str ], l_dict_data : List[ dict ] ) :
                        """ # 2023-08-28 23:34:17 
                        process data in a bucket, and return the output
                        """
                        l_dict_output = list( ) # initialize 'l_dict_output'
                        for id_anno, dict_data in zip( l_id_anno, l_dict_data ) : # iterate each bucket
                            dict_output = ( _process_gene_and_isoform_data if ( dict_data["annotation_type"] == "gene_and_isoform" ) else _process_misc_anno_data )( id_anno, dict_data ) # process the data in the the bucket and retrieve the output
                            l_dict_output.append( dict_output ) # collect the result
                        return l_dict_output # return the collected result
                    
                    l_newfile = list( dict_t_distribution_range_of_interest_to_newfile_df_count[ e ] for e in l_t_distribution_range_of_interest ) + [ newfile_df_analysis_statistics ] # compose a linear list of file handles for output TSV files
                    def _post_process_buckets( l_dict_output : List[ dict ] ) :
                        """ # 2023-08-30 12:22:32 
                        perform post-processing of the bucket analysis result
                        (1) write the result of the processed bucket and (2) update the summary metrics using the analysis result of the bucket.
                        """
                        for dict_output in l_dict_output : # for each output
                            dict_output_count = dict_output[ 'dict_t_distribution_range_of_interest_to_bio_newfile_df_count' ]
                            # export the output
                            for b, f in zip( 
                                list( dict_output_count[ e ] for e in l_t_distribution_range_of_interest ) + [ dict_output[ 'bio_newfile_df_analysis_statistics' ], ],
                                l_newfile,
                            ) : # map the list of output file handle and bytes streams in 'dict_output'
                                b.seek( 0 ) # rewind the tape
                                shutil.copyfileobj( b, f, length = 131072 ) # write the buffer to the file # using buffer size of '2 ** 17 = 131072'

                    def _write_results_from_offloaded_works( flag_wait_all : bool = False ) :
                        """ # 2023-08-28 23:20:22 
                        flag_wait_all : bool = False # if True, wait until all processes completed their works, if False, write results currently contained in the workers object.
                        """
                        for res in ( workers_for_bucket_processing.wait_all if flag_wait_all else workers_for_bucket_processing.collect_results )( flag_return_results = True ).values( ) : # 'flag_wait_all' is True, wait until all processes completed their works. # wait for all submitted works to be completed, and retrieve results for each work
                            _post_process_buckets( res ) # post-processing of the output
                    
                    def _initialize_a_batch_of_buckets( ) :
                        """ # 2023-08-28 23:15:05 
                        """
                        ns[ 'l_dict_data' ] = [ ] # a data container for bulk-processing
                        ns[ 'l_id_anno' ] = [ ] # a data container for bulk-processing
                        ns[ 'int_total_num_records_in_a_batch_of_buckets' ] = 0 # initialize the total number of records in a batch of buckets
                    _initialize_a_batch_of_buckets( ) # initialize the bucket

                    def _flush_the_current_batch_of_buckets( ) :
                        """ # 2023-08-27 18:11:36 
                        flush the current batch of buckets
                        
                        flag_wait_all : bool = False # if True, wait until all processes completed their works, if False, write results currently contained in the workers object.
                        """
                        if not workers_for_bucket_processing.is_worker_available : # if all workers are working, wait for a while until all workers are idle 
                            _write_results_from_offloaded_works( flag_wait_all = False ) # flush results from offloaded computations, without waiting all works to be completed 
                            _write_results_from_offloaded_works( flag_wait_all = True ) # flush results from offloaded computations, waiting all works to be completed 
                        elif workers_for_bucket_processing.int_num_completed_results >= int_max_num_batches_in_the_result_container_before_flushing : # if the result container became too large, empty the container
                            _write_results_from_offloaded_works( flag_wait_all = False ) # flush results from offloaded computations, without waiting all works to be completed 
                        workers_for_bucket_processing.submit_work( _process_buckets, args = ( ns[ 'l_id_anno' ], ns[ 'l_dict_data' ] ) ) # submit the work for offloading (append the list of pysam objects as the data associated with the work) # flush the current batch of the buckets
                        _initialize_a_batch_of_buckets( ) # initialize the next batch of the buckets
                        
                    def _empty_bucket( id_anno : str ):
                        """ # 2023-08-28 23:25:20 
                        empty bucket for the 'id_anno' for processing reads assigned to the 'id_anno'
                        """
                        ''' empty the bucket '''
                        dict_data = reads["data"].pop(id_anno) # retrieve data
                        reads[ "int_n_removed_elements" ] += 1  # update the number of elements deleted
                        
                        ''' append the bucket to the batch of buckets that are currently being collected '''
                        ns[ 'l_id_anno' ].append( id_anno )
                        ns[ 'l_dict_data' ].append( dict_data )
                        ns[ 'int_total_num_records_in_a_batch_of_buckets' ] += len( dict_data["l_read"] ) # update the number of reads in a batch
                        if ns[ 'int_total_num_records_in_a_batch_of_buckets' ] > int_max_num_records_in_a_batch_of_buckets : # if the number of records in the batch of buckets exceed the limit, flush the current batch
                            _flush_the_current_batch_of_buckets( ) # flush the current batch of the buckets
                            
                        """ recreate dictionary to avoid 'memory dictionary' if the number of 'pop' operations exceeds certain threshold """
                        if int_max_n_removed_elements < reads["int_n_removed_elements"]:
                            d = reads["data"]
                            reads["data"] = dict((e, d[e]) for e in d)
                            reads[
                                "int_n_removed_elements"
                            ] = 0  # reset the number of pop operations
                            
                    '''
                    iterate a portion of BAM file
                    '''
                    while True:
                        ins = pipe_receiver.recv()
                        if ins is None:
                            break
                        name_contig = ins  # parse input
                        # initialize the output
                        int_n_sam_record_count = 0
                        if verbose:
                            logger.info( f"processing of '{name_contig}' started." )
                        reads = __Initialize_Reads__()  # initialize an object that contains a dictionary to collect information about the analyzed reads for each gene, repeatmasker, regulatory element, and genomic region (bin) annotations to which the reads were assigned.
                        with pysam.AlignmentFile(
                            path_file_bam_input,
                            str_mode_sam,
                            reference_filename=path_file_fa_for_cram,
                        ) as samfile:
                            for r in samfile.fetch( contig = name_contig ) :
                                """retrieve attributes of the current alignment and filter the alignment based on mapq, flag, and position (ignore reads outside start/end sites, which often occurs due to samtools's algorithm)"""
                                (
                                    refstart,
                                    refend,
                                    qname,
                                    refname,
                                    seq,
                                    flag,
                                    cigartuples,
                                    int_mapq,
                                ) = (
                                    r.reference_start,
                                    r.reference_end,
                                    r.qname,
                                    r.reference_name,
                                    r.seq,
                                    r.flag,
                                    r.cigartuples,
                                    r.mapq,
                                )  # 0-based coordinates # retrieve mapping quality
                                if flag_is_5prime:
                                    if flag_is_paired_end:
                                        raise NotImplementedError(
                                            "5prime-paired-end currently not supported."
                                        )
                                    else:
                                        flag ^= (
                                            1 << 4
                                        )  # flip the read's direction if the technology is 5prime and single-end
                                # check whether the read was reverse complemented
                                flag_is_reverse_complemented = _check_binary_flags( flag, 4 ) 

                                """ filter aligned records using mapping quality (this filter out reads with invalid cigar string, too) """
                                if (
                                    int_mapq < int_min_mapq_unique_mapped
                                ):  # skip read whose mapq is below 'int_min_mapq_unique_mapped'
                                    continue

                                """ filter aligned records using flags """
                                # ignore records with flag indicating secondary alignment
                                if (
                                    flag_ignore_record_with_flag_secondary_alignment
                                    and (flag & (1 << 8))
                                ):
                                    continue
                                # ignore records with flag indicating the alignment is optical pcr duplicates
                                if (
                                    flag_ignore_record_with_flag_optical_or_pcr_duplicate
                                    and (flag & (1 << 10))
                                ):
                                    continue

                                """ assign 'int_pos_of_read_determining_feature_assignment' """
                                int_pos_of_read_determining_feature_assignment = (
                                    (refend - 1 if (flag & 1 << 4) else refstart)
                                    if flag_is_mode_scarab_count_atac
                                    else int((refstart + refend) / 2)
                                )  # the middle position will be used to identify genomic region that represent the current read most accurately

                                """ initialize read analysis """
                                int_n_sam_record_count += 1 # increase the counter
                                refname = __chromosome_name_remove_chr__(
                                    refname
                                )  # remove 'chr' prefix from the reference name
                                # retrieve length of sequence
                                len_seq = (
                                    -1 if seq is None else len(seq)
                                )  # if seq is not included in the record, write -1 as the length of the sequence
                                # retrieve a dictionary of SAM tags
                                dict_tags = dict(r.tags)
                                
                                # retrieve useful flags
                                flag_internal_polya_tract_primed = dict_tags[ str_name_bam_tag_internal_polya_length ] > 0 if str_name_bam_tag_internal_polya_length in dict_tags else False # if 'str_name_bam_tag_internal_polya_length' tag is available, classify the alignment as internal-polyA-tract-primed read if the detected length of the internal polyA tract is longer than 0. if the tag is not available, assume all read is not primed by internal polyA-tract

                                # initialize the dictionary of sam tags
                                for key in [
                                    str_name_bam_tag_cb_corrected,
                                    str_name_bam_tag_cb_uncorrected,
                                    str_name_bam_tag_umi_corrected,
                                    str_name_bam_tag_umi_uncorrected,
                                    str_name_bam_tag_internal_polya_length,
                                    "TX",
                                    "AN",
                                    "GX",
                                    "GN",
                                    "MM",
                                    "RE",
                                    "xf",
                                ]:
                                    if key not in dict_tags:
                                        dict_tags[key] = ""

                                """ "Flush" count data for the annotations that are no longer overlaps with the current genomic position (thus no reads overlapping with the genes will be encountered from this point, since the alingment is sorted by refstart position) """
                                """
                                'Flush' data
                                """
                                for id_anno in list(reads["data"]):
                                    dict_data = reads["data"][
                                        id_anno
                                    ]  # retrieve data
                                    if (
                                        dict_data["refname_anno"] != refname
                                        or dict_data["refend_anno"] < refstart
                                    ):  # flush data for annotation in the previous contigs, annotation that has been 'expired' during the iteration of refstart-sorted BAM file
                                        _empty_bucket( id_anno, )

                                """
                                Overall Structure of Scarab short-read

                                reads are classified, and 'reads' object initialized for each feature
                                at the end of the analysis of the read, the data of the read is appended to all the features of the 'reads' object, which will be converted to as a count matrix once each feature is flushed.

                                """
                                """ initialize """
                                int_flag_classification = (
                                    0  # initialize binary classification flag
                                )
                                (
                                    id_rpmk,
                                    int_max_num_base_pairs_overlap_with_rpmk,
                                    id_gene,
                                    id_promoter,
                                    int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id,
                                    int_base_gene_exon_count,
                                    int_base_filtered_rpmk_count,
                                    int_base_unfiltered_rpmk_count,
                                    int_base_reg_count,
                                    id_tx_assigned_by_minimap2,
                                    id_reg,
                                    id_bin_genome,
                                ) = ("", "", "", "", "", "", "", "", "", "", "", "")
                                res_overlapped_exons = None  # initialize variables
                                l_id_anno_exon_and_splice_junc = (
                                    []
                                )  # a list that will contain valie id_anno_exon and id_anno_splice_junc of the current read if the counting behavior has been enabled.
                                start_gene, end_gene = None, None
                                l_start_end_read = [refstart, refend]

                                """ retrieve mapped segments """
                                (
                                    l_seg,
                                    int_total_aligned_length,
                                ) = SAM.Retrive_List_of_Mapped_Segments(
                                    cigartuples, pos_start=refstart
                                )  # 0-based coordinates
                                str_l_seg = ",".join(
                                    list(
                                        str(t_seg[0] + 1) + "-" + str(t_seg[1])
                                        for t_seg in l_seg
                                    )
                                )  # retrieve string representation of mapped segments # (0-based > 1-based coordinates)

                                """ for ATAC-seq analysis only consider the cut site (more specifically, the base next to the cut site that is included in the read) """
                                if (
                                    flag_is_mode_scarab_count_atac
                                ):  # use the region containing the cutsite for 'l_start_end_read' and 'l_seg'
                                    l_start_end_read = [
                                        int_pos_of_read_determining_feature_assignment,
                                        int_pos_of_read_determining_feature_assignment
                                        + 1,
                                    ]
                                    l_seg = [l_start_end_read]
                                    int_total_aligned_length = 1

                                """
                                check overlap with promoter regions (ATAC mode specific)
                                """
                                if (
                                    flag_is_mode_scarab_count_atac
                                ):  # specific to ATAC-seq data
                                    float_time_start = (
                                        time.time()
                                    )  # record the start time
                                    # ignore reads if it does not overlap repeatmasker annotation
                                    set_t_interval_overlap_promoter = set(
                                        tuple(e)
                                        for e in __data_object_search_query(
                                            data_dict_it_promoter,
                                            refname,
                                            int_pos_of_read_determining_feature_assignment,
                                        )
                                    )  # retrieve set of promoter annotations overlapped with the current alignment
                                    if (
                                        len(set_t_interval_overlap_promoter) > 0
                                    ):  # when current read overlaps with a promoter region
                                        int_flag_classification ^= int_flag_class_overlap_with_promoter  # update read classification info.

                                        # identify promoter region whose TSS is closest to the current read's 'int_pos_of_read_determining_feature_assignment'
                                        (
                                            l_interval_promoter,
                                            int_max_num_base_pairs_overlap_with_promoter,
                                        ) = bk.DICTIONARY_Find_keys_with_min_value(
                                            dict(
                                                (
                                                    e,
                                                    e[1]
                                                    - int_pos_of_read_determining_feature_assignment
                                                    if e[2][1] == "+"
                                                    else int_pos_of_read_determining_feature_assignment
                                                    - e[0]
                                                    + 1,
                                                )
                                                for e in set_t_interval_overlap_promoter
                                            )
                                        )  # (1) retrieve the distance from the TSS to the 'int_pos_of_read_determining_feature_assignment' for each promter region overlapped with the current 'int_pos_of_read_determining_feature_assignment', and (2) select the promoter whose TSS is closest to the 'int_pos_of_read_determining_feature_assignment'
                                        if len(l_interval_promoter) == 1:
                                            (
                                                start_promoter,
                                                end_promoter,
                                                arr,
                                            ) = l_interval_promoter[
                                                0
                                            ]  # retrieve information about the assigned region
                                            (
                                                id_gene_promoter,
                                                strand_promoter,
                                            ) = arr  # retrieve information about the assigned promoter
                                            id_promoter = f"promoter|gene_id={id_gene_promoter}|pos={refname}:{start_promoter + 1}-{end_promoter}"  # 0>1 based coordinates # compose id_promoter

                                            """ initialize data for the new promoter annotation encountered """
                                            if (
                                                id_promoter not in reads["data"]
                                            ):  # if a annotation is newly encountered
                                                __Initialize_misc_anno_data__(
                                                    reads,
                                                    refname,
                                                    start_promoter,
                                                    end_promoter,
                                                    id_promoter,
                                                    "promoter",
                                                    str_mode_scarab_count,
                                                )

                                            """ increament the time passed to process the read to the total wall time for the gene """
                                            reads["data"][id_promoter][
                                                "wall_time"
                                            ] += (time.time() - float_time_start)
                                        else:
                                            int_flag_classification ^= int_flag_class_ambiguous_assignment_to_promoter  # update read classification info.

                                """ assign gene to the current read (GEX and ATAC data shares the same algorithm for gene_id assignment) """
                                if flag_use_gene_assignment_from_10x_cellranger_for_the_current_bam_file:
                                    """use gene assignment from 10x cellranger"""
                                    if len(dict_tags["GX"]) > 0:
                                        int_flag_classification ^= int_flag_class_overlap_with_gene_body  # set 'gene_body' flag
                                        if (
                                            ";" in dict_tags["GX"]
                                            or dict_tags["GX"]
                                            not in scidx["dict_index_df_gtf_gene"]
                                        ):
                                            """when 10x cellranger assigned multiple genes or the 10x cellranger-assigned gene_id does not exist in the given gtf gene annotations"""
                                            int_flag_classification ^= int_flag_class_ambiguous_assignment_to_gene
                                        else:
                                            """for valid id_gene assignment"""
                                            id_gene = dict_tags["GX"]
                                else:
                                    """identify overlapping gene annotations"""
                                    if refname in scidx["dict_it_gene"]:
                                        set_t_interval_overlap_gene = set(
                                            tuple(e)
                                            for e in scidx["dict_it_gene"][refname][
                                                int_pos_of_read_determining_feature_assignment
                                            ]
                                        )  # retrieve set of gene annotations overlapped with the current alignment, based on the middle position of the alignment
                                        if (
                                            len(set_t_interval_overlap_gene) > 0
                                        ):  # flag indicating whether gene is overlapped with the current read
                                            int_flag_classification ^= int_flag_class_overlap_with_gene_body  # set 'gene_body' flag
                                            if (
                                                len(set_t_interval_overlap_gene)
                                                == 1
                                            ):
                                                """assign gene_id if read overlaps with a single gene body"""
                                                id_gene = list(
                                                    set_t_interval_overlap_gene
                                                )[0][2]
                                            else:
                                                """assign gene_id if read overlaps with gene bodies of more then two genes by searching overlapping exons directly"""
                                                """ retrieve a list of 'id_gene' classifed as overlapping with the current read ('the middle position' of the aligned read) to filter exons that does not belong to the classified genes """
                                                set_id_gene_overlapped = set(
                                                    t[2]
                                                    for t in set_t_interval_overlap_gene
                                                )
                                                """ count the number of overlapped base pairs with exons of each gene overlaps with the current read """
                                                dict_id_gene_to_int_num_base_pairs_overlap_with_exons = (
                                                    dict()
                                                )

                                                res_overlapped_exons = __data_object_search_queries(
                                                    data_dict_it_exon,
                                                    refname,
                                                    list(
                                                        slice(start_seg, end_seg)
                                                        for start_seg, end_seg in l_seg
                                                    ),
                                                )  # for each segment of the current read, identify overlapped exons
                                                for i_seg in range(
                                                    len(l_seg)
                                                ):  # for each segment
                                                    start_seg, end_seg = l_seg[
                                                        i_seg
                                                    ]  # retrieve positions of the segment
                                                    res_overlapped_exons_of_a_seg = res_overlapped_exons[
                                                        i_seg
                                                    ]  # retrieve exons overlapped with each segment
                                                    for (
                                                        start_exon,
                                                        end_exon,
                                                        id_gene_of_the_current_exon,
                                                    ) in res_overlapped_exons_of_a_seg:  # for each overlapped exon
                                                        """ignore overlaps with exons of the gene whose gene body does not overlaps with 'int_pos_of_read_determining_feature_assignment' of the aligned read"""
                                                        if (
                                                            id_gene_of_the_current_exon
                                                            not in set_id_gene_overlapped
                                                        ):
                                                            continue
                                                        """ for each exon overlapping with the current segment, add the number of base pairs of overlap to the gene_id to which the current exon belongs to """
                                                        if (
                                                            id_gene_of_the_current_exon
                                                            not in dict_id_gene_to_int_num_base_pairs_overlap_with_exons
                                                        ):
                                                            dict_id_gene_to_int_num_base_pairs_overlap_with_exons[
                                                                id_gene_of_the_current_exon
                                                            ] = 0  # initialize base count for the gene_id
                                                        dict_id_gene_to_int_num_base_pairs_overlap_with_exons[
                                                            id_gene_of_the_current_exon
                                                        ] += bk.INTERVAL_Overlap(
                                                            [start_exon, end_exon],
                                                            [start_seg, end_seg],
                                                            flag_0_based_coordinate_system=True,
                                                        )  # add the number of base pairs of overlap to the gene_id to which the current exon belongs to
                                                """ (GEX mode specific) if the strand to which read was aligned is different from the gene annotation strand, filter out the gene_id from the possible list of gene_ids that can be assigned to the current read. if strand specific sequencing information is not available, does not filter possible list of genes using the information """
                                                if (
                                                    not ( flag_is_mode_scarab_count_atac or flag_include_read_aligned_to_opposite_strand ) # if 'flag_include_read_aligned_to_opposite_strand' is True, ignore the strand information of the read.
                                                ):
                                                    strand_read = (
                                                        "-"
                                                        if flag & (1 << 4)
                                                        else "+"
                                                    )  # retrieve 'strand' from which the current read was aligned
                                                    for (
                                                        id_gene_overlapped_with_read
                                                    ) in list(
                                                        dict_id_gene_to_int_num_base_pairs_overlap_with_exons
                                                    ):
                                                        (
                                                            seqname_anno,
                                                            source_anno,
                                                            feature_anno,
                                                            start_anno,
                                                            end_anno,
                                                            score_anno,
                                                            strand_anno,
                                                            frame_anno,
                                                            attribute_anno,
                                                            name_anno,
                                                        ) = scidx[
                                                            "arr_data_df_gtf_gene"
                                                        ][
                                                            scidx[
                                                                "dict_index_df_gtf_gene"
                                                            ][
                                                                id_gene_overlapped_with_read
                                                            ][
                                                                0
                                                            ]
                                                        ]  # retrieve information about the gene_id # 1-based
                                                        if (
                                                            strand_anno
                                                            != strand_read
                                                        ):
                                                            dict_id_gene_to_int_num_base_pairs_overlap_with_exons.pop(
                                                                id_gene_overlapped_with_read
                                                            )
                                                """ assign gene_id containing exons with the largest number of overlapped base pairs with the current read """
                                                (
                                                    l_id_gene_assigned,
                                                    int_num_base_pairs_overlap_with_exons,
                                                ) = bk.DICTIONARY_Find_keys_with_max_value(
                                                    dict_id_gene_to_int_num_base_pairs_overlap_with_exons
                                                )  # retrieve the gene_id containing the maximum number of base pairs overlapping with the read

                                                if len(l_id_gene_assigned) == 1:
                                                    """for valid id_gene assignment"""
                                                    int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id = int_num_base_pairs_overlap_with_exons  # update the metric
                                                    id_gene = l_id_gene_assigned[0]
                                                else:
                                                    """when no valid id_gene assignment is available"""
                                                    int_flag_classification ^= int_flag_class_ambiguous_assignment_to_gene

                                """
                                When the current read was assigned to a single gene unambiguously
                                """
                                if (
                                    int_flag_classification
                                    & int_flag_class_overlap_with_gene_body
                                    and not (
                                        int_flag_classification
                                        & int_flag_class_ambiguous_assignment_to_gene
                                    )
                                    and len(id_gene) != 0
                                ):  # when a single valid gene_id was assigned to the read
                                    float_time_start = (
                                        time.time()
                                    )  # record the start time
                                    # identify transcript
                                    (
                                        seqname_gene,
                                        source_gene,
                                        feature_gene,
                                        start_gene,
                                        end_gene,
                                        score_gene,
                                        strand_gene,
                                        frame_gene,
                                        attribute_gene,
                                        name_gene,
                                    ) = scidx["arr_data_df_gtf_gene"][
                                        scidx["dict_index_df_gtf_gene"][id_gene][0]
                                    ]  # retrieve start and end positions and other data of the gene_id # 1-based
                                    start_gene -= 1  # 1-based > 0-based coordinates

                                    """ initialize data for the new gene annotation encountered """
                                    if (
                                        id_gene not in reads["data"]
                                    ):  # if a gene newly encountered
                                        if (
                                            flag_is_mode_scarab_count_atac
                                        ):  # when processing ATAC data, ignore isoform information, and initialize data as 'miscellaneous annotation'
                                            __Initialize_misc_anno_data__(
                                                reads,
                                                refname,
                                                start_gene,
                                                end_gene,
                                                id_gene,
                                                "gene",
                                                str_mode_scarab_count,
                                            )
                                        else:
                                            """initialize gene and isoform data"""
                                            __Initialize_gene_and_isoform_data__(
                                                reads,
                                                refname,
                                                start_gene,
                                                end_gene,
                                                id_gene,
                                            )
                                    """ (GEX mode specific) """
                                    if (
                                        not flag_is_mode_scarab_count_atac
                                    ):  # analyze transcript information, including isoform assignment
                                        """
                                        if 'flag_use_isoform_assignment_from_10x_cellranger' is True, the splice junciton counting and exon counting will not be performed because the transcript annotations used in cellranger and Scarab might be different.
                                        """
                                        """ align current read to known transcript sequences from given fasta sequences and assign id_tx to the read (if isoform assignment from 10x cellranger is not used) """
                                        if (
                                            not flag_use_isoform_assignment_from_10x_cellranger
                                            and seq is not None
                                        ):
                                            ''' 
                                            retrieve alignments to transcripts 
                                            # align the current read to the transcript sequences # exhaust the iterator to avoid the potential memory leakage issue from minimap2 mappy
                                            '''
                                            if flag_filtering_alignment_to_transcript_during_realignment_based_on_structural_difference : # retrieve alignments to transcripts by aligning a sequence excluding softclipped regions
                                                ''' retrieve a read sequence excluding soft-clipped regions, in correct strand '''
                                                int_length_softclipped_left = cigartuples[ 0 ][ 1 ] if int_cigarop_S == cigartuples[ 0 ][ 0 ] else 0
                                                int_length_softclipped_right = cigartuples[ -1 ][ 1 ] if int_cigarop_S == cigartuples[ -1 ][ 0 ] else 0
                                                
                                                seq_excluding_soft_clipping = seq[ int_length_softclipped_left : len( seq ) - int_length_softclipped_right ] # retrieve a read sequence excluding softclipped, reverse complemented back to the original read if the read has been reverse complemented.
                                                seq_excluding_soft_clipping_correct_strand = SEQ.Reverse_Complement( seq_excluding_soft_clipping ) if flag_is_reverse_complemented else seq_excluding_soft_clipping # retrieve sequence of the correct strand
                                                len_seq_excluding_soft_clipping = len( seq_excluding_soft_clipping_correct_strand ) # retrieve length of sequence excluding soft clipping from the genomic alignment
                                                
                                                ''' retrieve alignments to transcripts '''
                                                l_aln_to_tx = list( reads["data"][id_gene]["am_tx"].map( seq_excluding_soft_clipping_correct_strand ) ) # retrieve alignments to the transcripts by aligning the sequence excluding soft-clipped portions, strand-corrected
                                            else :
                                                seq_correct_strand = SEQ.Reverse_Complement( seq ) if flag_is_reverse_complemented else seq # retrieve sequence of the correct strand
                                                l_aln_to_tx = list( reads["data"][id_gene]["am_tx"].map( seq_correct_strand ) ) # retrieve alignments to the transcripts by aligning an entire sequence, strand-corrected
                                                
                                            ''' filter out alignments that were reverse complemented during re-alignment, if read aligned to opposite strand will be excluded. '''
                                            if not flag_include_read_aligned_to_opposite_strand :
                                                l_aln_to_tx = list( hit for hit in l_aln_to_tx if hit.strand > 0 ) # filter out reverse complemented alignments
                                                
                                            if flag_filtering_alignment_to_transcript_during_realignment_based_on_structural_difference : # filter alignment to transcript based on softclipping, insertions, and deletions
                                                ''' filter alignments with extensive softclipping, insertions, and deletions '''
                                                l_aln_to_tx = list( hit for hit in l_aln_to_tx if ( max( hit.q_st, len_seq_excluding_soft_clipping - hit.q_en ) <= int_max_softclipping_and_indel_length_for_filtering_alignment_to_transcript_during_realignment ) and _mappy_detect_structural_difference( hit.cigar ) ) # filter alignments to the transcript if the length of the softclipping, insertions, and deletions exceed the limit
                                            
                                            if flag_enforce_transcript_end_site_matching_for_long_read_during_realignment and not flag_internal_polya_tract_primed :
                                                ''' filter transcript alignments by enforcing transcript end site (TES) matching. transcript alignment with the distance between the alignment end position and the actual end of the transcript larger than the threshold will be filtered, if the read is not internal-polyA-primed. '''
                                                l_aln_to_tx = list( hit for hit in l_aln_to_tx if ( hit.ctg_len - hit.r_en ) <= int_max_distance_from_transcript_end_for_tes_matching_during_realignment ) # filter alignments with transcript end site (TES) matching, and retain transcript alignment 
                                                
                                            dict_assignment_to_score = dict( ( ( hit.ctg, hit.r_st, hit.r_en ), hit.mlen + hit.mapq ) for hit in l_aln_to_tx if hit.mapq >= int_min_mapq_minimap2_tx_assignment ) # filter with mapping quality, # retrieve tx assignment - score mapping, where score is calculated as matched-length + mapping quality, since 'mapping quality' can be 0 value despite its longer aligned length.
                                            if len( dict_assignment_to_score ) > 0 : # if at least one transcript alignment is available
                                                id_tx_assigned_by_minimap2, start_in_transcript, end_in_transcript = _argmax( dict_assignment_to_score ) # find the transcript alignment with the best score (if more than two best scores are available, select one without specific criteria, for now)
                                                    
                                            if ( not flag_skip_exon_and_splice_junc_counting and len( id_tx_assigned_by_minimap2 ) > 0 ):  # if exon and splice_junction counting behavior is enabled and valid 'id_tx_assigned_by_minimap2' has been assigned.
                                                """ retrieve exons of the current transcript overlapped with the current read - using minimap2 realignment result """
                                                for (
                                                    e
                                                ) in __data_object_search_query(
                                                    data_dict_it_exon_transcriptome,
                                                    id_tx_assigned_by_minimap2,
                                                    slice(
                                                        start_in_transcript,
                                                        end_in_transcript,
                                                    ),
                                                ):
                                                    l_id_anno_exon_and_splice_junc.append(
                                                        id_gene
                                                        + "|tx_name="
                                                        + scidx[
                                                            "dict_id_tx_to_name_tx"
                                                        ][
                                                            id_tx_assigned_by_minimap2
                                                        ]
                                                        + "|tx_id="
                                                        + id_tx_assigned_by_minimap2
                                                        + "|exon_id="
                                                        + "{}:{}-{}.{}".format(
                                                            *e[2]
                                                        )
                                                        + "|realigned"
                                                    )
                                                """ retrieve splice junctions of the current transcript overlapped with the current read - using minimap2 realignment result """
                                                for (
                                                    e
                                                ) in __data_object_search_query(
                                                    data_dict_it_splice_junc_transcriptome,
                                                    id_tx_assigned_by_minimap2,
                                                    slice(
                                                        start_in_transcript,
                                                        end_in_transcript,
                                                    ),
                                                ):  # single exon genes lack splice junction annotations
                                                    l_id_anno_exon_and_splice_junc.append(
                                                        id_gene
                                                        + "|tx_name="
                                                        + scidx[
                                                            "dict_id_tx_to_name_tx"
                                                        ][
                                                            id_tx_assigned_by_minimap2
                                                        ]
                                                        + "|tx_id="
                                                        + id_tx_assigned_by_minimap2
                                                        + "|sj_id="
                                                        + "{}:{}-{}.{}".format(
                                                            *e[2]
                                                        )
                                                        + "|realigned"
                                                    )

                                        """ calculate proportion of non-exonic features (=intronic) in the reads and classify reads """
                                        # count bases overlapping exonic features
                                        int_base_gene_exon_count = 0
                                        pos_seg_end_previous_seg = (
                                            -1
                                        )  # initialize the end position of the previous segment (0-based coordinates)

                                        ''' for better performance, search for all segments at once '''
                                        if not (
                                            flag_use_isoform_assignment_from_10x_cellranger
                                            or flag_skip_exon_and_splice_junc_counting
                                        ):  # if exon and splice_junction counting behavior is enabled
                                            l_queries = list( slice(start_seg, end_seg) for start_seg, end_seg in l_seg ) # compose range queries
                                            if res_overlapped_exons is None:
                                                res_overlapped_exons = __data_object_search_queries(
                                                    data_dict_it_exon,
                                                    refname,
                                                    l_queries,
                                                )  # for each segment of the current read, identify overlapped exons

                                            res_overlapped_splice_donor_and_acceptor = __data_object_search_queries(
                                                data_dict_it_splice_donor_and_acceptor,
                                                refname,
                                                l_queries,
                                            )  # for each segment of the current read, identify overlapped splice donor and acceptors
                                        l_t_splice_junc_genome = (
                                            []
                                        )  # collect the list of splice junctions

                                        for i_seg in range(
                                            len(l_seg)
                                        ):  # for each segment
                                            pos_seg_start, pos_seg_end = l_seg[
                                                i_seg
                                            ]  # retrieve positions of the segment
                                            # count the number of exonic bases in the aligned segments
                                            ba_seg_mask_gene_exon = scidx[
                                                "dict_seqname_to_mask_gtf_exon"
                                            ][refname][pos_seg_start:pos_seg_end]
                                            int_base_gene_exon_count += (
                                                ba_seg_mask_gene_exon.count()
                                            )

                                            if not (
                                                flag_use_isoform_assignment_from_10x_cellranger
                                                or flag_skip_exon_and_splice_junc_counting
                                            ):  # if exon and splice_junction counting behavior is enabled
                                                """retrieve exons of the current gene overlapped with the current read - using initial genomic alignment result"""
                                                l_interval = list(
                                                    e
                                                    for e in res_overlapped_exons[
                                                        i_seg
                                                    ]
                                                    if e[2] == id_gene
                                                )  # retrieve exons of the currently assigned 'id_gene' that were overlapped with the current segment
                                                if (
                                                    len(l_interval) == 1
                                                ):  # when counting exon based on genome alignment, assign id_anno_exon only when a unique exon can be assigned to a segment
                                                    e = l_interval[
                                                        0
                                                    ]  # retrieve the interval of the assigned exon
                                                    l_id_anno_exon_and_splice_junc.append(
                                                        id_gene
                                                        + "|exon_id="
                                                        + f"{refname}:{e[ 0 ] + 1}-{e[ 1 ]}.{strand_gene}"
                                                    )  # 1>0-based coordinates

                                                """[Intron Retention Counting] retrieve splice donor and acceptor sites of the current gene overlapped with the current read - using initial genomic alignment result"""
                                                if not flag_skip_intron_retention_counting : # if intron retention events are counted
                                                    for e in res_overlapped_splice_donor_and_acceptor[ i_seg ] : # retrieve splice_donor_and_acceptors that were overlapped with the current segment
                                                        if e[2][0] != id_gene : # splice_donor_and_acceptors belonging to the currently assigned 'id_gene' will be analyzed
                                                            continue
                                                        if ( ( e[ 0 ] + int_min_length_intron_for_detecting_intron_retention_event < pos_seg_end ) if e[2][1] == 'L' else ( e[ 0 ] + int_min_length_intron_for_detecting_intron_retention_event >= pos_seg_start ) ) : # if splice donor and acceptor is located at the left exon (and thus intron is located at the right side), check whether sufficient length of intron is available on the right side to detect intron retention events.
                                                            l_id_anno_exon_and_splice_junc.append( id_gene + "|intron_retention@splice_donor_and_acceptor_id=" + f"{refname}:{e[ 0 ] + 1}.{e[2][1]}.{strand_gene}" )  # 1>0-based coordinates # assign an annotation

                                                """ retrieve splice junctions of the current transcript overlapped with the current read - using initial genomic alignment result """
                                                if (
                                                    pos_seg_end_previous_seg != -1
                                                ):  # when previous segment exists, check whether the splicing junction exists # skip when this segment is the first segment
                                                    l_t_splice_junc_genome.append(
                                                        (
                                                            refname,
                                                            pos_seg_end_previous_seg,
                                                            pos_seg_start,
                                                        )
                                                    )  # collect splice junctions of the read

                                                    flag_skip_intron_retention_counting
                                                pos_seg_end_previous_seg = pos_seg_end  # update the end of segment end coordinates for the next splicing events

                                        if not (
                                            flag_use_isoform_assignment_from_10x_cellranger
                                            or flag_skip_exon_and_splice_junc_counting
                                        ):  # if exon and splice_junction counting behavior is enabled
                                            dict_t_splice_junc_to_info_genome_subset = __data_object_subset(
                                                data_dict_t_splice_junc_to_info_genome,
                                                l_t_splice_junc_genome,
                                            )  # retrieve data of the collected splice junctions
                                            for (
                                                t_splice_junc_genome
                                            ) in dict_t_splice_junc_to_info_genome_subset:  # for each valid splice junction
                                                if (
                                                    len(
                                                        list(
                                                            1
                                                            for id_tx, strand in dict_t_splice_junc_to_info_genome_subset[
                                                                t_splice_junc_genome
                                                            ]
                                                            if scidx[
                                                                "dict_id_tx_to_id_gene"
                                                            ][id_tx]
                                                            == id_gene
                                                        )
                                                    )
                                                    > 0
                                                ):  # if the splice junction belongs to the current gene, assign the read to the 'splice_junc' feature type.
                                                    l_id_anno_exon_and_splice_junc.append(
                                                        id_gene
                                                        + "|sj_id="
                                                        + f"{refname}:{t_splice_junc_genome[ 1 ] + 1}-{t_splice_junc_genome[ 2 ]}.{strand_gene}"
                                                    )  # 0>1-based coordinates

                                        # add classification label based on the proportion of exons in the read
                                        # defualt: 'spanning_both_intron_and_exon'
                                        float_prop_exon = (
                                            int_base_gene_exon_count
                                            / int_total_aligned_length
                                        )  # calculate the proportion of exons contained in the read
                                        if float_prop_exon == 0:
                                            int_flag_classification ^= int_flag_class_completely_intronic  # set 'completely intronic' flag
                                        elif float_prop_exon < 0.90:
                                            pass  # classify 'spanning_both_intron_and_exon'
                                        else:
                                            int_flag_classification ^= int_flag_class_mostly_exonic  # set 'mostly exonic' flag

                                    """ increament the time passed to process the read to the total wall time for the gene """
                                    reads["data"][id_gene]["wall_time"] += (
                                        time.time() - float_time_start
                                    )

                                """
                                check overlap with repeatmasker annotations 
                                """
                                if (
                                    not (
                                        int_flag_classification
                                        & int_flag_class_overlap_with_gene_body
                                    )
                                    or int_flag_classification
                                    & int_flag_class_completely_intronic
                                ):  # when read is not overlapped with a gene (intergenic) or completely intronic, search overlaps with repeatmasker annotations # for atac data, all reads assigned to genes will be skipped from this analysis
                                    float_time_start = (
                                        time.time()
                                    )  # record the start time
                                    # ignore reads if it does not overlap repeatmasker annotation
                                    set_t_interval_overlap_rpmk = set(
                                        tuple(e)
                                        for e in __data_object_search_query(
                                            data_dict_it_rpmk,
                                            refname,
                                            int_pos_of_read_determining_feature_assignment,
                                        )
                                    )  # retrieve set of gene annotations overlapped with the current alignment
                                    if len(set_t_interval_overlap_rpmk) > 0:
                                        int_flag_classification ^= int_flag_class_overlap_with_filtered_rpmk_anno  # update read classification info.
                                        # identify the name of repeatmasker annotation with the largest overlap with the current read
                                        (
                                            l_interval_rpmk,
                                            int_max_num_base_pairs_overlap_with_rpmk,
                                        ) = bk.DICTIONARY_Find_keys_with_max_value(
                                            dict(
                                                (
                                                    e,
                                                    bk.INTERVAL_Overlap(
                                                        [e[0], e[1]],
                                                        l_start_end_read,
                                                        flag_0_based_coordinate_system=True,
                                                    ),
                                                )
                                                for e in set_t_interval_overlap_rpmk
                                            )
                                        )  # retrieve the proportion of overlap with the given feature
                                        if len(l_interval_rpmk) == 1:
                                            (
                                                start_rpmk,
                                                end_rpmk,
                                                id_rpmk,
                                            ) = l_interval_rpmk[
                                                0
                                            ]  # retrieve information about the assigned annotation

                                            """ initialize data for the new rpmk annotation encountered """
                                            if (
                                                id_rpmk not in reads["data"]
                                            ):  # if a annotation is newly encountered
                                                __Initialize_misc_anno_data__(
                                                    reads,
                                                    refname,
                                                    start_rpmk,
                                                    end_rpmk,
                                                    id_rpmk,
                                                    "rpmk",
                                                    str_mode_scarab_count,
                                                )

                                            """ (GEX mode specific) calculate proportion of filtered repeatmasker features of the read """
                                            if not flag_is_mode_scarab_count_atac:
                                                # count bases overlapping repeatmasker features
                                                int_base_filtered_rpmk_count = 0
                                                for (
                                                    pos_seg_start,
                                                    pos_seg_end,
                                                ) in l_seg:
                                                    ba_seg_mask = scidx[
                                                        "dict_seqname_to_mask_gtf_rpmk_filtered"
                                                    ][refname][
                                                        pos_seg_start:pos_seg_end
                                                    ]
                                                    int_base_filtered_rpmk_count += (
                                                        ba_seg_mask.count()
                                                    )

                                                # add classification label based on the proportion of filtered repeatmasker annotations in the read
                                                float_prop_filtered_rpmk = (
                                                    int_base_filtered_rpmk_count
                                                    / int_total_aligned_length
                                                )  # calculate the proportion of exons contained in the read
                                                if float_prop_filtered_rpmk == 1:
                                                    int_flag_classification ^= int_flag_class_complete_overlap_with_filtered_rpmk_anno

                                            """ increament the time passed to process the read to the total wall time for the gene """
                                            reads["data"][id_rpmk]["wall_time"] += (
                                                time.time() - float_time_start
                                            )
                                        else:
                                            int_flag_classification ^= int_flag_class_ambiguous_assignment_to_filtered_rpmk_anno

                                """
                                check overlap with regulatory elements
                                """
                                if (
                                    not (
                                        int_flag_classification
                                        & int_flag_class_overlap_with_gene_body
                                    )
                                    or int_flag_classification
                                    & int_flag_class_completely_intronic
                                ):  # when read is not overlapped with a gene body or is completely intronic, search overlaps with regulatory annotations
                                    float_time_start = (
                                        time.time()
                                    )  # record the start time
                                    # ignore reads if it does not overlap regulatory element annotation
                                    """ calculate proportion of unfiltered repeatmasker elements of the read """
                                    # count bases overlapping regulatory element features
                                    int_base_unfiltered_rpmk_count = 0
                                    for pos_seg_start, pos_seg_end in l_seg:
                                        ba_seg_mask = scidx[
                                            "dict_seqname_to_mask_gtf_rpmk_unfiltered"
                                        ][refname][pos_seg_start:pos_seg_end]
                                        int_base_unfiltered_rpmk_count += (
                                            ba_seg_mask.count()
                                        )

                                    # calculate the proportion of unfiltered repeatmasker annotations in the read
                                    float_prop_unfiltered_rpmk = (
                                        int_base_unfiltered_rpmk_count
                                        / int_total_aligned_length
                                    )

                                    set_t_interval_overlap_reg = set(
                                        tuple(e)
                                        for e in __data_object_search_query(
                                            data_dict_it_reg,
                                            refname,
                                            int_pos_of_read_determining_feature_assignment,
                                        )
                                    )  # retrieve set of gene annotations overlapped with the current alignment
                                    """ if the read overlaps with regulatory element and not overlapped with unfiltered repeatmasker counts """
                                    if len(set_t_interval_overlap_reg) > 0:
                                        int_flag_classification ^= int_flag_class_overlap_with_reg  # update read classification info.
                                        if (
                                            len(id_rpmk) > 0
                                        ):  # if the read has been already overlapped with id_rpmk, update the read classification info.
                                            int_flag_classification ^= int_flag_class_overlap_with_reg_and_rpmk  # update read classification info.
                                        if (
                                            float_prop_unfiltered_rpmk
                                            <= float_max_prop_unfiltered_rpmk
                                        ):
                                            int_flag_classification ^= int_flag_class_overlap_with_reg_not_overlap_with_unfiltered_rpmk_anno  # update read classification info.
                                            # identify the name of repeatmasker annotation with the largest overlap with the current read
                                            (
                                                l_interval_reg,
                                                int_max_num_base_pairs_overlap_with_reg,
                                            ) = bk.DICTIONARY_Find_keys_with_max_value(
                                                dict(
                                                    (
                                                        e,
                                                        bk.INTERVAL_Overlap(
                                                            [e[0], e[1]],
                                                            l_start_end_read,
                                                            flag_0_based_coordinate_system=True,
                                                        ),
                                                    )
                                                    for e in set_t_interval_overlap_reg
                                                )
                                            )  # retrieve the proportion of overlap with the given feature
                                            if len(l_interval_reg) == 1:
                                                (
                                                    start_reg,
                                                    end_reg,
                                                    id_reg,
                                                ) = l_interval_reg[
                                                    0
                                                ]  # retrieve information about the assigned annotation

                                                """ initialize data for the new regulatory element annotation encountered """
                                                if (
                                                    id_reg not in reads["data"]
                                                ):  # if a annotation is newly encountered
                                                    __Initialize_misc_anno_data__(
                                                        reads,
                                                        refname,
                                                        start_reg,
                                                        end_reg,
                                                        id_reg,
                                                        "regulatory_region",
                                                        str_mode_scarab_count,
                                                    )

                                                """ (GEX mode specific) calculate proportion of regulatory element of the read """
                                                if (
                                                    not flag_is_mode_scarab_count_atac
                                                ):
                                                    # count bases overlapping regulatory element features
                                                    int_base_reg_count = 0
                                                    for (
                                                        pos_seg_start,
                                                        pos_seg_end,
                                                    ) in l_seg:
                                                        ba_seg_mask = scidx[
                                                            "dict_seqname_to_mask_gtf_reg"
                                                        ][refname][
                                                            pos_seg_start:pos_seg_end
                                                        ]
                                                        int_base_reg_count += (
                                                            ba_seg_mask.count()
                                                        )

                                                    # add classification label based on the proportion of unfiltered repeatmasker annotations in the read
                                                    float_prop_reg = (
                                                        int_base_reg_count
                                                        / int_total_aligned_length
                                                    )  # calculate the proportion of exons contained in the read
                                                    if float_prop_reg == 1:
                                                        int_flag_classification ^= int_flag_class_complete_overlap_with_reg_not_overlap_with_unfiltered_rpmk_anno

                                                """ increament the time passed to process the read to the total wall time for the gene """
                                                reads["data"][id_reg][
                                                    "wall_time"
                                                ] += (
                                                    time.time() - float_time_start
                                                )
                                            else:
                                                int_flag_classification ^= int_flag_class_ambiguous_assignment_to_reg

                                """
                                Add data extracted from reads to each appropriate annotation
                                """
                                str_umi_corrected, str_umi_uncorrected = (
                                    dict_tags[str_name_bam_tag_umi_corrected],
                                    dict_tags[str_name_bam_tag_umi_uncorrected],
                                )
                                """ collect variant information """
                                l_id_anno_variant = (
                                    []
                                )  # a list of id_anno of the variant annotation type
                                if (
                                    flag_does_not_collect_variant_information
                                    or refname not in scidx["dict_fa_genome"]
                                ):  # does not retrieve variants if reference sequence is not available
                                    str_l_var = ""
                                else:
                                    l_var_name = SAM.Call_Variant(
                                        r,
                                        scidx["dict_fa_genome"],
                                        function_for_processing_reference_name=__chromosome_name_remove_chr__,
                                    )
                                    if (
                                        flag_filter_variant_using_predefined_set
                                    ):  # filter variants using predefined set
                                        l_var_name = list(
                                            e
                                            for e in l_var_name
                                            if e in set_var_name_valid
                                        )
                                    str_l_var = ";".join(
                                        l_var_name
                                    )  # compose 'str_l_var'

                                    """
                                    'variant' feature type
                                    """
                                    if (
                                        flag_filter_variant_using_predefined_set
                                        and refname
                                        in dict_it_pos_variant_predefined_set
                                    ):  # when filtering variants using predefined set, add 'variant' feature type
                                        # collect overlapped positions of the predefined set of variants
                                        set_overlapped_pos = set()
                                        for st, en in l_seg:
                                            set_overlapped_pos |= (
                                                dict_it_pos_variant_predefined_set[
                                                    refname
                                                ][st:en]
                                            )
                                        for (
                                            start,
                                            end,
                                            values,
                                        ) in (
                                            set_overlapped_pos
                                        ):  # for each overlapped position of predefined set of variants
                                            id_anno = f"variant|pos={refname}:{start + 1}"  # compose 'id_anno' of the reference position
                                            l_id_anno_variant.append(
                                                id_anno
                                            )  # collect 'id_anno'
                                            if (
                                                id_anno not in reads["data"]
                                            ):  # initialize the 'variant' feature
                                                __Initialize_misc_anno_data__(
                                                    reads,
                                                    refname,
                                                    start,
                                                    end,
                                                    id_anno,
                                                    "variant",
                                                    str_mode_scarab_count,
                                                )  # for variant annotation, start and end position of the annotation will be set as the reference position of the variant

                                """ 
                                structure of 'l_data' : 
                                front) existing annotation from the input BAM file
                                rear) new annotations added by the current program
                                """
                                if flag_is_mode_scarab_count_atac:
                                    l_data = [
                                        qname,
                                        int_mapq,
                                        refname,
                                        flag,
                                        refstart + 1,
                                        refend,
                                        str_l_seg,
                                        dict_tags[str_name_bam_tag_cb_corrected],
                                        dict_tags[str_name_bam_tag_cb_uncorrected],
                                        int_flag_classification,
                                        id_rpmk,
                                        id_gene,
                                        id_reg,
                                        id_promoter,
                                        str_l_var,
                                        int_total_aligned_length,
                                    ]  # refstart in 1-based coordinate
                                else:
                                    l_data = [
                                        qname,
                                        int_mapq,
                                        refname,
                                        flag,
                                        refstart + 1,
                                        refend,
                                        str_l_seg,
                                        dict_tags[str_name_bam_tag_cb_corrected],
                                        dict_tags[str_name_bam_tag_cb_uncorrected],
                                        str_umi_corrected
                                        if len(str_umi_corrected) > 0
                                        else str_umi_uncorrected,
                                        str_umi_uncorrected,
                                        dict_tags["TX"],
                                        dict_tags["AN"],
                                        dict_tags["GX"],
                                        dict_tags["GN"],
                                        dict_tags["MM"],
                                        dict_tags["RE"],
                                        dict_tags["xf"],
                                        int_flag_classification,
                                        id_rpmk,
                                        int_max_num_base_pairs_overlap_with_rpmk,
                                        id_gene,
                                        int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id,
                                        int_base_gene_exon_count,
                                        int_base_filtered_rpmk_count,
                                        id_reg,
                                        int_base_unfiltered_rpmk_count,
                                        int_base_reg_count,
                                        id_tx_assigned_by_minimap2,
                                        str_l_var,
                                        int_total_aligned_length,
                                        flag_internal_polya_tract_primed,
                                    ]  # refstart in 1-based coordinate # use 'str_umi_uncorrected' if 'str_umi_corrected' does not exist
                                l_data_for_counting = list(
                                    l_data[index_col]
                                    for index_col in l_index_col_for_counting
                                )  # retrieve partial data for counting

                                """ initialize feature data object for exon and splice junctions """
                                for (
                                    id_anno
                                ) in (
                                    l_id_anno_exon_and_splice_junc
                                ):  # for ATAC data, 'l_id_anno_exon_and_splice_junc' should be empty
                                    if (
                                        id_anno not in reads["data"]
                                    ):  # if the exon or splice_juc anno has not been initialized, initialize the feature
                                        __Initialize_misc_anno_data__(
                                            reads,
                                            refname,
                                            start_gene,
                                            end_gene,
                                            id_anno,
                                            "exon_and_splice_junc",
                                            str_mode_scarab_count,
                                        )  # for exon and splice_junc annotation, start and end position of the annotation will be set as that of the gene to which read has been assigned to.
                                l_id_anno_valid = list(
                                    id_anno
                                    for id_anno in [
                                        id_gene,
                                        id_rpmk,
                                        id_reg,
                                        id_promoter,
                                    ]
                                    + l_id_anno_exon_and_splice_junc
                                    + l_id_anno_variant
                                    if len(id_anno) > 0
                                )  # retrieve a list of valid 'id_anno'

                                """ append data to each feature """
                                for id_anno in l_id_anno_valid:
                                    if len(id_anno) > 0:  # if 'id_anno' is valid
                                        reads["data"][id_anno]["l_read"].append(
                                            l_data_for_counting
                                        )

                                """
                                - A Wildcard - 
                                Catch-All (binning for each genomic region)
                                """
                                if (
                                    not flag_turn_off_catching_all_reads_by_binning
                                    and (
                                        len(l_id_anno_valid) == 0
                                        or not flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning
                                    )
                                ):  # if conditions is met # if 'flag_exclude_reads_assigned_to_features_and_count_every_read_in_bam_during_binning' is False, collect the read in a genomic bin regardless of whether the read has been assigned to existing features.
                                    """retrieve genomie region (bin)"""
                                    (
                                        start_bin_genome,
                                        end_bin_genome,
                                    ) = __Get_Genomic_Region__(
                                        int_pos_of_read_determining_feature_assignment
                                    )  # retrieve 0-based coordinates of the current genomic region (bin)
                                    id_bin_genome = f"genomic_region|pos={refname}:{start_bin_genome + 1}-{end_bin_genome}"  # 0>1 based coordinates

                                    """ initialize data for the current genomic region encountered """
                                    if (
                                        id_bin_genome not in reads["data"]
                                    ):  # if a annotation is newly encountered
                                        __Initialize_misc_anno_data__(
                                            reads,
                                            refname,
                                            start_bin_genome,
                                            end_bin_genome,
                                            id_bin_genome,
                                        )
                                    reads["data"][id_bin_genome]["l_read"].append(
                                        l_data_for_counting
                                    )

                                """ write the analysis result to the output BAM file """
                                if (
                                    not flag_skip_read_analysis_summary_output_bam_file
                                ):
                                    """
                                    write annotated record as a tublar or SAM record
                                    """
                                    """
                                    SAM record structure:
                                    (index : Field )
                                    0      : QNAME
                                    1      : FLAG
                                    2      : RNAME
                                    3      : POS
                                    4      : MAPQ
                                    5      : CIGAR
                                    6      : RNEXT
                                    7      : PNEXT
                                    8      : TLEN
                                    9      : SEQ
                                    10     : QUAL
                                    """
                                    """ retrieve sam record as a string """
                                    l_sam = r.tostring().split("\t")
                                    """ remove chr prefix from the chromosome names (if present) """
                                    l_sam[2] = __chromosome_name_remove_chr__(
                                        l_sam[2]
                                    )
                                    l_sam[6] = __chromosome_name_remove_chr__(
                                        l_sam[6]
                                    )
                                    """ delete sequence and quality scores if 'flag_does_not_delete_sequence_and_sequence_qual' is False """
                                    if (
                                        not flag_does_not_delete_sequence_and_sequence_qual
                                    ):
                                        l_sam[9] = "*"
                                        l_sam[10] = "*"

                                    for name_col, data in zip(
                                        l_name_col_newanno,
                                        l_data[-len(l_name_col_newanno) :],
                                    ):
                                        (
                                            name_tag,
                                            str_type,
                                        ) = dict_name_col_newanno_to_sam_tag_name[
                                            name_col
                                        ]
                                        # retrieve string representation of data
                                        str_data = (
                                            str(int(data))
                                            if str_type == "i"
                                            and isinstance(data, int)
                                            else data
                                        )
                                        # only add tag if data contains a valid value
                                        if len(str_data) > 0:
                                            l_sam.append(
                                                name_tag
                                                + ":"
                                                + str_type
                                                + ":"
                                                + str_data
                                            )  # add a tag value
                                    newsamfile.write(
                                        pysam.AlignedSegment.fromstring(
                                            "\t".join(l_sam), samfile_header
                                        )
                                    )  # write a record to the output bam file

                                """ write the analysis result to the output TSV file """
                                if (
                                    not flag_skip_read_analysis_summary_output_tsv_file
                                ):
                                    newfile.write(
                                        (
                                            "\t".join(list(map(str, l_data))) + "\n"
                                        ).encode()
                                    )  # write a record to the output tabular text file

                        """ Flush all remaining data once all iteration has been completed """
                        """
                        'Flush' all remaining data
                        """
                        for id_anno in list(reads["data"]):
                            _empty_bucket( id_anno, ) # flush the data
                        _flush_the_current_batch_of_buckets( ) # flush the last batch of the buckets
                        _write_results_from_offloaded_works( flag_wait_all = True ) # wait for all works to be completed, and flush results from offloaded computations
                            
                        pipe_sender.send(
                            int_n_sam_record_count
                        )  # report the number of processed sam records

                    ''' close output files '''
                    if not flag_skip_read_analysis_summary_output_bam_file:
                        newsamfile.close()
                    if not flag_skip_read_analysis_summary_output_tsv_file:
                        newfile.close()
                    for e in dict_t_distribution_range_of_interest_to_newfile_df_count :
                        dict_t_distribution_range_of_interest_to_newfile_df_count[ e ].close()
                    newfile_df_analysis_statistics.close()
                    pipe_sender.send( 'completed' ) # report the worker has completed all works
                    if verbose:
                        logger.info(f"[Completed] ({str_uuid})")

                ns = dict()  # define a namespace
                ns[
                    "int_num_read_currently_processed"
                ] = 0  # initialize total number of reads processed by the algorithm

                def post_process_batch(res):
                    # parse received result
                    int_n_sam_record_count = res
                    ns["int_num_read_currently_processed"] += int_n_sam_record_count
                    logger.info(
                        f"[{path_file_bam_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed."
                    )  # report

                """
                Analyze a Barcoded BAM file
                """
                if verbose:
                    logger.info(
                        f"[{path_file_bam_input}] the analysis pipeline will be run with {n_threads} number of threads"
                    )
                bk.Multiprocessing_Batch_Generator_and_Workers(
                    gen_batch=iter( set( SAM.Get_contig_names_from_bam_header( path_file_bam_input ) ).difference( set_seqname_to_skip ) ), # analyze the pre-processed BAM file for each chromosome # exclude the chromosomes in the given list of sequence names to exclude in the analysis
                    process_batch=process_batch,
                    post_process_batch=post_process_batch,
                    int_num_threads=n_threads
                    + 2,  # one thread for generating batch, another thread for post-processing of the batch
                )
                
                """ 
                Export Gene, Transcript, and Misc. Annotation Read Count Matrix 
                """

                def export_count_matrix():  # off-loading a single-core work
                    logger.info(
                        f"[{path_file_bam_input}] Exporting count matrix started."
                    )

                    ''' combine count files '''
                    l_col_df_count = ["barcode", "feature", "id_feature", "read_count"]
                    for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest'
                        str_suffix = f"count.size_distribution__{_t_distribution_range_of_interest_to_str( t_distribution_range_of_interest )}.tsv.gz" # compose the suffix of the files belonging to the 't_distribution_range_of_interest'
                        bk.OS_FILE_Combine_Files_in_order(
                            glob.glob(
                                f"{path_folder_temp}*.{str_suffix}"
                            ), # retrieve the list of input file paths
                            f"{path_folder_output}df_{str_suffix}",
                            header="\t".join(l_col_df_count) + "\n",
                            delete_input_files=True,
                            overwrite_existing_file=True,
                        )
                    ''' combine other files '''
                    l_col_df_stat = ["str_uuid", "id_anno", "wall_time", "int_n_reads"]
                    bk.OS_FILE_Combine_Files_in_order(
                        glob.glob(f"{path_folder_temp}*.analysis_statistics.tsv.gz"),
                        f"{path_folder_output}analysis_statistics.tsv.gz",
                        header="\t".join(l_col_df_stat) + "\n",
                        delete_input_files=True,
                        overwrite_existing_file=True,
                    )
                    # combine results into a single output file (initial read analysis)
                    if not flag_skip_read_analysis_summary_output_tsv_file:
                        if str_mode_scarab_count == "atac":
                            l_col = [
                                "qname",
                                "mapq",
                                "refname",
                                "flag",
                                "refstart",
                                "refend",
                                "str_l_seg",
                                "CB",
                                "CR",
                                "int_flag_classification",
                                "id_rpmk",
                                "id_gene",
                                "id_reg",
                                "id_promoter",
                                "l_name_variant",
                                'int_total_aligned_length',
                            ]
                        else:
                            l_col = [
                                "qname",
                                "mapq",
                                "refname",
                                "flag",
                                "refstart",
                                "refend",
                                "str_l_seg",
                                "CB",
                                "CR",
                                "UB",
                                "UR",
                                "TX",
                                "AN",
                                "GX",
                                "GN",
                                "MM",
                                "RE",
                                "xf",
                                "int_flag_classification",
                                "id_rpmk",
                                "int_max_num_base_pairs_overlap_with_rpmk",
                                "id_gene",
                                "int_num_base_pairs_overlap_with_exons_of_the_assigned_gene_id",
                                "int_base_gene_exon_count",
                                "int_base_filtered_rpmk_count",
                                "id_reg",
                                "int_base_unfiltered_rpmk_count",
                                "int_base_reg_count",
                                "id_tx_assigned_by_minimap2",
                                "l_name_variant",
                                'int_total_aligned_length',
                                'flag_internal_polya_tract_primed',
                            ]
                            
                        bk.OS_FILE_Combine_Files_in_order(
                            glob.glob(f"{path_folder_temp}*.analysis.*.tsv.gz"),
                            f"{path_folder_output}analysis.tsv.gz",
                            header="\t".join(l_col) + "\n",
                            delete_input_files=True,
                            overwrite_existing_file=True,
                        )

                    ''' export_count_matrix '''
                    for t_distribution_range_of_interest in l_t_distribution_range_of_interest : # for each 't_distribution_range_of_interest'
                        flag_output_dtype_is_integer = t_distribution_range_of_interest is None # output dtype is integer if the raw count matrix is exported.
                        str_t_distribution_range_of_interest = _t_distribution_range_of_interest_to_str( t_distribution_range_of_interest ) # retrieve string representation of 't_distribution_range_of_interest'
                        str_suffix = f"count.size_distribution__{str_t_distribution_range_of_interest}.tsv.gz" # compose the suffix of the files belonging to the 't_distribution_range_of_interest'
                        path_folder_mtx = f"{path_folder_output}mtx/{str_t_distribution_range_of_interest}/" # compose the output folder
                        os.makedirs( path_folder_mtx, exist_ok = True ) # create the output folder
                        Convert_df_count_to_MTX_10X(
                            path_file_df_count=f"{path_folder_output}df_{str_suffix}",
                            path_folder_mtx_10x_output=f"{path_folder_mtx}raw_feature_bc_matrix/",
                            path_folder_mtx_10x_filtered_output=f"{path_folder_mtx}filtered_feature_bc_matrix/",
                            chunksize=1000000,
                            int_min_count_features_for_filtering_barcodes=int_min_count_features_for_filtering_barcodes,
                            flag_output_dtype_is_integer = flag_output_dtype_is_integer,
                        )  # export count matrix as a 10X MTX object
                        
                    # write a flag indicating the count matrix is exported
                    os.mknod(f"{path_folder_output}count_matrix.export_completed.txt")
                    release_lock()  # release the lock
                    logger.info(
                        f"[{path_file_bam_input}] Count matrix in 10X MTX format was exported."
                    )

                workers.submit_work(export_count_matrix)

                """ 
                Export BAM Output
                """
                # combine the output bam files
                if not flag_skip_read_analysis_summary_output_bam_file:
                    def output_bam_file():  # off-loading a single-core work
                        logger.info(
                            f"[{path_file_bam_input}] Combining BAM Output started."
                        )
                        l_path_file_bam_to_be_merged = glob.glob(
                            f"{path_folder_temp}/*.analysis.*.bam"
                        )
                        pysam.merge(
                            f"{path_folder_output}analysis.bam",
                            *l_path_file_bam_to_be_merged,
                            "--threads",
                            str(n_threads),
                        )
                        for path_file_bam in l_path_file_bam_to_be_merged:
                            os.remove(path_file_bam)
                        # index the output bam file
                        bk.OS_Run(
                            ["samtools", "index", f"{path_folder_output}analysis.bam"]
                        )
                        os.mknod(f"{path_folder_output}bam_output.export_completed.txt")
                        release_lock()  # release the lock
                        logger.info(
                            f"[{path_file_bam_input}] BAM Output has been completed."
                        )

                    workers.submit_work(output_bam_file)

            """ for multiome sample, rename feature names (adding '|mode=atac' suffix at the end of feature name) and combine count matrices """
            if str_mode_scarab_count_for_the_current_sample == "multiome":

                def combine_count_matrices_for_multiome():  # off-loading a single-core work
                    path_folder_mtx_gex = f"{path_folder_output_for_the_current_sample}gex/filtered_feature_bc_matrix/"
                    path_folder_mtx_atac = f"{path_folder_output_for_the_current_sample}atac/filtered_feature_bc_matrix/"
                    path_folder_mtx_multiome = f"{path_folder_output_for_the_current_sample}filtered_feature_bc_matrix/"

                    """ if an multiome count matrix folder already exists, check whether the output is complete """
                    if os.path.exists(path_folder_mtx_multiome):
                        """check whether the pipeline has been completed"""
                        if os.path.exists(
                            f"{path_folder_mtx_multiome}matrix.mtx.gz"
                        ) and (
                            not os.path.exists(f"{path_folder_mtx_multiome}matrix.mtx")
                        ):
                            logger.warning(
                                f"[Output Matrix Already Exists] the output folder {path_folder_mtx_multiome} contains a valid count matrix file, skipping"
                            )
                            return  # skip if the pipeline has been completed for the output folder
                        else:
                            """if required output files does not exist or the an intermediate file exists, remove the entire output folder, and rerun the process"""
                            logger.warning(
                                f"[Output Matrix Already Exists] the output folder {path_folder_mtx_multiome} does not contain a valid count matrix file. The output folder will be deleted and the pipeline will be continued from there."
                            )
                            shutil.rmtree(path_folder_mtx_multiome)

                    # add suffix to the ATAC output to distinguish them from GEX outputs
                    SC.MTX_10X_Feature_add_prefix_or_suffix(
                        f"{path_folder_mtx_atac}features.tsv.gz",
                        id_feature_suffix="|mode=atac",
                        name_feature_suffix="|mode=atac",
                    )
                    # combine matrices
                    SC.MTX_10X_Combine(
                        path_folder_mtx_multiome,
                        path_folder_mtx_gex,
                        path_folder_mtx_atac,
                    )
                    os.mknod(
                        f"{path_folder_output_for_the_current_sample}count_matrix.export_completed.txt"
                    )
                    release_lock()  # release the lock
                    logger.info( f'[Combined Count Matrix] A combined count matrix (GEX + ATAC) for "{l_path_file_bam_input_for_the_current_sample}" was exported.' )

                workers.submit_work(combine_count_matrices_for_multiome)

            release_lock()  # release the lock

        # wait all the single-core works offloaded to the workers to be completed.
        workers.wait_all()
        logger.info(
            f"[Pipeline Completion] Forked Pipeline (id={str_uuid_pipeline}) Completed."
        )

    for _ in range(
        int_num_samples_analyzed_concurrently
    ):  # run 'int_num_samples_analyzed_concurrently' number of pipelines
        pipelines.submit_work(run_pipeline)

    # wait all pipelines to be completed
    pipelines.wait_all()
    logger.info(f"[Scarab-Count Start] Program Completed.")
    return scidx  # return the loaded index object

def ourotools(str_mode=None, **dict_args):
    """
    Package settings
    """
    name_package = "ourotools"
    path_remote = "https://github.com/ahs2202/ouro/raw/main/ouro-tools/"  # remote directory from which datafiles will be downloaded
    path_folder_ouro = f"{pkg_resources.resource_filename( name_package, '' )}/"  # directory of the current installed package

    """ check whether the program is called from the command-line interface or from an interactive Python programming environment """
    str_name_program = sys.argv[0]
    if "/" in str_name_program:
        str_name_program = str_name_program.rsplit("/", 1)[1]
    flag_usage_from_command_line_interface = (
        str_name_program[: len("ourotools")] == "ourotools"
    )
    if flag_usage_from_command_line_interface:
        str_mode = sys.argv[1]
        if str_mode == "LongFilterNSplit":
            LongFilterNSplit(flag_usage_from_command_line_interface=True)
        elif str_mode == "LongExtractBarcodeFromBAM":
            LongExtractBarcodeFromBAM(flag_usage_from_command_line_interface=True)
    else:
        if str_mode == "LongFilterNSplit":
            LongFilterNSplit(**dict_args)
        elif str_mode == "LongExtractBarcodeFromBAM":
            LongExtractBarcodeFromBAM(**dict_args)

""" functions that are currently not supported in the command line (only accessible in python environment, including jupyter notebook) """
# functions independent of the main ourotools pipeline (utility functions)
def SplitBAM( path_file_bam_input : str, path_folder_output : str, dict_cb_to_name_clus : dict, name_tag_cb : str = 'CB' ) :
    """ # 2023-09-07 22:18:43 
    A standalone pipeline employing multiprocessing for faster splitting of a barcoded BAM file to multiple BAM files, each containing records of cells belonging to a single 'name_cluster' (a single cell type)
    
    path_file_bam_input : str # an input Barcoded BAM file to split
    path_folder_output : str # the output folder where splitted Barcoded BAM files will be exported
    dict_cb_to_name_clus : dict # a dictionary containing corrected cell barcode to 'name_clus' (name of cluster) mapping. if the name is not compatible with Windoe OS path, incompatible characters will be replaced with spaceholder character.
    name_tag_cb : str = 'CB' # name of the SAM tag containing the corrected cell barcode
    
    """
    # import packages
    import multiprocessing as mp
    import pysam
    import os
    # define functions
    def To_window_path_compatible_str(a_string):
        """
            replace following characters to '_' so that a given string will be compatible for Window file system :
        : (colon)    " (double quote)    / (forward slash)    \ (backslash)    | (vertical bar or pipe)    ? (question mark)    * (asterisk)
            Also, replace new line character into '_'
        """
        return a_string.replace("\n", "_").replace(":", "_").replace('"', "_").replace("/", "_").replace("\\", "_").replace("|", "_").replace("?", "_").replace("*", "_")
    
    # convert 'name_cluster' to window-compatible file name
    dict_convert = dict( ( v, To_window_path_compatible_str( v ).replace( ' ', '_' ) ) for v in set( dict_cb_to_name_clus.values( ) ) )
    dict_cb_to_name_clus = dict( ( e, dict_convert[ dict_cb_to_name_clus[ e ] ] ) for e in dict_cb_to_name_clus )
    l_name_clus = list( dict_convert.values( ) ) # retrieve list of 'name_clus' (after correction)

    # create the output folder
    os.makedirs( path_folder_output, exist_ok = True )
    
    # read the header of the input BAM file    
    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
        sam_header = samfile.header

    def _write_bam_of_name_clus( p_in, p_out ) :
        ''' # 2023-09-07 21:38:10 
        for writing bam file for each 'name_clus'
        '''
        name_clus = p_in.recv( ) # retrieve the name of the cluster
        path_file_bam = f'{path_folder_output}{name_clus}.bam'
        with pysam.AlignmentFile( path_file_bam, 'wb', header = sam_header ) as newsamfile :
            while True :
                batch = p_in.recv( ) # receive a record
                if batch is None :
                    break
                for str_r in batch : # parse the batch
                    r = pysam.AlignedSegment.fromstring( str_r, sam_header ) # compose a pysam record
                    newsamfile.write( r )
        pysam.index( path_file_bam ) # index the given bam file
        p_out.send( 'completed' ) # indicate the work has been completed

    dict_name_clus_p_to_writers = dict( )
    l_p_from_writers = [ ]
    l_process = [ ] # list of processes
    for name_clus in l_name_clus : # for each 'name_clus', recruite a worker
        pm2w, pw4m = mp.Pipe( )
        pw2m, pm4w = mp.Pipe( )
        p = mp.Process( target = _write_bam_of_name_clus, args = ( pw4m, pw2m ) )
        dict_name_clus_p_to_writers[ name_clus ] = pm2w 
        l_p_from_writers.append( pm4w )
        p.start( )
        pm2w.send( name_clus ) # initialize the worker with 'name_clus'
        l_process.append( p ) # collect the process

    # internal setting
    int_max_num_record_in_a_batch = 100
    dict_name_clus_to_batch = dict( ( name_clus, [ ] ) for name_clus in l_name_clus ) # initialize 'dict_name_clus_to_batch'

    def _flush_batch( name_clus : str ) :
        """ # 2023-09-07 22:00:29 
        flush the batch
        """
        batch = dict_name_clus_to_batch[ name_clus ]
        if len( batch ) > 0 : # if batch is not empty
            dict_name_clus_p_to_writers[ name_clus ].send( batch ) # send the batch to the writer
            dict_name_clus_to_batch[ name_clus ] = [ ] # empty the batch    

    def _write_record( name_clus : str, str_r : str ) :
        """ # 2023-09-07 22:00:38 
        write a record
        """
        dict_name_clus_to_batch[ name_clus ].append( str_r ) # add the record
        if len( dict_name_clus_to_batch[ name_clus ] ) >= int_max_num_record_in_a_batch :
            _flush_batch( name_clus ) # flush the batch

    # read file and write the record
    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            dict_tags = dict( r.get_tags( ) ) # retrieve tags of the read
            if name_tag_cb in dict_tags :
                str_cb = dict_tags[ name_tag_cb ]
                if str_cb in dict_cb_to_name_clus : 
                    name_clus = dict_cb_to_name_clus[ str_cb ] # retrieve name of the cluster
                    str_r = r.tostring( ) # convert samtools record to a string
                    _write_record( name_clus, str_r ) # write the record

    # flush remaining records
    for name_clus in l_name_clus :
        _flush_batch( name_clus )

    # notify all works in the main process has been completed
    for name_clus in l_name_clus :
        dict_name_clus_p_to_writers[ name_clus ].send( None )

    # wait for all workers to complete their jobs
    for p in l_p_from_writers :
        p.recv( ) # receive a signal indicating the worker has dismissed itself
    # pipeline completed
    return

def SplitBAMs( dict__path_file_bam_input__to__dict_cb_to_name_clus : dict, path_folder_output : str, name_tag_cb : str = 'CB', int_num_threads : int = 5 ) :
    """ # 2023-09-11 16:52:39 
    A pipeline employing multiprocessing for faster splitting of barcoded BAM files of a single cell dataset to multiple BAM files, each containing records of cells belonging to a single 'name_cluster' (a single cell type)
    
    dict__path_file_bam_input__to__dict_cb_to_name_clus : dict # a dictionary with key = 'path_file_bam_input' and value = 'dict_cb_to_name_clus'.
        * dict_cb_to_name_clus : dict # a dictionary containing corrected cell barcode to 'name_clus' (name of cluster) mapping. if the name is not compatible with Windoe OS path, incompatible characters will be replaced with spaceholder character.
        * path_file_bam_input : str # an input BAM file to split
        
    path_folder_output : str # the output folder where splitted BAM files will be exported
    
    name_tag_cb : str = 'CB' # name of the SAM tag containing the corrected cell barcode
    
    int_num_threads : int = 5 # the number of threads for merging BAM files of the same cluster name (cell type)
    """
    # create the output folder
    os.makedirs( path_folder_output, exist_ok = True )
    
    ''' Initiate pipelines for processing individual BAM file separately. '''
    logger.info(f"Started.")
    pipelines = bk.Offload_Works( None )  # no limit for the number of works that can be submitted.
    
    ''' run 'SplitBAM' for individual BAM files '''
    for path_file_bam_input in dict__path_file_bam_input__to__dict_cb_to_name_clus : # run 'SplitBAM' for each input BAM file separately
        str_uuid_file_bam_input = bk.UUID( ) # retrieve UUID of the input bam file
        pipelines.submit_work( SplitBAM, kwargs = { 'path_file_bam_input' : path_file_bam_input, 'path_folder_output' : f'{path_folder_output}temp_{str_uuid_file_bam_input}/', 'dict_cb_to_name_clus' : dict__path_file_bam_input__to__dict_cb_to_name_clus[ path_file_bam_input ], 'name_tag_cb' : name_tag_cb } )
                    
    # wait all pipelines to be completed
    pipelines.wait_all()
    
    ''' combine temporary output BAM files '''
    df_file_bam = bk.GLOB_Retrive_Strings_in_Wildcards( f'{path_folder_output}*/*.bam' ) # retrieve a list of temporary output BAM files
    for name_clus, _df in df_file_bam[ [ 'wildcard_1', 'path' ] ].groupby( 'wildcard_1' ) : # for each 'name_clus'
        l_path_file_output_temp = _df.path.values # retrieve list of temporary output files
        path_file_bam_output = f"{path_folder_output}{name_clus}.bam"
        pysam.merge( '--threads', str( min( int_num_threads, 10 ) ), '-c', '-p', path_file_bam_output, * l_path_file_output_temp ) # merge splitted BAM files into a single output BAM file
        for path_file in l_path_file_output_temp : # delete the temporary output files
            os.remove( path_file )
        pysam.index( path_file_bam_output ) # index the output file
        
    # remove temporary output files
    for path_folder in bk.GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_output}temp_*/" ).path.values :
        shutil.rmtree( path_folder )
        
    logger.info(f"Completed.")

def DeduplicateBAM( 
    path_file_bam_input : str, # an input Barcoded BAM file to split
    path_folder_output : str, # the output folder where splitted BAM files will be exported
    name_tag_cb : str = 'CB', 
    name_tag_umi : str = 'UB',
    name_tag_length : str = 'LE', # length of molecule aligned to the genome
    int_num_processes : int = 8, # number of processes to use for processing each chunk
    int_num_threads_for_sorting_bam : int = 5, # the number of threads for sorting the BAM file 
) :
    """ # 2023-09-19 00:50:48 
    A standalone pipeline for de-duplicating a Barcoded BAM file based on the CB-UMI attachment genomic position and CB-UMI pairs, generated by ourotools.LongExtractBarcodeFromBAM method. Assumes reads are aligned so that read alignment provides strand-specific information, and alignment direction will be used to infer CB-UMI attachment genomic position.
    
    path_file_bam_input : str, # an input Barcoded BAM file to split
    path_folder_output : str, # the output folder where de-duplicated BAM file will be exported
    name_tag_cb : str = 'CB', 
    name_tag_umi : str = 'UB',
    name_tag_length : str = 'LE', # length of molecule aligned to the genome
    int_num_processes : int = 8, # number of processes to use for processing each chunk
    int_num_threads_for_sorting_bam : int = 5, # the number of threads for sorting the BAM file 
    
    """
    ''' initialize '''
    # define folders
    path_folder_temp = f'{path_folder_output}temp/'

    # create the output folder
    for path_folder in [ path_folder_temp, path_folder_output ] :
        os.makedirs( path_folder, exist_ok = True ) # create the parent folder where the output BAM file will reside

    # read the header of the input BAM file    
    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
        sam_header = samfile.header

    # internal settings
    int_max_num_bucket_deleted = 100000

    def _check_binary_flags( flags : int, int_bit_flag_position : int ) :
        """ # 2023-08-08 22:47:02 
        check a flag in the binary flags at the given position
        """
        return ( flags & ( 1 << int_bit_flag_position ) ) > 0 

    def _generate_batch( ) :
        """ # 2023-09-15 20:31:45 
        generate batch from the input BAM file
        """
        ns = dict( ) # create a namespace
        ns[ 'int_num_buckets_deleted' ] = 0 # initialize 'int_num_buckets_deleted'
        ns[ 'dict_t_id_to_bucket' ] = dict( ) # a dictionary containing batches
        reference_name_current = None
        reference_start_current = None

        def _flush_bucket( t_id ) :
            """ # 2023-09-19 00:27:31 
            """
            bucket = ns[ 'dict_t_id_to_bucket' ].pop( t_id )
            ns[ 'int_num_buckets_deleted' ] += 1
            if ns[ 'int_num_buckets_deleted' ] >= int_max_num_bucket_deleted : # if the number of pop operations exceed the limit, recreate the dictionary
                data = ns[ 'dict_t_id_to_bucket' ]
                ns[ 'dict_t_id_to_bucket' ] = dict( ( k, data[ k ] ) for k in data )
            return bucket

        def _add_record( t_id, r ) :
            """ # 2023-09-19 00:27:25 
            """
            dict_tags = dict( r.get_tags( ) )
            if name_tag_cb not in dict_tags or name_tag_umi not in dict_tags or name_tag_length not in dict_tags : # ignore invalid reads
                return

            if t_id not in ns[ 'dict_t_id_to_bucket' ] : # initialize the bucket for 't_id'
                ns[ 'dict_t_id_to_bucket' ][ t_id ] = [ ]
            ns[ 'dict_t_id_to_bucket' ][ t_id ].append( [ r.to_string( ), dict_tags[ name_tag_cb ], dict_tags[ name_tag_umi ], dict_tags[ name_tag_length ] ] )

        # read file and write the record
        with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
            for r in samfile.fetch( ) :
                # check whether the read was reverse complemented
                flags, reference_name, reference_start, reference_end = r.flag, r.reference_name, r.reference_start, r.reference_end # retrieve read properties

                ''' process reads for each 'bucket' (reads with the same poly CB-UMI attachment sites) '''
                ''' when the contig has changed, empty all buckets '''
                if reference_name_current != reference_name :
                    for t_id in list( ns[ 'dict_t_id_to_bucket' ] ) : # retrieve list of 't_id'
                        yield _flush_bucket( t_id )
                    reference_name_current = reference_name # update 'reference_name_current'

                ''' when the position has changed, detect buckets that should be emptied '''
                if reference_start_current != reference_start :
                    ''' determine whether to empty bucket or not, based on the current position on the sorted BAM file '''
                    for t_id in list( ns[ 'dict_t_id_to_bucket' ] ) : # retrieve list of 't_id'
                        ''' regardlesss of whether CB-UMI attachment site is located at the left or the right side of the read, when the current position passes the poly A site, process the bucket '''
                        flag_is_reverse_complemented, pos = t_id # parse 't_id'
                        if pos < reference_start :
                            yield _flush_bucket( t_id )
                    reference_start_current = reference_start # update 'reference_start_current'

                ''' compose 't_id' '''
                flag_is_reverse_complemented = _check_binary_flags( flags, 4 ) # retrieve a flag indicating whether the read has been reverse-complemented
                t_id = ( flag_is_reverse_complemented, reference_start if flag_is_reverse_complemented else reference_end ) # compose 't_id'
                _add_record( t_id, r ) # add record

        # flush remaining data
        for t_id in list( ns[ 'dict_t_id_to_bucket' ] ) : # retrieve list of 't_id'
            yield _flush_bucket( t_id )


    def _process_batch( p_in, p_out ) :
        """ # 2023-09-15 20:32:17 
        """
        str_uuid_process = bk.UUID( ) # create uuid of the process
        path_file_bam_unsorted = f"{path_folder_temp}{str_uuid_process}.bam"
        path_file_bam_sorted = f"{path_folder_temp}{str_uuid_process}.sorted.bam"
        with pysam.AlignmentFile( path_file_bam_unsorted, 'wb', header = sam_header ) as newsamfile :
            while True :
                batch = p_in.recv( ) # receive a record
                if batch is None :
                    break

                # compose a dataframe containing reads for the bucket
                df = pd.DataFrame( batch, columns = [ 'str_r', 'cb', 'umi', 'length' ] )
                df.sort_values( 'length', inplace = True, ascending = False )
                df.drop_duplicates( subset = [ 'cb', 'umi' ], keep = 'first', inplace = True ) # drop cb-umi duplicates, while keep the longest molecule

                for str_r in df.str_r.values : # for each de-duplicated record
                    r = pysam.AlignedSegment.fromstring( str_r, sam_header ) # compose a pysam record
                    newsamfile.write( r )
                p_out.send( 'completed' )
        # process the output BAM file
        pysam.sort( '-o', path_file_bam_sorted, '-@', str( min( int_num_threads_for_sorting_bam, 5 ) ), path_file_bam_unsorted )
        os.remove( path_file_bam_unsorted ) # remove the temporary file
        pysam.index( path_file_bam_sorted ) # index the given bam file
        p_out.send( 'completed' ) # indicate the work has been completed


    bk.Multiprocessing_Batch_Generator_and_Workers(
        gen_batch = _generate_batch( ),
        process_batch = _process_batch,
        int_num_threads = int_num_processes
        + 2,  # one thread for generating batch, another thread for post-processing of the batch
        flag_wait_for_a_response_from_worker_after_sending_termination_signal = True, # wait until all worker exists before resuming works in the main process
    )

    """ combine results into a single output BAM file """
    path_file_bam_output = f"{path_folder_output}deduplicated_barcoded.bam"
    l_path_file = glob.glob( f"{path_folder_temp}*.sorted.bam" ) # retrieve a list of BAM files to combine
    pysam.merge( '--threads', str( min( int_num_threads_for_sorting_bam, 10 ) ), '-c', '-p', path_file_bam_output, * l_path_file ) # merge output BAM files
    for path_file in l_path_file : # delete the temporary files
        os.remove( path_file )
    pysam.index( path_file_bam_output ) # index the input BAM file
    
    # delete the temporary folder
    shutil.rmtree( path_folder_temp )
    
def FilterInternalPolyAPrimedReadFromBAM( path_file_bam_input : str, path_folder_output : str, int_min_length_internal_polyA_tract : int = 8, name_tag_ia : str = 'IA' ) :
    """ # 2023-10-11 19:56:11 
    Filter internal PolyA primed reads from BAM file
    
    path_file_bam_input : str # an input Barcoded BAM file to filter internal-polyA-region primed reads
    path_folder_output : str # the output folder where splitted Barcoded BAM files will be exported
    int_min_length_internal_polyA_tract : int = 8 # minimum length of an internal poly A/T tract to classify a read as 'internal poly A/T tract'
    name_tag_ia : str = 'IA' # name of the SAM tag containing the length of internal polyA tract. 
    
    """
    # import packages
    import multiprocessing as mp
    import pysam
    import os
    # define functions

    # create the output folder
    os.makedirs( path_folder_output, exist_ok = True )
    
    # read the header of the input BAM file    
    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
        sam_header = samfile.header

    l_name_file = [ 'internal_polyA_primed_reads', 'without_internal_polyA_primed_reads' ] # define the name of the files
    def _write_bam_of_name_file( p_in, p_out ) :
        ''' # 2023-09-07 21:38:10 
        for writing bam file for each 'name_file'
        '''
        name_file = p_in.recv( ) # retrieve the name of file
        path_file_bam = f'{path_folder_output}{name_file}.bam'
        with pysam.AlignmentFile( path_file_bam, 'wb', header = sam_header ) as newsamfile :
            while True :
                batch = p_in.recv( ) # receive a record
                if batch is None :
                    break
                for str_r in batch : # parse the batch
                    r = pysam.AlignedSegment.fromstring( str_r, sam_header ) # compose a pysam record
                    newsamfile.write( r )
        pysam.index( path_file_bam ) # index the given bam file
        p_out.send( 'completed' ) # indicate the work has been completed

    dict_name_file_p_to_writers = dict( )
    l_p_from_writers = [ ]
    l_process = [ ] # list of processes
    for name_file in l_name_file : # for each 'name_file', recruite a worker
        pm2w, pw4m = mp.Pipe( )
        pw2m, pm4w = mp.Pipe( )
        p = mp.Process( target = _write_bam_of_name_file, args = ( pw4m, pw2m ) )
        dict_name_file_p_to_writers[ name_file ] = pm2w 
        l_p_from_writers.append( pm4w )
        p.start( )
        pm2w.send( name_file ) # initialize the worker with 'name_file'
        l_process.append( p ) # collect the process

    # internal setting
    int_max_num_record_in_a_batch = 100
    dict_name_file_to_batch = dict( ( name_file, [ ] ) for name_file in l_name_file ) # initialize 'dict_name_file_to_batch'

    def _flush_batch( name_file : str ) :
        """ # 2023-09-07 22:00:29 
        flush the batch
        """
        batch = dict_name_file_to_batch[ name_file ]
        if len( batch ) > 0 : # if batch is not empty
            dict_name_file_p_to_writers[ name_file ].send( batch ) # send the batch to the writer
            dict_name_file_to_batch[ name_file ] = [ ] # empty the batch    

    def _write_record( name_file : str, str_r : str ) :
        """ # 2023-09-07 22:00:38 
        write a record
        """
        dict_name_file_to_batch[ name_file ].append( str_r ) # add the record
        if len( dict_name_file_to_batch[ name_file ] ) >= int_max_num_record_in_a_batch :
            _flush_batch( name_file ) # flush the batch

    # read file and write the record
    with pysam.AlignmentFile( path_file_bam_input, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            dict_tags = dict( r.get_tags( ) ) # retrieve tags of the read
            if name_tag_ia in dict_tags :
                str_r = r.tostring( ) # convert samtools record to a string
                _write_record( 'internal_polyA_primed_reads' if dict_tags[ name_tag_ia ] >= int_min_length_internal_polyA_tract else 'without_internal_polyA_primed_reads', str_r ) # write the record
                
    # flush remaining records
    for name_file in l_name_file :
        _flush_batch( name_file )

    # notify all works in the main process has been completed
    for name_file in l_name_file :
        dict_name_file_p_to_writers[ name_file ].send( None )

    # wait for all workers to complete their jobs
    for p in l_p_from_writers :
        p.recv( ) # receive a signal indicating the worker has dismissed itself
    # pipeline completed
    return
    
if __name__ == "__main__":
    ourotools()  # run ouro at the top-level environment
