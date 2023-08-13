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
from io import StringIO
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
_version_ = "0.0.1"
_scelephant_version_ = _version_
_last_modified_time_ = "2023-08-10 16:40:52"

str_release_note = [
    """
    # %% RELEASE NOTE %%
    # 2023-08-10 16:41:24 
    draft version of 'LongFilterNSplit' function completed

    ##### Future implementations #####

    """
]

str_description = """
 ______     __  __     ______     ______     ______   ______     ______     __         ______    
/\  __ \   /\ \/\ \   /\  == \   /\  __ \   /\__  _\ /\  __ \   /\  __ \   /\ \       /\  ___\   
\ \ \/\ \  \ \ \_\ \  \ \  __<   \ \ \/\ \  \/_/\ \/ \ \ \/\ \  \ \ \/\ \  \ \ \____  \ \___  \  
 \ \_____\  \ \_____\  \ \_\ \_\  \ \_____\    \ \_\  \ \_____\  \ \_____\  \ \_____\  \/\_____\ 
  \/_____/   \/_____/   \/_/ /_/   \/_____/     \/_/   \/_____/   \/_____/   \/_____/   \/_____/ 
                                                                                                 
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
):
    """# 2023-01-06 23:46:20
    convert df_count (ouro output) to 10X MTX (matrix market) format in a memory-efficient manner.

    path_file_df_count : str, # file path to 'df_count'
    path_folder_mtx_10x_output : str, # a folder containing 10x output matrix (unfiltered)
    path_folder_mtx_10x_filtered_output : str, # a folder containing 10x output matrix (filtered)
    chunksize : int = 500000,
    int_min_count_features_for_filtering_barcodes : int = 50, # the minimum number of features in a barcode to be included in the filtered output
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
    with gzip.open(f"{path_folder_mtx_10x_output}matrix.mtx.gz", "wb") as newfile:
        newfile.write(
            f"""%%MatrixMarket matrix coordinate integer general\n%\n{int_num_features} {int_num_barcodes} {int_num_records}\n""".encode()
        )
        with gzip.open(
            f"{path_folder_mtx_10x_filtered_output}matrix.mtx.gz", "wb"
        ) as newfile_filtered:
            newfile_filtered.write(
                f"""%%MatrixMarket matrix coordinate integer general\n%\n{int_num_features} {int_num_barcodes_filtered} {int_num_records_filtered}\n""".encode()
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

def LongFilterNSplit(
    flag_usage_from_command_line_interface: bool = False,
    path_file_minimap_index_genome: Union[str, None] = None,
    l_path_file_minimap_index_unwanted: List[str] = [ ],
    l_path_file_fastq_input: Union[list, None] = None,
    l_path_folder_output: [list[str], None] = None,
    n_threads: int = 32,
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_reads_in_a_batch : int = 10_000, # the number of reads in a batch
    float_memory_in_GiB: float = 50,
    verbose: bool = True,
    str_minimap_aligner_preset : str = 'splice', # preset of the minimap2 aligner
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome
    int_size_window_for_searching_poly_a_tail : int = 15, # the size of the window from the end of the alignment to search for poly A tail.
    int_max_size_intervening_sequence_between_alignment_end_and_poly_A : int = 20, # max size of the intervening sequence between alignment end position and poly A tract. it will be applied for both internal poly A or external (enzymatically attached) poly A.
    float_min_A_frequency_for_identifying_poly_A : float = 0.75, # the minimum frequency to determine a sequence contains a poly A tract
    int_min_size_intervening_sequence_for_splitting : int = 150, # the minimum length of intervening sequence between alignments for splitting the reads
    int_max_intron_size_for_determining_chimeric_molecule : int = 200000, # the maximum allowed intron size for classifying considering the molecule as a intra-chromosomal chimeric read
    am_genome = None, # mappy aligner for genome (optional. if given, will override 'path_file_minimap_index_genome' argument)
    l_am_unwanted : Union[ None, List ] = None, # mappy aligner for unwanted sequences (optional. if given, will override 'l_path_file_minimap_index_unwanted' argument)
) -> None :
    """# 2023-08-10 16:30:50 
    
    flag_usage_from_command_line_interface: bool = False,
    path_file_minimap_index_genome: Union[str, None] = None, # required for identifying valid regions of a read and identify chimeric transcripts
    l_path_file_minimap_index_unwanted: List[str] = [ ], # a list of minimap indices of sequences to which unwanted reads can originate (e.g., mt-dna, rRNA repeat, etc.) in a decreasing order of priority
    l_path_file_fastq_input: Union[list, None] = None,
    l_path_folder_output: [list[str], None] = None,
    n_threads: int = 32,
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_reads_in_a_batch : int = 10_000, # the number of reads in a batch
    float_memory_in_GiB: float = 50,
    str_minimap_aligner_preset = 'splice', # preset of the minimap2 aligner
    verbose: bool = True,
    int_min_mapq = 1, # minimum mapping quality of the alignment to consider a read (or parts of a read)  were aligned to the genome
    int_size_window_for_searching_poly_a_tail : int = 16, # the size of the window from the end of the alignment to search for poly A tail.
    int_max_size_intervening_sequence_between_alignment_end_and_poly_A : int = 20, # max size of the intervening sequence between alignment end position and poly A tract. it will be applied for both internal poly A or external (enzymatically attached) poly A.
    float_min_A_frequency_for_identifying_poly_A : float = 0.75, # the minimum frequency to determine a sequence contains a poly A tract
    int_max_intron_size_for_determining_chimeric_molecule : int = 200000, # the maximum allowed intron size for classifying considering the molecule as a intra-chromosomal chimeric read
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
            "--int_num_reads_in_a_batch",
            help="(default: 10000) the number of reads in a batch.",
            default=10_000,
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
        
        args = parser.parse_args()

        flag_usage_from_command_line_interface = args.flag_usage_from_command_line_interface
        path_file_minimap_index_genome = args.path_file_minimap_index_genome
        l_path_file_fastq_input = args.l_path_file_fastq_input
        l_path_folder_output = args.l_path_folder_output
        n_threads = args.n_threads
        int_num_samples_analyzed_concurrently = args.int_num_samples_analyzed_concurrently
        float_memory_in_GiB = args.float_memory_in_GiB
        verbose = args.verbose
        int_num_reads_in_a_batch = args.int_num_reads_in_a_batch
        l_path_file_minimap_index_unwanted = args.l_path_file_minimap_index_unwanted
        int_min_mapq = args.int_min_mapq
        str_minimap_aligner_preset = args.str_minimap_aligner_preset
        int_size_window_for_searching_poly_a_tail = args.int_size_window_for_searching_poly_a_tail
        int_max_size_intervening_sequence_between_alignment_end_and_poly_A = args.int_max_size_intervening_sequence_between_alignment_end_and_poly_A
        float_min_A_frequency_for_identifying_poly_A = args.float_min_A_frequency_for_identifying_poly_A
        int_min_size_intervening_sequence_for_splitting = args.int_min_size_intervening_sequence_for_splitting
        int_max_intron_size_for_determining_chimeric_molecule = args.int_max_intron_size_for_determining_chimeric_molecule

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
    Pipeline specific functions
    """
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
    
    def _classify_read( seq : str, hit ) :
        """ # 2023-08-03 23:15:49 
        classlfy read by detecting poly A and unreferenced Gs.
        
        hit # a Minimap2 mappy alignment record. multiple hits can be given as a list of hits. 
        
        # possible labels
        'no_poly_A', 'external_poly_A__unrefenced_G', 'internal_poly_A__unrefenced_G', 'external_poly_A__no_unrefenced_G', 'internal_poly_A__no_unrefenced_G'
        
        return the direction of the read and the classification label
        """
        # retrieve 'q_st' and 'q_en'
        if isinstance( hit, list ) : # if a list of hits were given, 
            q_st, q_en = min( h.q_st for h in hit ), max( h.q_en for h in hit ) # retrieve smallest q_st and largest q_en for q_st and q_en of the segment if multiple hits were given
        else : # when single hit was given
            q_st, q_en = hit.q_st, hit.q_en
        
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
                    return 'internal_poly_A__unrefenced_G', '-'
                else :
                    return 'internal_poly_A__no_unrefenced_G', '-'
            else :
                if flag_unreferenced_Gs :
                    return 'internal_poly_A__unrefenced_G', '+'
                else :
                    return 'internal_poly_A__no_unrefenced_G', '+'
        
        ''' handles external poly A '''
        if flag_external_poly_A_flipped :
            if flag_unreferenced_Gs_flipped :
                return 'external_poly_A__unrefenced_G', '-'
            else :
                return 'external_poly_A__no_unrefenced_G', '-'
        else :
            if flag_unreferenced_Gs :
                return 'external_poly_A__unrefenced_G', '+'
            else :
                return 'external_poly_A__no_unrefenced_G', '+'
    
    def _initialize_dict_arr_dist( ) :
        """ # 2023-08-03 11:49:26 
        initialize 'dict_arr_dist'
        """
        return {
            'aligned_to_unwanted_sequence' : None,
            'cannot_aligned_to_genome' : None,
            'aligned_to_genome' : None,

            'aligned_to_genome__non_chimeric__no_poly_A' : None,
            'aligned_to_genome__non_chimeric__external_poly_A__unrefenced_G' : None,
            'aligned_to_genome__non_chimeric__internal_poly_A__unrefenced_G' : None,
            'aligned_to_genome__non_chimeric__external_poly_A__no_unrefenced_G' : None,
            'aligned_to_genome__non_chimeric__internal_poly_A__no_unrefenced_G' : None,

            'aligned_to_genome__chimeric__no_poly_A' : None,
            'aligned_to_genome__chimeric__external_poly_A__unrefenced_G' : None,
            'aligned_to_genome__chimeric__internal_poly_A__unrefenced_G' : None,
            'aligned_to_genome__chimeric__external_poly_A__no_unrefenced_G' : None,
            'aligned_to_genome__chimeric__internal_poly_A__no_unrefenced_G' : None,
        }
        
    def run_pipeline():
        """# 2023-07-30 17:20:19 
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
                "int_num_reads_in_a_batch" : int_num_reads_in_a_batch,
                'str_minimap_aligner_preset' : str_minimap_aligner_preset,
                'int_min_mapq' : int_min_mapq,
                'int_size_window_for_searching_poly_a_tail' : int_size_window_for_searching_poly_a_tail,
                'float_min_A_frequency_for_identifying_poly_A' : float_min_A_frequency_for_identifying_poly_A,
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
                int_read_counter = 0 # count read
                l_r_for_a_batch = [ ] # a list of records for a batch
                for r in bk.FASTQ_Iterate( path_file_fastq_input ) : # iterate through the input FASTQ file
                    l_r_for_a_batch.append( r ) # add the record
                    int_read_counter += 1 # increase the counter
                    if int_read_counter % int_num_reads_in_a_batch == 0 : # if the batch is full, yield the batch
                        yield l_r_for_a_batch
                        l_r_for_a_batch = [ ] # initialize the next batch
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
                dict_newfile_fastq_output = {
                    'aligned_to_unwanted_sequences' : gzip.open( f"{path_folder_temp}{str_uuid}.aligned_to_unwanted_sequences.fastq.gz", "wb", ), 
                    'cannot_aligned_to_genome' : gzip.open( f"{path_folder_temp}{str_uuid}.cannot_aligned_to_genome.fastq.gz", "wb", ), 

                    'aligned_to_genome__non_chimeric__no_poly_A' : gzip.open( f"{path_folder_temp}{str_uuid}.aligned_to_genome__non_chimeric__no_poly_A.fastq.gz", "wb", ), 
                    'aligned_to_genome__non_chimeric__poly_A__plus_strand' : gzip.open( f"{path_folder_temp}{str_uuid}.aligned_to_genome__non_chimeric__poly_A__plus_strand.fastq.gz", "wb", ), # main fastq output file

                    'aligned_to_genome__chimeric__no_poly_A' : gzip.open( f"{path_folder_temp}{str_uuid}.aligned_to_genome__chimeric__no_poly_A.fastq.gz", "wb", ), 
                    'aligned_to_genome__chimeric__poly_A__plus_strand' : gzip.open( f"{path_folder_temp}{str_uuid}.aligned_to_genome__chimeric__poly_A__plus_strand.fastq.gz", "wb", ), 
                }
                
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
                    def _process_molecule( qname : str, seq : str, qual : str, hit, st : Union[ None, int ] = None, en : Union[ None, int ] = None, type_molecule : Literal[ 'non_chimeric', 'chimeric' ] = 'non_chimeric' ) :
                        """ # 2023-08-03 10:28:48 
                        process a molecule (segment) of a sequence
                        
                        hit, # a Minimap2 mappy alignment record. multiple hits can be given as a list of hits. 
                        type_molecule : Literal[ 'non_chimeric', 'chimeric' ] = 'non_chimeric' # type of the molecule
                        """
                        ''' Retrieve the segment and the length of the segment. Also, compose suffix to the qname '''
                        if st is None and en is None : # analyze an entire molecule
                            seg, qual_seg, qname_suffix = seq, qual, '' 
                        elif st is None :
                            seg, qual_seg, qname_suffix = seq[ : en ], qual[ : en ], '_to' + str( en )
                        elif en is None :
                            seg, qual_seg, qname_suffix = seq[ st : ], qual[ st : ], '_' + str( st ) + 'to'
                        else :
                            seg, qual_seg, qname_suffix = seq[ st : en ], qual[ st : en ], '_' + str( st ) + 'to' + str( en )
                        len_seg = len( seq )
                
                        label, direction = _classify_read( seq, hit ) # classify the segment 

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
                        flag_is_segment_chimeric, q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev = False, None, None, None, None, None # initialize a flag indicating the segment is chimeric or non-chimeric
                        l_hit_segment = [ ] # initialize a list of hit of a segment
                        for q_st, hit in arr_algn : # for each q_st and hit
                            q_en = hit.q_en # retrieve q_en
                            ''' initialize the segment (if it was not initialized) '''
                            if q_st_segment is None :
                                q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev = 0, hit.q_en, hit.ctg, hit.r_st, hit.r_en # use the start of the molecule as 'q_en_segment'
                                l_hit_segment.append( hit )
                                continue # continue to the next hit
                            
                            """ split the read """
                            int_size_gap = q_st - q_en_segment # calculate the size of the intervening sequence
                            if int_size_gap >= int_min_size_intervening_sequence_for_splitting : # if the size of the intervening sequences is larger then the threshold, split the read
                                int_size_flanking = min( int_size_gap, int_min_size_intervening_sequence_for_splitting ) # retrieve the size of the flanking sequence to include in the segment. the size of the flanking sequence cannot be larger than the intervening sequence
                                _process_molecule( qname, seq, qual, l_hit_segment, q_st_segment, q_en_segment + int_size_flanking, 'chimeric' if flag_is_segment_chimeric else 'non_chimeric' ) # process chimeric segment # add 'int_size_flanking' to q_en_segment to include a flanking sequence
                                # initialize the next segment
                                flag_is_segment_chimeric, q_st_segment, q_en_segment, ctg_prev, r_st_prev, r_en_prev = False, q_st - int_size_flanking, hit.q_en, hit.ctg, hit.r_st, hit.r_en # set q_st_segment as q_st - int_size_flanking to include a flanking sequence
                                l_hit_segment = [ hit ] 
                                continue
                                
                            """ concatenate genomic alignments and determine whether the segment is chimeric or not """
                            # identify chimeric molecule
                            if ctg_prev != hit.ctg :
                                flag_is_segment_chimeric = True
                            elif max( r_st_prev - hit.r_en, hit.r_st - r_en_prev ) > int_max_intron_size_for_determining_chimeric_molecule : # if the distance between alignment is longer than maximum intron size, consider reads as an intra-chromosomal chimeric molecule
                                flag_is_segment_chimeric = True
                            # extend segment
                            l_hit_segment.append( hit )
                            q_en_segment = max( q_en_segment, q_en ) # update 'q_en_segment' ('q_st_segment' will not change)
                            ctg_prev, r_st_prev, r_en_prev = hit.ctg, hit.r_st, hit.r_en # update the previous alignment
                                
                        if len( l_hit_segment ) > 0 : # if a segment is remaining, process the segment
                            _process_molecule( qname, seq, qual, l_hit_segment, q_st_segment, None, 'chimeric' if flag_is_segment_chimeric else 'non_chimeric' ) # process chimeric segment # use the end of the molecule as 'q_en_segment'
                    
                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : int_total_num_records_for_a_batch,
                        'dict_arr_dist' : dict_arr_dist,
                    } )  # report the number of processed records
                    
                """ close output files """
                for name_type in dict_newfile_fastq_output :
                    dict_newfile_fastq_output[ name_type ].close( )
                    
                """ report the worker has completed all works """
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")
                pipe_sender.send( 'completed' )  

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
                for name_file in [ 
                    "aligned_to_unwanted_sequences.fastq.gz",
                    "cannot_aligned_to_genome.fastq.gz",
                    "aligned_to_genome__non_chimeric__no_poly_A.fastq.gz",
                    "aligned_to_genome__non_chimeric__poly_A__plus_strand.fastq.gz",
                    "aligned_to_genome__chimeric__no_poly_A.fastq.gz",
                    "aligned_to_genome__chimeric__poly_A__plus_strand.fastq.gz",
                ] :
                    bk.OS_Run( [ 'cat' ] + glob.glob( f"{path_folder_temp}*.{name_file}" ), stdout_binary = True, path_file_stdout = f"{path_folder_output}{name_file}" )
                    

                # draw plots 
                # plot settings
                int_max_molecule_size_plot = 10000
                for name_cat_dist in _initialize_dict_arr_dist( ) : # for each category
                    if ns[ 'dict_arr_dist' ][ name_cat_dist ] is not None :
                        len_max_molecule_size_data = len( ns[ 'dict_arr_dist' ][ name_cat_dist ] ) # retrieve max molecule size 
                        plt.plot( np.arange( min( int_max_molecule_size_plot, len_max_molecule_size_data ) ), ns[ 'dict_arr_dist' ][ name_cat_dist ] if len_max_molecule_size_data <= int_max_molecule_size_plot else ns[ 'dict_arr_dist' ][ name_cat_dist ][ : int_max_molecule_size_plot ] )
                        plt.title( f"{name_cat_dist} ({ns[ 'dict_arr_dist' ][ name_cat_dist ].sum( )} molecules)" )
                        bk.MPL_SAVE( f"{name_cat_dist}.distribution", folder = path_folder_graph, l_format=['.pdf', '.png'] )
                    
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
    int_num_reads_in_a_batch : int = 10_000, # the number of reads in a batch
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to filter read with low alignment quality
    float_memory_in_GiB : float = 50, # expected memory usage of the pipeline
    float_error_rate : float = 0.2, # maximum error rate to consider when searching adaptor sequence in the read
    int_length_cb : int = 16, # the length of the cell barcode
    int_length_umi : int = 12, # the length of the UMI (unique molecular identifier)
    str_seq_r1 : str = 'CTACACGACGCTCTTCCGATCT', # the sequence of R1 adaptor (in 10x GEX v3 kit, located upstream of CB and UMI)
    str_seq_tso : str = 'AAGCAGTGGTATCAACGCAGAG', # the sequence of TSO adaptor (in 10x GEX v3 kit, located at 5' end of the molecule)
    path_file_valid_cb : str = None, # (required argument) the path to tsv file of whitelist barcodes. For more details, please see 10x cellranger references.
    n_cell_expected : int = 1000, # the number of expected cells
    int_len_sliding_window_internal_polyT : int = 10, # the length of sliding window for searching internal poly T (poly A) tract. (When poly-A tailed read is reverse complemented, R1 adaptor become situated in the forward direction
    int_len_window_internal_polyT : int = 30, # the size of window for searching for internal poly T
    float_min_T_fraction : float = 0.8, # the minimum T fraction for identifying the stretch of poly T tract
    int_min_n_overlap_kmer_for_clustering_umi : int = 1, # the minimum number of overlapped kmer for initiating UMI clustering 
    int_len_kmer_for_clustering_umi : int = 7, # the length of kmer for clustering UMI
    float_min_proportion_read_to_select_kmer_representative_for_clustering_umi : float = 0.75, # if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.
    verbose: bool = True,
) -> None :
    """# 2023-07-30 16:10:37 
    
    flag_usage_from_command_line_interface: bool = False, # a flag indicating the usage in the command line
    l_path_file_bam_input: Union[list, None] = None, # list of input BAM files
    l_path_folder_output: [list[str], None] = None, # list of output folders
    n_threads: int = 32, # the number of threads to use
    int_num_samples_analyzed_concurrently : int = 2, # the number of samples that can be analyzed concurrently to reduce bottlenecks due to processing of very large chunks.
    int_num_reads_in_a_batch : int = 10_000, # the number of reads in a batch
    int_min_mapq : int = 1, # minimum mapping quality of the alignment to filter read with low alignment quality
    float_memory_in_GiB: float = 50,
    float_error_rate : float = 0.2, # maximum error rate to consider when searching adaptor sequence in the read
    int_length_cb : int = 16, # the length of the cell barcode
    int_length_umi : int = 12, # the length of the UMI (unique molecular identifier)
    str_seq_r1 : str = 'CTACACGACGCTCTTCCGATCT', # the sequence of R1 adaptor (in 10x GEX v3 kit, located upstream of CB and UMI)
    str_seq_tso : str = 'AAGCAGTGGTATCAACGCAGAG', # the sequence of TSO adaptor (in 10x GEX v3 kit, located at 5' end of the molecule)
    path_file_valid_cb : str = None, # (required argument) the path to tsv file of whitelist barcodes. For more details, please see 10x cellranger references.
    n_cell_expected : int = 1000, # the number of expected cells
    int_len_sliding_window_internal_polyT : int = 10, # the length of sliding window for searching internal poly T (poly A) tract. (When poly-A tailed read is reverse complemented, R1 adaptor become situated in the forward direction
    int_len_window_internal_polyT : int = 30, # the size of window for searching for internal poly T
    float_min_T_fraction : float = 0.8, # the minimum T fraction for identifying the stretch of poly T tract
    int_min_n_overlap_kmer_for_clustering_umi : int = 1, # the minimum number of overlapped kmer for initiating UMI clustering 
    int_len_kmer_for_clustering_umi : int = 7, # the length of kmer for clustering UMI
    float_min_proportion_read_to_select_kmer_representative_for_clustering_umi : float = 0.75, # if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.
    verbose: bool = True,
    
    returns None
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
            "--int_num_reads_in_a_batch",
            help="(default: 10000) the number of reads in a batch.",
            default=10_000,
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
        arg_grp_cb_correction.add_argument( "-N", "--n_cell_expected", help = "(default: 1000) the number of expected cells", default = 1000, type = int )

        arg_grp_internal_polyt = parser.add_argument_group("Internal Poly(A) Tract-Primed Read Identification")
        arg_grp_internal_polyt.add_argument( "-S", "--int_len_sliding_window_internal_polyT", help = "(default: 10) the length of sliding window for searching internal poly T (poly A) tract. (When poly-A tailed read is reverse complemented, R1 adaptor become situated in the forward direction", type = int, default = 10 )
        arg_grp_internal_polyt.add_argument( "-w", "--int_len_window_internal_polyT", help = "(default: 30) the size of window for searching for internal poly T", type = int, default = 30 )
        arg_grp_internal_polyt.add_argument( "-F", "--float_min_T_fraction", help = "(default: 0.8) the minimum T fraction for identifying the stretch of poly T tract", type = float, default = 0.8 )
                    
        arg_grp_umi_clustering = parser.add_argument_group("UMI Clustering")
        arg_grp_umi_clustering.add_argument( "-O", "--int_min_n_overlap_kmer_for_clustering_umi", help = "(default: 1) the minimum number of overlapped kmer for initiating UMI clustering ", type = int, default = 1 )
        arg_grp_umi_clustering.add_argument( "-L", "--int_len_kmer_for_clustering_umi", help = "(default: 7) the length of kmer for clustering UMI", type = int, default = 7 )
        arg_grp_umi_clustering.add_argument( "-P", "--float_min_proportion_read_to_select_kmer_representative_for_clustering_umi", help = "(default: 0.75) if the given proportion of UMI contains the kmer, include the kmer in a set of kmers representing the UMI clusters.", type = float, default = 0.75 )
        
        args = parser.parse_args( )

        l_path_file_bam_input = args.l_path_file_bam_input
        l_path_folder_output = args.l_path_folder_output
        n_threads = args.n_threads
        int_num_reads_in_a_batch = args.int_num_reads_in_a_batch
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
        n_cell_expected = args.n_cell_expected
        int_len_sliding_window_internal_polyT = args.int_len_sliding_window_internal_polyT
        int_len_window_internal_polyT = args.int_len_window_internal_polyT
        float_min_T_fraction = args.float_min_T_fraction
        int_min_n_overlap_kmer_for_clustering_umi = args.int_min_n_overlap_kmer_for_clustering_umi
        int_len_kmer_for_clustering_umi = args.int_len_kmer_for_clustering_umi
        float_min_proportion_read_to_select_kmer_representative_for_clustering_umi = args.float_min_proportion_read_to_select_kmer_representative_for_clustering_umi

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
    def _check_binary_flags( flags : int, int_bit_flag_position : int ) :
        """ # 2023-08-08 22:47:02 
        check a flag in the binary flags at the given position
        """
        return ( flags & ( 1 << int_bit_flag_position ) ) > 0 

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

    
    def run_pipeline():
        """# 2023-07-30 17:20:19 
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
                "int_num_reads_in_a_batch" : int_num_reads_in_a_batch,
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
                    try :
                        r = next( gen_r ) # retrieve the first read
                    except StopIteration : # if the bam file is emtpy, end the generator
                        return
                    # initialize the batch
                    ns_batch = { 'int_num_reads_encountered_for_a_batch' : 1, 'start__reference_name' : r.reference_name, 'start__reference_start' : r.reference_start, } # initialize the dictionary containing information about the batch # counts of reads in a batch
                    
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
                        cigartuples, flags = r.cigartuples, r.flag # retrieve attributes
                        if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                            continue 
                        if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                            continue
                            
                        """ when contig changes """
                        if r.reference_name != ns_batch[ 'start__reference_name' ] :
                            yield ns_batch # yield the last batch for the last contig
                            # initialize the next batch
                            ns_batch = { 'int_num_reads_encountered_for_a_batch' : 0 } # initialize the counter
                            ns_batch[ 'start__reference_name' ] = r.reference_name
                            ns_batch[ 'start__reference_start' ] = r.reference_start
                            
                        ns_batch[ 'int_num_reads_encountered_for_a_batch' ] += 1 # increase the read count
                        if int_num_reads_in_a_batch <= ns_batch[ 'int_num_reads_encountered_for_a_batch' ] : # once the batch is full, yield the batch and consume remaining reads starting at the reference start position, so that the reads of the same reference start position are processed together.
                            # update batch information
                            ns_batch[ 'end__reference_start' ] = r.reference_start
                            while True :
                                """ retrieve a read """ 
                                try :
                                    r = next( gen_r )
                                except StopIteration : # once all reads were analyzed, exit the loop
                                    break
                                    
                                if ns_batch[ 'end__reference_start' ] != r.reference_start : # when the 'reference_start' position changes, finish the batch
                                    break
                                
                                ns_batch[ 'int_num_reads_encountered_for_a_batch' ] += 1 # increase the counter
                            yield ns_batch # yield the batch
                            ns_batch = { 'int_num_reads_encountered_for_a_batch' : 1 } # initialize the counter
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
                        for r in samfile.fetch( start__reference_name, start__reference_start, end__reference_start + 1 ) if end__reference_start is not None else samfile.fetch( start__reference_name, start__reference_start ) : # include the end position by adding +1
                            
                            """ filter read """
                            if r.mapq < int_min_mapq : # filter out reads with low mapq
                                continue
                            if r.seq is None : # consider only the primary alignment
                                continue
                            seq, cigartuples, flags = r.seq, r.cigartuples, r.flag # retrieve attributes
                            if int_cigarop_H == cigartuples[ 0 ][ 0 ] or int_cigarop_H == cigartuples[ -1 ][ 0 ] : # skip hard-clipped reads
                                continue 
                            if _check_binary_flags( flags, 10 ) or _check_binary_flags( flags, 8 ) : # filter out optical duplicates or secondary alignments
                                continue
                                
                            ''' if the batch has been completed, exit the loop '''
                            if end__reference_start is not None and r.reference_end > end__reference_start : 
                                break
                            
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
                            l_tags = [ ( 'RX', res_r1[ 'num_errors' ], 'i' ), ( 'AX', res_tso[ 'num_errors' ], 'i' ) ] # add the number of errors from R1 and TSO adaptor search results as tags
                            
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
                                l_tags += [ ("CX", seq_cb_umi, 'Z'), ('PX', int_length_internal_polyT, 'i') ] # add uncorrected cb and umi sequence as a tag # add the identified poly T length as a tag
                                # collect data 
                                l_cb_umi.append( seq_cb_umi ) # collect 'seq_cb_umi'

                            ''' write the SAM record ''' 
                            int_total_num_records_processed += 1
                            r.set_tags( l_tags ) # set tags
                            newsamfile.write( r ) # write the record to the output BAM file
                            
                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : int_total_num_records_processed, # record the actual number of records processed for the batch
                        'l_cb_umi' : l_cb_umi,
                    } )  # report the number of processed records

                """ close output files """
                newsamfile.close( )
                # index the resulting BAM file
                pysam.index( path_file_bam_preprocessed )
                
                """ report the worker has completed all works """
                pipe_sender.send( 'completed' )  
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")

            ns = { 'int_num_read_currently_processed' : 0, 'int_num_records_with_cb_umi' : 0, 'l_cb_umi' : [ ] }  # define a namespace # initialize total number of reads processed by the algorithm

            def post_process_batch(res):
                # update data using the received result
                ns["int_num_read_currently_processed"] += res[ 'int_total_num_records_for_a_batch' ]
                ns["int_num_records_with_cb_umi"] += len( res["l_cb_umi"] ) # update ns["int_num_records_with_cb_umi"]
                ns["l_cb_umi"] += res["l_cb_umi"]
                logger.info(
                    f"[{path_file_bam_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed. CB/UMI sequence identification rate is {np.round(ns['int_num_records_with_cb_umi'] / ns['int_num_read_currently_processed'], 2 ) if ns['int_num_read_currently_processed'] > 0 else np.nan}"
                )  # report
                    
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
            
            bk.PICKLE_Write( f"{path_folder_output}l_cb_umi.pickle", ns["l_cb_umi"] )# ❤️
            
            """ combine results into a single output BAM file """
            path_file_bam_preprocessed = f"{path_folder_temp}preprocessed.bam"
            pysam.merge( '--threads', str( min( n_threads, 10 ) ), '-c', '-p', path_file_bam_preprocessed, * glob.glob( f"{path_folder_temp}*.preprocessed.bam" ) ) # merge output BAM files
            pysam.index( path_file_bam_preprocessed ) # index the input BAM file

            """ 
            Correct and Assign Cell Barcodes to Each Read
            """
            # internal settings
            n_droplets = 100000 # number of droplets generated in 10X Chromium instruments (number of barcoded beads) - with sufficient margin
            n_minimum_count_cb_before_correction = 2 # threshold for filtering cell barcodes before correction
            ''' retrieve list of potentially valid barcodes '''
            l_cb_umi = ns["l_cb_umi"]
            # retrieve all cb sequence by its normal length, and count each unique cb
            s_cb_count = pd.Series( bk.COUNTER( list( e[ : int_length_cb ] for e in l_cb_umi ) ) )
            # filtering uncorrected cb
            s_cb_count = s_cb_count[ s_cb_count >= n_minimum_count_cb_before_correction  ]
            # sort uncorrected cb by their counts and only retrieve 'n_droplets' number of cell barcodes
            if len( s_cb_count ) > n_droplets :
                s_cb_count = s_cb_count.sort_values( ascending = False ).iloc[ : n_droplets ]
            # drop invalid barcodes
            s_cb_count = bk.Series_Subset( s_cb_count, set_valid_cb )
            # retrieve a list of valid barcodes
            l_cb_valid = s_cb_count.index.values

            ''' 1) use pre-computed error-correcting dictionary to identify cell barcodes with 1 error '''
            # build a list of possible errors for each base
            dict_base_to_l_error = dict( )
            for str_base in 'ATGC' :
                l = [ '' ]
                for str_base_error in 'ATGC' :
                    l.append( str_base + str_base_error )
                    if str_base != str_base_error :
                        l.append( str_base_error )
                dict_base_to_l_error[ str_base ] = l
            # retrieve mapping between cb with error to error-free cb
            dict_cb_with_error_to_cb = dict( ) 
            for cb in l_cb_valid :
                for pos in range( int_length_cb ) :
                    str_base = cb[ pos ]
                    for error in dict_base_to_l_error[ str_base ] :
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

            ''' define a function for correcting CB sequences retrieved from reads using the different levels of dictionaries '''
            def _correct_cell_barcode( cb_umi_padded : str ) :
                """ # 2023-08-09 22:03:25 
                correct a single cell barcode 
                """
                cb_corrected = np.nan # initialized to 'cb not found'
                ''' 1) use pre-computed error-correcting dictionary to identify cell barcodes with 1 error '''
                for cb_from_read in [ cb_umi_padded[ : int_length_cb ], cb_umi_padded[ : int_length_cb - 1 ], cb_umi_padded[ : int_length_cb + 1 ] ] :
                    if cb_from_read in dict_cb_with_error_to_cb :
                        cb_corrected = dict_cb_with_error_to_cb[ cb_from_read ]
                        break
                if not isinstance( cb_corrected, float ) : # if the correct cell barcode was assigned, return the result
                    return cb_corrected
                ''' 2) Using varying number of kmer to identify cell barcodes with many number of errors '''
                cb_from_read = cb_umi_padded[ : int_length_cb ] # only consider 'int_length_cb' number of bases
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
                """ # 2023-08-12 18:44:33 
                cluster UMI sequences of a single cell barcode
                """
                ''' handle simple cases '''
                if len( l_umi_for_clustering ) <= 1 : # does not perform UMI clustering when the number of given UMI sequences is <= 1
                    return l_umi_for_clustering
                
                ''' perform UMI clustering '''
                # cluster umi with extracted k-mer
                dict_cluster = SEQ.Cluster_with_Kmer( dict_seq_count = bk.COUNTER( l_umi_for_clustering ), int_min_n_overlap_kmer = int_min_n_overlap_kmer_for_clustering_umi, len_kmer = int_len_kmer_for_clustering_umi, float_min_proportion_read_to_select_kmer_representative = float_min_proportion_read_to_select_kmer_representative_for_clustering_umi )

                ''' retrieve a representative UMI sequence for each cluster based on the frequency, and replace the UMI sequences of the cluster with the representative sequence of the cluster '''
                dict_seq_umi_to_str_name_umi_cluster = dict( ) # retrieve mapping
                for id_c in dict_cluster : # for each cluster
                    c = dict_cluster[ id_c ] # retrieve the cluster
                    # set the seq_umi with the largest count as a name of current umi cluster
                    str_name_umi_cluster, _ = bk.DICTIONARY_Find_Max( c[ 'seq_count' ] )
                    # assign the current umi cluster name to the umi sequences belonging to the current cluster
                    for seq_umi in c[ 'seq_count' ] :
                        dict_seq_umi_to_str_name_umi_cluster[ seq_umi ] = str_name_umi_cluster

                l_umi_corrected = list( dict_seq_umi_to_str_name_umi_cluster[ umi_uncorrected ] for umi_uncorrected in l_umi_for_clustering )
                return l_umi_corrected
            
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
                workers_for_bucket_processing = bk.Offload_Works( int_num_workers_for_bucket_processing )  # 💛no limit for the number of works that can be submitted.
                int_min_num_reads_in_a_bucket_for_parallel_processing = 100 # minimum number of reads in a bucket for offloaded in a worker processs
                    
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
                    ns = { 'int_total_num_records_processed' : 0, 'int_bucket_deletion_count' : 0 } # create a namespace for buckets # a counter counting the number of bucket deleted from 'dict_poly_a_site_to_l_l'. if the number exceed
                    ns[ 'dict_poly_a_site_to_l_l' ] = dict( ) # a container to collect reads for alignment end position
                    reference_name_current = None # initialize the current contig name
                    
                    set_e = set( ) # ❤️
                    
                    def _process_bucket( l_seq_cb_umi ) :
                        """ # 2023-08-12 15:34:56 
                        process reads in a bucket, and return the reads
                        """
                        """ Correct the cell barcode sequence, and collect the read for UMI clustering """
                        l_cb = [ ]
                        for seq_cb_umi in l_seq_cb_umi : # retrieve cb_umi sequence
                            cb_corrected = _correct_cell_barcode( seq_cb_umi ) # correct cell barcode
                            l_cb.append( cb_corrected )
                            
                        ''' retrieve UMI sequence (uncertain boundaries, might contain more or less bases than the actual sequenced UMI sequence) '''
                        l_umi_uncorrected = [ ]
                        for seq_cb_umi, cb_assigned in zip( l_seq_cb_umi, l_cb ) :
                            l_umi_uncorrected.append( np.nan if isinstance( cb_assigned, float ) else seq_cb_umi[ int_length_cb : ] )
                        
                        
                        if len( l_umi_uncorrected ) <= 1 : # if there is less than 2 umi sequences, clustering is not required.
                            ''' handle simple cases '''
                            l_umi_corrected = l_umi_uncorrected
                        else :
                            l_umi_uncorrected = np.array( l_umi_uncorrected, dtype = object ) # convert to numpy array
                            l_umi_corrected = np.zeros( len( l_umi_uncorrected ), dtype = object ) # initialize an empty array
                            # l_umi = l_umi_uncorrected
                            ''' cluster UMI sequences for each cell barcodes '''
                            dict_index = _index_array( l_cb ) # index cell barcode values
                            for cb in dict_index :
                                l_index = dict_index[ cb ] # retrieve indices of the entries for the current cell barcode
                                l_umi_corrected[ l_index ] = _cluster_umi( l_umi_uncorrected[ l_index ] )
                        return l_cb, l_umi_corrected, l_umi_uncorrected

                    def _write_processed_bucket( res, l_l ) :
                        """ # 2023-08-12 16:11:02 
                        write the result of the processed bucket 
                        """
                        l_cb, l_umi_corrected, l_umi_uncorrected = res # parse the result
                        ns[ 'int_total_num_records_processed' ] += len( l_l ) # update 'int_total_num_records_processed'
                        for str_cb, str_umi_corrected, str_umi_uncorrected, t_r_and_tags in zip( l_cb, l_umi_corrected, l_umi_uncorrected, l_l ) : # save the sam records to the file
                            # parse records
                            r, dict_tags_existing = t_r_and_tags
                            if isinstance( str_cb, str ) : # if valid cell barcode was assigned, add tags to the reads
                                r.set_tags( [ ( 'CB', str_cb, 'Z' ), ( 'UB', str_umi_corrected, 'Z' ), ( 'UR', str_umi_uncorrected, 'Z' ) ] )
                            newsamfile.write( r ) 

                    def _write_results_from_offloaded_works( ) :
                        """ # 2023-08-12 21:57:27 
                        """
                        for res in workers_for_bucket_processing.wait_all( flag_return_results = True ).values( ) : # wait for all submitted works to be completed, and retrieve results for each work
                            _write_processed_bucket( res[ 'result' ], res[ 'associated_data' ] ) # save the sam records to the file
                            
                    def _empty_bucket( t_poly_a_site ) :
                        """ # 2023-08-10 21:21:29 
                        empty bucket for the 't_poly_a_site' by clustering UMI of the reads of the bucket and write the reads to the output BAM file
                        """
                        l_l = ns[ 'dict_poly_a_site_to_l_l' ].pop( t_poly_a_site ) # remove the bucket
                        ns[ 'int_bucket_deletion_count' ] += 1 # increase the counter
                        
                        if t_poly_a_site in set_e : # ❤️
                            logger.warn( f"{t_poly_a_site} already processed!" ) # ❤️
                        set_e.add( t_poly_a_site ) # ❤️
                        
                        ''' if the number of deletion count exceed the deletion count, re-initialize the bucket container '''
                        if ns[ 'int_bucket_deletion_count' ] > int_max_bucket_deletion_count_before_reinitialize :
                            dict_poly_a_site_to_l_l = dict( ) # create a new dictionary for reinitialization
                            for e in ns[ 'dict_poly_a_site_to_l_l' ] : # copy the dictionary
                                dict_poly_a_site_to_l_l[ e ] = ns[ 'dict_poly_a_site_to_l_l' ][ e ]
                            ns[ 'dict_poly_a_site_to_l_l' ] = dict_poly_a_site_to_l_l # reinitialize the container
                            ns[ 'int_bucket_deletion_count' ] = 0 # reset the counter

                        ''' correct CB / cluster UMI sequences and save results to the file '''
                        l_seq_cb_umi = list( dict_tags_existing[ 'CX' ] for r, dict_tags_existing in l_l ) # retrieve 'l_seq_cb_umi'
                        int_num_reads_in_a_bucket = len( l_l ) # retrieve the number of reads in a bucket
                        if int_num_reads_in_a_bucket < int_min_num_reads_in_a_bucket_for_parallel_processing : # 
                            _write_processed_bucket( _process_bucket( l_seq_cb_umi ), l_l ) # process a bucket in the current process, without offloading, and save the sam records to the file
                        else :
                            if not workers_for_bucket_processing.is_worker_available : # if all workers are working, wait for a while until workers are idle
                                _write_results_from_offloaded_works( ) # flush results from offloaded computations
                            workers_for_bucket_processing.submit_work( _process_bucket, args = ( l_seq_cb_umi, ), associated_data = l_l ) # submit the work for offloading (append 'l_r' as the data associated with the work)
                    
                    with pysam.AlignmentFile( path_file_bam_preprocessed, 'rb' ) as samfile :
                        for r in samfile.fetch( name_contig ) : # analyze all reads (since the pre-processed BAM file only contains valid reads that are already filtered) for the given chromosome
                            seq, cigartuples, flags, reference_start, reference_name = r.seq, r.cigartuples, r.flag, r.reference_start, r.reference_name # retrieve attributes
                            ''' process reads for each 'bucket' (reads with the same poly A tail attachment sites) '''
                            ''' when the contig has changed, empty all buckets '''
                            if reference_name_current != reference_name :  # retrieve a flag for emptying the buckets (when the contig changes)
                                for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                                    _empty_bucket( t_poly_a_site )
                                reference_name_current = reference_name # update the current contig name
                                
                            ''' determine whether to empty bucket or not, based on the current position on the sorted BAM file '''
                            for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                                ''' whether poly A site is located at the left or the right side of the read, when the current position passes the poly A site, process the bucket '''
                                flag_is_reverse_complemented, pos = t_poly_a_site # parse 't_poly_a_site'
                                if pos < reference_start :
                                    _empty_bucket( t_poly_a_site )
                                
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
                            if 'CX' in dict_tags_existing : # if cb_umi sequence is present
                                ''' retrieve poly (A) tail attachment site (specific definition: the alignment 'end' position that are closer to the poly (A) tail) '''
                                t_poly_a_site = ( flag_is_reverse_complemented, ( r.reference_start if flag_is_reverse_complemented else r.reference_end ) ) # retrieve a tuple indicating the aligned direction and poly A tail attachment position (alignment end position closer to the identified poly A tail)
                                if t_poly_a_site not in ns[ 'dict_poly_a_site_to_l_l' ] : # initialize 'dict_poly_a_site_to_l_l' for 't_poly_a_site'
                                    ns[ 'dict_poly_a_site_to_l_l' ][ t_poly_a_site ] = [ ]
                                ns[ 'dict_poly_a_site_to_l_l' ][ t_poly_a_site ].append( [ r, dict_tags_existing ] )
                            else :
                                ''' write the SAM record (record that does not contain the cell barcode - UMI sequence) ''' 
                                ns[ 'int_total_num_records_processed' ] += 1
                                newsamfile.write( r ) # write the record to the output BAM file
                            
                        ''' when all reads of the contig were read, empty all buckets '''
                        for t_poly_a_site in list( ns[ 'dict_poly_a_site_to_l_l' ] ) : # retrieve list of 't_poly_a_site'
                            _empty_bucket( t_poly_a_site )
                        _write_results_from_offloaded_works( ) # flush results from offloaded computations

                    """ report a batch has been completed """
                    pipe_sender.send( { 
                        'int_total_num_records_for_a_batch' : ns[ 'int_total_num_records_processed' ], # record the actual number of records processed for the batch
                    } )  # report the number of processed records

                """ close output files """
                newsamfile.close( )
                # sort the output sam file
                path_file_bam_barcoded_sorted = f"{path_folder_temp}{str_uuid}.barcoded.sorted.bam"
                pysam.sort( "-o", path_file_bam_barcoded_sorted, '-@', str( min( n_threads, 5 ) ), path_file_bam_barcoded )
                # index the resulting BAM file
                pysam.index( path_file_bam_barcoded_sorted )
                
                """ report the worker has completed all works """
                pipe_sender.send( 'completed' )  
                if verbose:
                    logger.info(f"[Completed] all works completed (worker_id={str_uuid})")

            ns = { 'int_num_read_currently_processed' : 0 }  # define a namespace 

            def post_process_batch(res):
                # update data using the received result
                ns["int_num_read_currently_processed"] += res[ 'int_total_num_records_for_a_batch' ]
                logger.info(
                    f"[{path_file_bam_input}] total {ns[ 'int_num_read_currently_processed' ]} number of reads has been processed."
                )  # report
                    
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

            """ 
            post-processing
            """

            def post_processing():  # off-loading a single-core work
                logger.info(
                    f"[{path_file_bam_input}] post-processing started"
                )
                
                # combine results into a single output file (initial read analysis)
                """ combine results into a single output BAM file """
                pysam.merge( '--threads', str( min( n_threads, 10 ) ), '-c', '-p', f"{path_folder_output}barcoded.bam", * glob.glob( f"{path_folder_temp}*.barcoded.sorted.bam" ) ) # merge output BAM files
                pysam.index( f"{path_folder_output}barcoded.bam" ) # index the input BAM file
                
                # write a flag indicating that the processing has been completed
                with open( f"{path_folder_output}pipeline_completed.txt", 'w' ) as newfile :
                    newfile.write( 'completed' )

                # delete temporary files
                # shutil.rmtree( path_folder_temp, ignore_errors = True )
                    
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


if __name__ == "__main__":
    ourotools()  # run ouro at the top-level environment
