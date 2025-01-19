from .biobookshelf import *
from . import biobookshelf as bk

def Get_contig_names_from_bam_header( path_file_bam : str, mode : str = 'rb' ) : 
    """ # 2023-08-10 22:13:00 
    retrieve contig names from the BAM header. To read a SAM file, change 'mode' to 'r'
    """
    try:
        import pysam
    except ImportError as e:
        e.add_note( f"Please install `pysam` and try again." )
        raise
    with pysam.AlignmentFile( path_file_bam, mode ) as samfile :
        l_contig_names = list( e[ 'SN' ] for e in samfile.header.to_dict( )[ 'SQ' ] )
    return l_contig_names

def Generate_Base_and_Qual(r):
    """# 2021-08-27 22:09:31
    a generator function yielding the base and quality of every matched/mismatched positions of a given pysam read record
    """
    pos_read, pos_ref = (
        0,
        0,
    )  # current position in the extracted reference sequence and the current read
    (
        int_oper_M,
        int_oper_I,
        int_oper_D,
        int_oper_N,
        int_oper_S,
        int_oper_H,
        int_oper_P,
        int_oper_equal,
        int_oper_X,
    ) = list(range(9))
    for int_oper, n_bases in r.cigartuples:
        if int_oper == int_oper_M:
            for _ in range(n_bases):
                str_base_read, str_qual_read = (
                    r.seq[pos_read],
                    r.query_qualities[pos_read],
                )
                yield r.reference_start + pos_ref, str_base_read, str_qual_read
                pos_ref += 1
                pos_read += 1
        elif int_oper == int_oper_N:
            pos_ref += n_bases
        elif int_oper == int_oper_S:
            pos_read += n_bases
        elif int_oper == int_oper_I:
            pos_read += n_bases
        elif int_oper == int_oper_D:
            pos_ref += n_bases
        elif int_oper == int_oper_H:
            pass
        elif int_oper == int_oper_P:
            pass
        elif int_oper == int_oper_equal:
            str_base_read, str_qual_read = r.seq[pos_read], r.query_qualities[pos_read]
            yield r.reference_start + pos_ref, str_base_read, str_qual_read
            pos_ref += n_bases
            pos_read += n_bases
        elif int_oper == int_oper_X:
            for _ in range(n_bases):
                str_base_read, str_qual_read = (
                    r.seq[pos_read],
                    r.query_qualities[pos_read],
                )
                yield r.reference_start + pos_ref, str_base_read, str_qual_read
                pos_ref += 1
                pos_read += 1


def Retrieve_List_of_Mapped_Segments(
    cigartuples,
    pos_start,
    return_1_based_coordinate=False,
    flag_pos_start_0_based_coord=True,
    flag_return_splice_junction=False,
    flag_is_cigartuples_from_mappy=False,
):
    """# 2024-01-06 20:04:15 
    return l_seq and int_total_aligned_length for given cigartuples (returned by pysam cigartuples) and 'pos_start' (0-based coordinates, assuming pos_start is 0-based coordinate)
    'return_1_based_coordinate' : return 1-based coordinate, assuming 'pos_start' is 0-based coordinate (pysam returns 0-based coordinate)
    'flag_return_splice_junction' : additionally return list of splice junction tuples
    'flag_is_cigartuples_from_mappy' : if cigartuples are from mappy.Alignment, set this flag to True
    """
    l_seg, start, int_aligned_length, int_total_aligned_length = list(), pos_start, 0, 0
    if return_1_based_coordinate and flag_pos_start_0_based_coord:  # 0-based > 1-based
        start += 1
    for operation, length in cigartuples:
        if (
            flag_is_cigartuples_from_mappy
        ):  # if flag_is_cigartuples_from_mappy, swap the two values
            temp = operation
            operation = length
            length = temp
        if operation in {0, 2, 7, 8}:  # 'MD=X'
            int_aligned_length += length
        elif (
            operation == 3
        ):  # 'N' if splice junction appears, split the region and make a record
            l_seg.append(
                (
                    start,
                    (start + int_aligned_length - 1)
                    if return_1_based_coordinate
                    else (start + int_aligned_length),
                )
            )  # set the end position
            start = start + int_aligned_length + length  # set the next start position
            int_total_aligned_length += (
                int_aligned_length  # update total aligned length
            )
            int_aligned_length = 0
    if int_aligned_length > 0:
        l_seg.append(
            (
                start,
                (start + int_aligned_length - 1)
                if return_1_based_coordinate
                else (start + int_aligned_length),
            )
        )
        int_total_aligned_length += int_aligned_length
    if flag_return_splice_junction:
        # collect splice junction tuples from l_seg
        l_sj = []
        for i in range(len(l_seg) - 1):
            l_sj.append((l_seg[i][1], l_seg[i + 1][0]))
        return l_seg, int_total_aligned_length, l_sj
    else:
        return l_seg, int_total_aligned_length
Retrive_List_of_Mapped_Segments = Retrieve_List_of_Mapped_Segments

def Call_Variant(r, dict_fasta_genome, function_for_processing_reference_name=None):
    """# 2023-01-07 18:58:43
    perform vcf-style variant calling of a single aligned read (pysam read object) using a given genome
    return a list of mutations and corrected read sequence

    'r' : AlignedSegment object returned by pysam
    'dict_fasta_genome' : dictionary of genome sequences
    'function_for_processing_reference_name' : a function that can be applied to reference name to make reference_name used in r consistent with those used in 'dict_fasta_genome' and the id_mut that will be returned by the current function. (e.g. a function that can remove 'chr' prefix from the reference name, if chromosome name without the 'chr' prefix is being used)

    --- returns ---
    'l_mut' : list of variants with the following nomenclature
    {refname}:{refpos}:{ref}>{alt}, where refpos, ref, alt follows nomenclature of VCF file
    """
    pos_read, pos_ref = (
        0,
        0,
    )  # current position in the extracted reference sequence and the current read

    # parse a sam record
    arr_qual = r.query_qualities
    seq = r.seq
    ref_name = r.reference_name
    if function_for_processing_reference_name is not None:
        ref_name = function_for_processing_reference_name(ref_name)
    ref_start = r.reference_start
    cigartuples = r.cigartuples
    alen = r.alen

    # retrieve a part of the reference sequence where the current read was aligned
    seq_ref = dict_fasta_genome[ref_name][ref_start : ref_start + alen]

    l_mut = []

    """
    # define interger representation of the CIGAR operations used in BAM files
    
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    """
    int_cigarop_M = 0
    int_cigarop_I = 1
    int_cigarop_D = 2
    int_cigarop_N = 3
    int_cigarop_S = 4
    int_cigarop_H = 5
    int_cigarop_P = 6
    int_cigarop_equal = 7
    int_cigarop_X = 8

    ns = dict()  # create a namespace
    # initialilze the namespace
    ns[
        "pos_ref_variant_start"
    ] = None  # 0-based coordinate of the start of the alternative allele on the reference
    ns[
        "pos_ref_variant_end"
    ] = None  # 0-based coordinate of the end of the alternative allele on the reference
    ns["alt"] = ""

    def __update_alt(pos_ref: int, len_bases_ref: int = 0, bases_alt: str = ""):
        """# 2023-01-07 17:51:14"""
        # initialize variant
        if ns["pos_ref_variant_start"] is None:
            ns["pos_ref_variant_start"] = pos_ref
            ns["pos_ref_variant_end"] = pos_ref
            ns["alt"] = ""
        ns["pos_ref_variant_end"] += len_bases_ref
        ns["alt"] += bases_alt

    def __flush_alt():
        """# 2023-01-07 17:24:20
        flush variant allele
        """
        if ns["pos_ref_variant_start"] is not None:
            if (
                len(ns["alt"]) == 0
                or ns["pos_ref_variant_end"] == ns["pos_ref_variant_start"]
            ):  # for simple insertion/deletion variants, add a single base pair of the reference before insertion or deletion to record the variant
                ns["pos_ref_variant_start"] -= 1
                ns["alt"] = seq_ref[ns["pos_ref_variant_start"]] + ns["alt"]
            l_mut.append(
                f"{ref_name}:{ref_start + 1 + ns[ 'pos_ref_variant_start' ]}:{seq_ref[ ns[ 'pos_ref_variant_start' ] : ns[ 'pos_ref_variant_end' ] ]}>{ns[ 'alt' ]}"
            )  # append the variant
            ns["pos_ref_variant_start"] = None
            ns["pos_ref_variant_end"] = None
            ns["alt"] = ""

    for int_oper, n_bases in cigartuples:
        if int_oper == int_cigarop_M:
            for _ in range(n_bases):
                str_base_read, str_base_ref, str_qual_read = (
                    seq[pos_read],
                    seq_ref[pos_ref],
                    arr_qual[pos_read],
                )
                if str_base_read == str_base_ref:
                    __flush_alt()  # flush a variant
                else:
                    __update_alt(pos_ref, len_bases_ref=1, bases_alt=str_base_read)
                pos_ref += 1
                pos_read += 1
        elif int_oper == int_cigarop_N:
            __flush_alt()  # flush a variant
            pos_ref += n_bases
        elif int_oper == int_cigarop_S:
            __flush_alt()  # flush a variant
            pos_read += n_bases
        elif int_oper == int_cigarop_I:
            __update_alt(pos_ref, bases_alt=seq[pos_read : pos_read + n_bases])
            pos_read += n_bases
        elif int_oper == int_cigarop_D:
            __update_alt(pos_ref, len_bases_ref=n_bases)
            pos_ref += n_bases
        elif int_oper == int_cigarop_H:
            __flush_alt()  # flush a variant
            pass
        elif int_oper == int_cigarop_P:
            __flush_alt()  # flush a variant
            pass
        elif int_oper == int_cigarop_equal:
            __flush_alt()  # flush a variant
            pos_ref += n_bases
            pos_read += n_bases
        elif int_oper == int_cigarop_X:
            for _ in range(n_bases):
                str_base_read, str_qual_read = seq[pos_read], arr_qual[pos_read]
                __update_alt(pos_ref, 1, str_base_read)
                pos_ref += 1
                pos_read += 1
    __flush_alt()  # flush a variant
    return l_mut


def Retrieve_Variant(
    r,
    dict_fasta_genome,
    set_mut_filtered=None,
    return_corrected_read_sequence=False,
    flag_ignore_indel=False,
    flag_return_as_string=True,
    flag_return_matched=False,
    function_for_processing_reference_name=None,
):
    """# 2022-01-15 00:37:28
    perform variant calling of a single aligned read (pysam read object) using a given genome
    return a list of mutations and corrected read sequence

    'r' : AlignedSegment object returned by pysam
    'dict_fasta_genome' : dictionary of genome sequences
    'set_mut_filtered' : only consider mutations in the given 'set_mut_filtered'. only valid when 'return_corrected_read_sequence' = True
    'return_corrected_read_sequence' : return aligned read sequence after hard-clipping and filtering variants using 'set_mut_filtered'
    'flag_ignore_indel' : ignore insertion and deletion mutation calls
    'flag_return_as_string' : return a called mutation as a string in the following format: f"{reference_name}:{start}_{cigar_operation}_{mutation}". If set to False, return the called mutation as a tuple (reference_name, start, int_operation, mutation, quality)
    'flag_return_matched' : return 'id_mut' even when there is no mutation. the 'int_operation' for matched base is int_cigarop_M
    'function_for_processing_reference_name' : a function that can be applied to reference name to make reference_name used in r consistent with those used in 'dict_fasta_genome' and the id_mut that will be returned by the current function. (e.g. a function that can remove 'chr' prefix from the reference name)

    --- returns ---
    'l_mut' : list of variants with the following nomenclature
    # substitution:
    {refname}:{refpos}_{refbase}>{mutbase}

    # insertion:
    {refname}:{refpos_after_the_insertion}_ins_{inserted_bases}

    # deletion:
    {refname}:{refpos_of_the_first_base_of_the_deletion}_del_{number_of_deleted_bases}
    """
    pos_read, pos_ref = (
        0,
        0,
    )  # current position in the extracted reference sequence and the current read

    # parse a sam record
    arr_qual = r.query_qualities
    seq = r.seq
    ref_name = r.reference_name
    if function_for_processing_reference_name is not None:
        ref_name = function_for_processing_reference_name(ref_name)
    ref_start = r.reference_start
    cigartuples = r.cigartuples
    alen = r.alen

    # retrieve a part of the reference sequence where the current read was aligned
    seq_ref = dict_fasta_genome[ref_name][ref_start : ref_start + alen]

    l_mut = []

    """
    # define interger representation of the CIGAR operations used in BAM files
    
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    """
    int_cigarop_M = 0
    int_cigarop_I = 1
    int_cigarop_D = 2
    int_cigarop_N = 3
    int_cigarop_S = 4
    int_cigarop_H = 5
    int_cigarop_P = 6
    int_cigarop_equal = 7
    int_cigarop_X = 8

    if return_corrected_read_sequence:  # correction mode
        str_seq_corrected_read = (
            ""  # corrected read sequence after hard clipping and filtering variants
        )
        for int_oper, n_bases in cigartuples:
            if int_oper == int_cigarop_M:
                for _ in range(n_bases):
                    str_base_read, str_base_ref, str_qual_read = (
                        seq[pos_read],
                        seq_ref[pos_ref],
                        arr_qual[pos_read],
                    )
                    str_base_read_corrected = (
                        str_base_ref  # default corrected read base = ref
                    )
                    if str_base_read != str_base_ref:
                        id_mut = (
                            f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}"
                            if flag_return_as_string
                            else (
                                ref_name,
                                ref_start + pos_ref,
                                "X",
                                str_base_read,
                                str_qual_read,
                            )
                        )  # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        if set_mut_filtered is None or id_mut in set_mut_filtered:
                            l_mut.append(id_mut)
                            str_base_read_corrected = str_base_read
                    elif flag_return_matched:
                        id_mut = (
                            f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}"
                            if flag_return_as_string
                            else (
                                ref_name,
                                ref_start + pos_ref,
                                "M",
                                str_base_read,
                                str_qual_read,
                            )
                        )  # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        l_mut.append(id_mut)
                    str_seq_corrected_read += (
                        str_base_read_corrected  # add corrected read
                    )
                    pos_ref += 1
                    pos_read += 1
            elif int_oper == int_cigarop_N:
                pos_ref += n_bases
            elif int_oper == int_cigarop_S:
                pos_read += n_bases
            elif int_oper == int_cigarop_I:
                id_mut = (
                    f"{ref_name}:{ref_start + 1 + pos_ref}_ins_{seq[ pos_read : pos_read + n_bases ]}"
                    if flag_return_as_string
                    else (
                        ref_name,
                        ref_start + pos_ref,
                        "I",
                        seq[pos_read : pos_read + n_bases],
                        arr_qual[pos_read : pos_read + n_bases],
                    )
                )  # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and (
                    set_mut_filtered is None or id_mut in set_mut_filtered
                ):
                    l_mut.append(id_mut)
                    str_seq_corrected_read += seq[
                        pos_read : pos_read + n_bases
                    ]  # add inserted sequence (if insertion is valid)
                pos_read += n_bases
            elif int_oper == int_cigarop_D:
                id_mut = (
                    f"{ref_name}:{ref_start + 1 + pos_ref}_del_{n_bases}"
                    if flag_return_as_string
                    else (ref_name, ref_start + pos_ref, "D", n_bases, None)
                )  # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and (
                    set_mut_filtered is None or id_mut in set_mut_filtered
                ):
                    l_mut.append(id_mut)
                else:
                    str_seq_corrected_read += seq_ref[
                        pos_ref : pos_ref + n_bases
                    ]  # add deleted reference sequence (if deletion is invalid)
                pos_ref += n_bases
            elif int_oper == int_cigarop_H:
                pass
            elif int_oper == int_cigarop_P:
                pass
            elif int_oper == int_cigarop_equal:
                if flag_return_matched:
                    str_base_read, str_qual_read = seq[pos_read], arr_qual[pos_read]
                    id_mut = (
                        f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}"
                        if flag_return_as_string
                        else (
                            ref_name,
                            ref_start + pos_ref,
                            "M",
                            str_base_read,
                            str_qual_read,
                        )
                    )  # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append(id_mut)
                str_seq_corrected_read += seq[
                    pos_read : pos_read + n_bases
                ]  # add read sequences
                pos_ref += n_bases
                pos_read += n_bases
            elif int_oper == int_cigarop_X:
                for _ in range(n_bases):
                    str_base_read, str_base_ref, str_qual_read = (
                        seq[pos_read],
                        seq_ref[pos_ref],
                        arr_qual[pos_read],
                    )
                    str_base_read_corrected = (
                        str_base_ref  # default corrected read base = ref
                    )
                    id_mut = (
                        f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}"
                        if flag_return_as_string
                        else (
                            ref_name,
                            ref_start + pos_ref,
                            "X",
                            str_base_read,
                            str_qual_read,
                        )
                    )  # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    if set_mut_filtered is None or id_mut in set_mut_filtered:
                        l_mut.append(id_mut)
                        str_base_read_corrected = str_base_read
                    str_seq_corrected_read += (
                        str_base_read_corrected  # add corrected read
                    )
                    pos_ref += 1
                    pos_read += 1
        return l_mut, str_seq_corrected_read
    else:  # normal mutation calling mode
        for int_oper, n_bases in cigartuples:
            if int_oper == int_cigarop_M:
                for _ in range(n_bases):
                    str_base_read, str_base_ref, str_qual_read = (
                        seq[pos_read],
                        seq_ref[pos_ref],
                        arr_qual[pos_read],
                    )
                    if str_base_read != str_base_ref:
                        l_mut.append(
                            f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}"
                            if flag_return_as_string
                            else (
                                ref_name,
                                ref_start + pos_ref,
                                "X",
                                str_base_read,
                                str_qual_read,
                            )
                        )  # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    elif flag_return_matched:
                        l_mut.append(
                            f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}"
                            if flag_return_as_string
                            else (
                                ref_name,
                                ref_start + pos_ref,
                                "M",
                                str_base_read,
                                str_qual_read,
                            )
                        )  # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
            elif int_oper == int_cigarop_N:
                pos_ref += n_bases
            elif int_oper == int_cigarop_S:
                pos_read += n_bases
            elif int_oper == int_cigarop_I:
                if not flag_ignore_indel:
                    l_mut.append(
                        f"{ref_name}:{ref_start + 1 + pos_ref}_ins_{seq[ pos_read : pos_read + n_bases ]}"
                        if flag_return_as_string
                        else (
                            ref_name,
                            ref_start + pos_ref,
                            "I",
                            seq[pos_read : pos_read + n_bases],
                            arr_qual[pos_read : pos_read + n_bases],
                        )
                    )  # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_read += n_bases
            elif int_oper == int_cigarop_D:
                if not flag_ignore_indel:
                    l_mut.append(
                        f"{ref_name}:{ref_start + 1 + pos_ref}_del_{n_bases}"
                        if flag_return_as_string
                        else (ref_name, ref_start + pos_ref, "D", n_bases, None)
                    )  # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_ref += n_bases
            elif int_oper == int_cigarop_H:
                pass
            elif int_oper == int_cigarop_P:
                pass
            elif int_oper == int_cigarop_equal:
                if flag_return_matched:
                    str_base_read, str_qual_read = seq[pos_read], arr_qual[pos_read]
                    id_mut = (
                        f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}"
                        if flag_return_as_string
                        else (
                            ref_name,
                            ref_start + pos_ref,
                            "M",
                            str_base_read,
                            str_qual_read,
                        )
                    )  # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append(id_mut)
                pos_ref += n_bases
                pos_read += n_bases
            elif int_oper == int_cigarop_X:
                for _ in range(n_bases):
                    str_base_read, str_qual_read = seq[pos_read], arr_qual[pos_read]
                    l_mut.append(
                        f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}"
                        if flag_return_as_string
                        else (
                            ref_name,
                            ref_start + pos_ref,
                            "X",
                            str_base_read,
                            str_qual_read,
                        )
                    )  # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
        return l_mut
