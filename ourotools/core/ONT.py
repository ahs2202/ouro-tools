# load internal module
from .biobookshelf import *
from . import biobookshelf as bk
from typing import Union, List, Literal, Dict, Callable, Set, Iterable, Tuple


def Minimap2_Align(
    path_file_fastq,
    path_file_minimap2_index="/node210data/shared/ensembl/Mus_musculus/index/minimap2/Mus_musculus.GRCm38.dna.primary_assembly.k_14.idx",
    path_folder_minimap2_output=None,
    n_threads=20,
    verbose=True,
    drop_unaligned=False,
    return_bash_shellscript=False,
    n_threads_for_sort=10,
    flag_use_split_prefix: bool = False,
    path_file_junc_bed: Union[
        None, str
    ] = None,  # if given, the bed file will be used for prioritizing known splice sites.
    path_file_gtf: Union[
        None, str
    ] = None,  # path to gene and exon annotation files, required if 'path_file_junc_bed' is given but the file does not exist
):
    """
    # 2023-04-23 01:18:58
    align given fastq file of nanopore reads using minimap2 and write an output as a bam file
    'path_file_fastq' : input fastq or fasta file (gzipped or uncompressed file is accepted)
    'path_file_minimap2_index' : minimap2 index file
    'path_folder_minimap2_output' : minimap2 output folder
    'drop_unaligned' : a flag indicating whether reads not aligned to the reference ('SAM flag == 4') are included in the output bam file
    'return_bash_shellscript' : return shellscript instead of running minimap2 using the subprocess module
    'flag_use_split_prefix' = False # for large index, split-prefix should be used
    """
    path_folder_fastq, name_file_fastq = path_file_fastq.rsplit("/", 1)
    if (
        path_folder_minimap2_output is None
    ):  # default output folder is a subdirectory of the folder containing the input fastq file
        path_folder_minimap2_output = f"{path_folder_fastq}/minimap2/"
    if (
        path_folder_minimap2_output[-1] != "/"
    ):  # add '/' at the end of the output directory if it does not exist
        path_folder_minimap2_output += "/"
    os.makedirs(
        path_folder_minimap2_output, exist_ok=True
    )  # create folder if it does not exist

    path_file_sam = (
        f"{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.sam"
    )
    path_file_bam = (
        f"{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.bam"
    )
    # if index file of the output BAM file exists, exit
    if os.path.exists( f"{path_file_bam}.bai" ) :
        return

    l_bash_shellscript = []

    """ perform minimap2 alignment """
    l_arg = [
        "minimap2",
        "-t",
        str(int(n_threads)),
        "-ax",
        "splice",
        "-o",
        path_file_sam,
    ]

    # for large index, split-prefix should be used
    if flag_use_split_prefix:
        l_arg += [f"--split-prefix={path_folder_minimap2_output}{UUID( )}"]

    if path_file_junc_bed is not None:
        if (
            not os.path.exists(path_file_junc_bed) and path_file_gtf is not None
        ):  # if the bed file does not exist, create the bed file using paftools.js, packaged with the minimap2 executable
            l_args_for_creating_junc_bed = ["paftools.js", "gff2bed", path_file_gtf]
            if (
                return_bash_shellscript
            ):  # perform minimap2 alignment using subprocess module
                l_bash_shellscript.append(
                    " ".join(l_args_for_creating_junc_bed + [">", path_file_junc_bed])
                )
            else:
                bk.OS_Run(
                    l_args_for_creating_junc_bed,
                    path_file_stdout=path_file_junc_bed,
                    stdout_binary=False,
                )
        if os.path.exists(path_file_junc_bed):
            l_arg += ["--junc-bed", path_file_junc_bed]

    if drop_unaligned:
        l_arg += ["--sam-hit-only"]
    l_arg += [path_file_minimap2_index, path_file_fastq]
    if return_bash_shellscript:  # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append(" ".join(l_arg))
    else:
        run = subprocess.run(l_arg, capture_output=True)
        with open(
            f"{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.out", "w"
        ) as file:
            file.write(run.stdout.decode())
        if verbose:
            print("minimap2 completed")

    """ sort output SAM file """
    l_arg = [
        "samtools",
        "sort",
        "-@",
        str(int(min(n_threads_for_sort, 10))),
        "-O",
        "BAM",
        "-o",
        path_file_bam,
        path_file_sam,
    ]
    if return_bash_shellscript:  # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append(" ".join(l_arg))
        l_bash_shellscript.append(" ".join(["rm", "-f", path_file_sam]))
    else:
        run = subprocess.run(l_arg, capture_output=False)
        os.remove(path_file_sam)  # remove sam file

    """ index resulting BAM file """
    l_arg = ["samtools", "index", path_file_bam]
    if return_bash_shellscript:  # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append(" ".join(l_arg))
    else:
        run = subprocess.run(l_arg, capture_output=False)
        if verbose:
            print("samtools bam file compressing and indexing completed")

    if return_bash_shellscript:  # retrun bash shell scripts
        return " && ".join(l_bash_shellscript)


def Minimap2_Index(path_file_fasta, path_file_minimap2_index=None, verbose=False):
    """
    # 2021-03-24 00:44:51
    index given fasta file for nanopore reads alignment
    'path_file_fasta' : input reference fasta file
    'path_file_minimap2_index' : minimap2 index file
    """
    path_folder_fastq, name_file_fasta = path_file_fasta.rsplit("/", 1)
    if (
        path_file_minimap2_index is None
    ):  # set the default directory of the minimap index
        path_file_minimap2_index = (
            f"{path_folder_fastq}/index/minimap2/{name_file_fasta}.ont.mmi"
        )
    path_folder_minimap2_index, name_file_index = path_file_minimap2_index.rsplit(
        "/", 1
    )
    path_folder_minimap2_index += "/"
    os.makedirs(
        path_folder_minimap2_index, exist_ok=True
    )  # create folder if it does not exist
    if os.path.exists(path_file_minimap2_index):  # exit if an index file already exists
        return
    # build minimap2 index
    run = subprocess.run(
        ["minimap2", "-x", "map-ont", "-d", path_file_minimap2_index, path_file_fasta],
        capture_output=True,
    )

    with open(
        f"{path_folder_minimap2_index}{name_file_index}.minimap2_index.out", "w"
    ) as file:
        file.write(run.stdout.decode())
    if verbose:
        print("minimap2 indexing completed")
