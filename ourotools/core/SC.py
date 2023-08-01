from .biobookshelf import *
from . import biobookshelf as bk


def MTX_10X_Feature_add_prefix_or_suffix(
    path_file_features_input,
    path_file_features_output=None,
    id_feature_prefix="",
    id_feature_suffix="",
    name_feature_prefix="",
    name_feature_suffix="",
):
    """# 2022-04-21 21:10:21
    Add prefix or suffix to the id_feature and name_feature of a given 'features.tsv.gz' file
    'path_file_features_output' : default: None. by default, the input 'path_file_features_input' file will be overwritten with the modified features
    """
    flag_replace_input_file = (
        path_file_features_output is None
    )  # retrieve a flag indicating the replacement of original input file with modified input file
    if flag_replace_input_file:
        path_file_features_output = (
            f"{path_file_features_input}.writing.tsv.gz"  # set a temporary output file
        )
    newfile = gzip.open(path_file_features_output, "wb")  # open an output file
    with gzip.open(path_file_features_input, "rb") as file:
        while True:
            line = file.readline()
            if len(line) == 0:
                break
            id_feature, name_feature, type_feature = (
                line.decode().strip().split("\t")
            )  # parse a feature
            id_feature_new = (
                id_feature_prefix + id_feature + id_feature_suffix
            )  # compose new id_feature
            name_feature_new = (
                name_feature_prefix + name_feature + name_feature_suffix
            )  # compose new name_feature
            newfile.write(
                (
                    "\t".join([id_feature_new, name_feature_new, type_feature]) + "\n"
                ).encode()
            )  # write a new feature
    newfile.close()  # close the output file
    # if the output file path was not given, replace the original file with modified file
    if flag_replace_input_file:
        os.remove(path_file_features_input)
        os.rename(path_file_features_output, path_file_features_input)


dict_id_entry_to_index_entry = dict()


def __MTX_10X_Combine__renumber_barcode_or_feature_index_mtx_10x__(
    path_file_input, path_folder_mtx_10x_output, flag_renumber_feature_index
):
    """
    internal function for MTX_10X_Combine
    # 2022-04-21 12:10:53

    'flag_renumber_feature_index' : if True, assumes barcodes are not shared between matrices and renumber features only. If False, assumes features are not shared between matrices and renumber barcodes only.
    """
    global dict_id_entry_to_index_entry
    for (
        path_folder_mtx_10x,
        int_total_n_entries_of_previously_written_matrices,
        index_mtx_10x,
    ) in pd.read_csv(path_file_input, sep="\t").values:
        # directly write matrix.mtx.gz file without header
        with gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        ) as newfile:
            arr_id_entry = pd.read_csv(
                f"{path_folder_mtx_10x}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
                sep="\t",
                header=None,
            ).values[
                :, 0
            ]  # retrieve a list of id_feature for the current dataset
            with gzip.open(
                f"{path_folder_mtx_10x}matrix.mtx.gz", "rb"
            ) as file:  # retrieve a list of features
                line = file.readline().decode()  # read the first line
                # if the first line of the file contains a comment line, read all comment lines and a description line following the comments.
                if len(line) > 0 and line[0] == "%":
                    # read comment and the description line
                    while True:
                        if line[0] != "%":
                            break
                        line = file.readline().decode()  # read the next line
                    line = (
                        file.readline().decode()
                    )  # discard the description line and read the next line
                # process entries
                while True:
                    if len(line) == 0:
                        break
                    index_row, index_col, int_value = tuple(
                        map(int, line.strip().split())
                    )  # parse each entry of the current matrix

                    newfile.write(
                        (
                            " ".join(
                                tuple(
                                    map(
                                        str,
                                        (
                                            [
                                                dict_id_entry_to_index_entry[
                                                    arr_id_entry[index_row - 1]
                                                ],
                                                index_col
                                                + int_total_n_entries_of_previously_written_matrices,
                                            ]
                                            if flag_renumber_feature_index
                                            else [
                                                index_row
                                                + int_total_n_entries_of_previously_written_matrices,
                                                dict_id_entry_to_index_entry[
                                                    arr_id_entry[index_col - 1]
                                                ],
                                            ]
                                        )
                                        + [int_value],
                                    )
                                )
                            )
                            + "\n"
                        ).encode()
                    )  # translate indices of the current matrix to that of the combined matrix
                    line = file.readline().decode()  # read the next line


def MTX_10X_Combine(
    path_folder_mtx_10x_output,
    *l_path_folder_mtx_10x_input,
    int_num_threads=15,
    flag_split_mtx=True,
    flag_split_mtx_again=False,
    int_max_num_entries_for_chunk=10000000,
    flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs=None,
    flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs=None,
    verbose=False,
):
    """
    # 2022-02-22 00:38:36
    Combine 10X count matrix files from the given list of folders and write combined output files to the given output folder 'path_folder_mtx_10x_output'
    If there are no shared cells between matrix files, a low-memory mode will be used. The output files will be simply combined since no count summing operation is needed. Only feature matrix will be loaded and updated in the memory.
    'id_feature' should be unique across all features

    'int_num_threads' : number of threads to use when combining datasets. multiple threads will be utilized only when datasets does not share cells and thus can be safely concatanated.
    'flag_split_mtx' : split the resulting mtx file so that the contents in the output mtx file can be processed in parallel without ungzipping the mtx.gz file and spliting the file.
    'flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs' : a flag for entering low-memory mode when there is no shared cells between given input matrices. By default (when None is given), matrices will be examined and the flag will be set automatically by the program. To reduce running time and memory, this flag can be manually set by users.
    'flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs' : a flag for entering low-memory mode when there is no shared features between given input matrices. By default (when None is given), matrices will be examined and the flag will be set automatically by the program. To reduce running time and memory, this flag can be manually set by users.
    """

    # create an output folder
    os.makedirs(path_folder_mtx_10x_output, exist_ok=True)

    if flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs is None:
        """retrieve cell barcodes of all 10X matrices and check whether cell barcodes are not shared between matrices"""
        int_total_n_barcodes_of_previously_written_matrices = (
            0  # follow the number of barcodes that are previously written
        )
        l_int_total_n_barcodes_of_previously_written_matrices = (
            []
        )  # calculate the number of barcodes of the previous dataset in the combined mtx.
        set_barcode = set()  # update a set of unique barcodes
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            arr_barcode = (
                pd.read_csv(
                    f"{path_folder_mtx_10x}barcodes.tsv.gz", sep="\t", header=None
                )
                .squeeze("columns")
                .values
            )  # retrieve a list of features
            set_barcode.update(arr_barcode)  # update a set of barcodes
            l_int_total_n_barcodes_of_previously_written_matrices.append(
                int_total_n_barcodes_of_previously_written_matrices
            )
            int_total_n_barcodes_of_previously_written_matrices += len(
                arr_barcode
            )  # update the number of barcodes
        """ check whether there are shared cell barcodes between matrices and set a flag for entering a low-memory mode """
        flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs = (
            len(set_barcode) == int_total_n_barcodes_of_previously_written_matrices
        )  # update flag

    if (
        not flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs
        and flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs is None
    ):
        """retrieve features of all 10X matrices and check whether features are not shared between matrices"""
        int_total_n_features_of_previously_written_matrices = (
            0  # follow the number of features that are previously written
        )
        l_int_total_n_features_of_previously_written_matrices = (
            []
        )  # calculate the number of features of the previous dataset in the combined mtx.
        set_feature = set()  # update a set of unique features
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            arr_feature = (
                pd.read_csv(
                    f"{path_folder_mtx_10x}features.tsv.gz",
                    sep="\t",
                    header=None,
                    usecols=[0],
                )
                .squeeze("columns")
                .values
            )  # retrieve a list of features
            set_feature.update(arr_feature)  # update a set of features
            l_int_total_n_features_of_previously_written_matrices.append(
                int_total_n_features_of_previously_written_matrices
            )
            int_total_n_features_of_previously_written_matrices += len(
                arr_feature
            )  # update the number of features
        """ check whether there are shared features between matrices and set a flag for entering a low-memory mode """
        flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs = (
            len(set_feature) == int_total_n_features_of_previously_written_matrices
        )  # update flag

    flag_low_memory_mode = (
        flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs
        or flag_low_memory_mode_because_there_is_no_shared_feature_between_mtxs
    )  # retrieve flag for low-memory mode
    if flag_low_memory_mode:
        """low-memory mode"""
        flag_renumber_feature_index = flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs  # retrieve a flag for renumbering features
        if verbose:
            print(
                f"entering low-memory mode, re-numbering {'features' if flag_renumber_feature_index else 'barcodes'} index because {'barcodes' if flag_renumber_feature_index else 'features'} are not shared across the matrices."
            )

        """ write a combined barcodes/features.tsv.gz - that are not shared between matrices """
        OS_FILE_Combine_Files_in_order(
            list(
                f"{path_folder_mtx_10x}{'barcodes' if flag_renumber_feature_index else 'features'}.tsv.gz"
                for path_folder_mtx_10x in l_path_folder_mtx_10x_input
            ),
            f"{path_folder_mtx_10x_output}{'barcodes' if flag_renumber_feature_index else 'features'}.tsv.gz",
            overwrite_existing_file=True,
        )

        """ collect a set of unique entries and a list of entries for each 10X matrix """
        set_t_entry = set()  # update a set unique id_entry (either id_cell or id_entry)
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            set_t_entry.update(
                list(
                    map(
                        tuple,
                        pd.read_csv(
                            f"{path_folder_mtx_10x}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
                            sep="\t",
                            header=None,
                        ).values,
                    )
                )
            )  # update a set of feature tuples

        """ write a combined features/barcodes.tsv.gz - that are shared between matrices """
        l_t_entry = list(set_t_entry)  # convert set to list
        with gzip.open(
            f"{path_folder_mtx_10x_output}{'features' if flag_renumber_feature_index else 'barcodes'}.tsv.gz",
            "wb",
        ) as newfile:
            for t_entry in l_t_entry:
                newfile.write(("\t".join(t_entry) + "\n").encode())

        """ build a mapping of id_entry to index_entry, which will be consistent across datasets - for features/barcodes that are shared between matrices """
        global dict_id_entry_to_index_entry  # use global variable for multiprocessing
        dict_id_entry_to_index_entry = dict(
            (t_entry[0], index_entry + 1)
            for index_entry, t_entry in enumerate(l_t_entry)
        )  # 0>1 based index
        PICKLE_Write(
            f"{path_folder_mtx_10x_output}dict_id_entry_to_index_entry.pickle",
            dict_id_entry_to_index_entry,
        )  # save id_feature to index_feature mapping as a pickle file

        """ collect the number of records for each 10X matrix """
        int_total_n_records = 0
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            with gzip.open(
                f"{path_folder_mtx_10x}matrix.mtx.gz", "rb"
            ) as file:  # retrieve a list of features
                file.readline(), file.readline()
                int_total_n_records += int(
                    file.readline().decode().strip().split()[2]
                )  # update the total number of entries

        """ write a part of a combined matrix.mtx.gz for each dataset using multiple processes """
        # compose inputs for multiprocessing
        df_input = pd.DataFrame(
            {
                "path_folder_input_mtx_10x": l_path_folder_mtx_10x_input,
                "int_total_n_barcodes_of_previously_written_matrices": (
                    l_int_total_n_barcodes_of_previously_written_matrices
                    if flag_renumber_feature_index
                    else l_int_total_n_features_of_previously_written_matrices
                ),
                "index_mtx_10x": np.arange(
                    len(l_int_total_n_barcodes_of_previously_written_matrices)
                    if flag_renumber_feature_index
                    else len(l_int_total_n_features_of_previously_written_matrices)
                ),
            }
        )
        Multiprocessing(
            df_input,
            __MTX_10X_Combine__renumber_barcode_or_feature_index_mtx_10x__,
            int_num_threads,
            global_arguments=[path_folder_mtx_10x_output, flag_renumber_feature_index],
        )
        #         os.remove( f'{path_folder_mtx_10x_output}dict_id_entry_to_index_entry.pickle' ) # remove pickle file

        """ combine parts and add the MTX file header to compose a combined mtx file """
        df_file = GLOB_Retrive_Strings_in_Wildcards(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
        )
        df_file.wildcard_0 = df_file.wildcard_0.astype(int)
        df_file.sort_values("wildcard_0", inplace=True)
        # if 'flag_split_mtx' is True, does not delete the split mtx files
        OS_FILE_Combine_Files_in_order(
            df_file.path.values,
            f"{path_folder_mtx_10x_output}matrix.mtx.gz",
            delete_input_files=not flag_split_mtx,
            header=f"%%MatrixMarket matrix coordinate integer general\n%\n{len( l_t_entry ) if flag_renumber_feature_index else len( set_feature )} {len( set_barcode ) if flag_renumber_feature_index else len( l_t_entry )} {int_total_n_records}\n",
        )  # combine the output mtx files in the given order
        # write a flag indicating that the current output directory contains split mtx files
        with open(f"{path_folder_mtx_10x_output}matrix.mtx.gz.split.flag", "w") as file:
            file.write("completed")

    else:
        """normal operation mode perfoming count merging operations"""
        l_df_mtx, l_df_feature = [], []
        for path_folder_mtx_10x in l_path_folder_mtx_10x_input:
            df_mtx, df_feature = MTX_10X_Read(path_folder_mtx_10x)
            l_df_mtx.append(df_mtx), l_df_feature.append(df_feature)

        # combine mtx
        df_mtx = pd.concat(l_df_mtx)
        df_mtx = df_mtx.groupby(["barcode", "id_feature"]).sum()
        df_mtx.reset_index(drop=False, inplace=True)

        # combine features
        df_feature = pd.concat(l_df_feature)
        df_feature.drop_duplicates(inplace=True)

        MTX_10X_Write(df_mtx, df_feature, path_folder_mtx_10x_output)

        # split a matrix file into multiple files
        MTX_10X_Split(
            path_folder_mtx_10x_output,
            int_max_num_entries_for_chunk=int_max_num_entries_for_chunk,
        )


def MTX_10X_Split(
    path_folder_mtx_10x_output,
    int_max_num_entries_for_chunk=10000000,
    flag_split_mtx=True,
    flag_split_mtx_again=False,
):
    """# 2022-04-28 01:16:35
    split input mtx file into multiple files and write a flag file indicating the splitting has been completed.
    return the list of split mtx files

    'flag_split_mtx' : if 'flag_split_mtx' is True, split input mtx file into multiple files. if False, does not split the input matrix, and just return the list containing a single path pointing to the input matrix. This flag exists for the compatibility with single-thread operations
    'flag_split_mtx_again' : split the input matrix again even if it has beem already split. It will remove previously split files.
    """
    # 'flag_split_mtx' : if False, does not split the input matrix, and just return the list containing a single path pointing to the input matrix
    if not flag_split_mtx:
        return [f"{path_folder_mtx_10x_output}matrix.mtx.gz"]

    """ if 'flag_split_mtx_again' flag is True, remove previously split files """
    path_file_flag = f"{path_folder_mtx_10x_output}matrix.mtx.gz.split.flag"
    if flag_split_mtx_again:
        os.remove(path_file_flag)  # remove the flag
        # remove previously split files
        for path_file in glob.glob(f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"):
            os.remove(path_file)

    """ split input matrix file """
    if not os.path.exists(path_file_flag):  # check whether the flag exists
        index_mtx_10x = 0
        newfile = gzip.open(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", "wb"
        )
        l_path_file_mtx_10x = [
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz"
        ]
        int_num_entries_written_for_the_current_chunk = 0
        with gzip.open(f"{path_folder_mtx_10x_output}matrix.mtx.gz", "rb") as file:
            while True:
                line = file.readline()  # binary string
                if len(line) == 0:
                    newfile.close()  # close the output file
                    break
                """ write the line to the current chunk and update the number of entries written for the current chunk """
                newfile.write(line)
                int_num_entries_written_for_the_current_chunk += 1
                """ initialize the next chunk if a sufficient number of entries were written """
                if (
                    int_num_entries_written_for_the_current_chunk
                    >= int_max_num_entries_for_chunk
                ):
                    newfile.close()  # close the output file
                    index_mtx_10x += 1
                    int_num_entries_written_for_the_current_chunk = 0
                    newfile = gzip.open(
                        f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz",
                        "wb",
                    )
                    l_path_file_mtx_10x.append(
                        f"{path_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz"
                    )
        with open(path_file_flag, "w") as file:
            file.write("completed")
    else:
        """retrieve the list of split mtx files"""
        df = GLOB_Retrive_Strings_in_Wildcards(
            f"{path_folder_mtx_10x_output}matrix.mtx.gz.*.gz"
        )
        df.wildcard_0 = df.wildcard_0.astype(int)
        df.sort_values("wildcard_0", ascending=True, inplace=True)
        l_path_file_mtx_10x = df.path.values
    return l_path_file_mtx_10x
