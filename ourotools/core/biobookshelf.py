from typing import Union, List, Literal, Dict
import numpy as np
import pandas as pd
import os
from io import StringIO  # for converting a string to a file-like stream
import time
import uuid
import glob
import datetime
import pickle
import matplotlib.pyplot as plt
from uuid import uuid4
import gzip
import subprocess
import shutil
import csv
import intervaltree
from multiprocessing import (
    Pool,
    get_context,
    set_start_method,
)  # for multiple processing  # with get_context("spawn").Pool() as pool:
import multiprocessing
import multiprocessing as mp
import collections
from copy import copy, deepcopy

# set default parameters
int_max_num_batches_in_a_queue_for_each_worker = 2 # 2 batches distributed to each process should be optimal, while preventing pipe buffer overloading. However, when deadlocks frequently occurs, changing this argument to 1 will prevent the occurrence of deadlocks.

def Wide(int_percent_html_code_cell_width=95):
    """
    # 20210224
    widen jupyter notebook cell
    """
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:100% !important; }</style>"))
    display(
        HTML(
            "<style>div.cell{width:"
            + str(int_percent_html_code_cell_width)
            + "%;margin-left:"
            + str(100 - int_percent_html_code_cell_width)
            + "%;margin-right:auto;}</style>"
        )
    )
    
def Search_list_of_strings(
    list_of_strings,
    query="cancer",
    return_mask_matched=False,
    return_location_matched=False,
):
    """search list of strings to find strings that contains query string and return the result as a list. if 'return_mask_matched' is True,
    return list_mask for locations of matched entries (return np.array( search_result, dtype = object ), list_mask_matched)
    """
    search_result, list_mask_matched = list(), list()
    list_of_strings = (
        list_of_strings.values
        if type(list_of_strings) is pd.Series
        else list_of_strings
    )  # if type of list_of_strings is pd.Series, even though iterating through pandas series is just fine, to be safe and fast, convert pd.Series to a numpy array
    for (
        string
    ) in (
        list_of_strings
    ):  # search list of strings to find strings that contains query string and return the result as a list
        if not isinstance(string, (float, int)) and query in string:
            search_result.append(string)
            list_mask_matched.append(True)
        else:
            list_mask_matched.append(False)
    if return_mask_matched:
        return np.array(search_result, dtype=object), np.array(
            list_mask_matched, dtype=bool
        )
    elif return_location_matched:
        return np.where(np.array(list_mask_matched, dtype=bool))[0]
    else:
        return search_result

def Search_list_of_strings_with_multiple_query(
    l_str,
    *l_query,
    flag_ignore_case: bool = True,
    return_mask=False,
    return_position=False,
):
    """# 2023-04-28 01:04:38
    Search list of strings with multiple query. for negative query, add '-' in front of the query

    flag_ignore_case : bool = True # ignore cases by default
    """

    arr_mask_matched = np.ones(len(l_str), dtype=bool)
    l_str_input = l_str  # save the reference to the original object
    if flag_ignore_case:
        l_str = deepcopy(l_str)  # copy l_str 'before' converting to lower cases
        for i in range(len(l_str)):
            l_str[i] = l_str[i].lower()
    for query in l_query:
        bool_query_positive = True
        if query[0] == "-":
            bool_query_positive, query = False, query[1:]
        if flag_ignore_case:
            query = query.lower()  # convert the query to lower case
        l_mask_matched_for_a_query = (
            list(True if query in entry else False for entry in l_str)
            if bool_query_positive
            else list(False if query in entry else True for entry in l_str)
        )
        arr_mask_matched = arr_mask_matched & np.array(
            l_mask_matched_for_a_query, dtype="bool"
        )
    if return_position:
        return np.where(arr_mask_matched)[0]
    return (
        arr_mask_matched
        if return_mask
        else np.array(l_str_input, dtype=object)[arr_mask_matched]
    )  # return a subset of the list of input strings

def Search_list_of_strings_Return_mask(data, query, is_negative_query=False):
    if is_negative_query:
        return np.array(
            list(False if query in entry else True for entry in data), dtype=bool
        )
    else:
        return np.array(
            list(True if query in entry else False for entry in data), dtype=bool
        )

def MPL_SAVE(fig_name, l_format=[".pdf", ".png"], close_fig=True, **dict_save_fig):
    """With the given 'fig_name', save fiqures in both svg and png format
    'l_format' : list of image extensions for saving files
    """
    for str_format in l_format:
        MATPLOTLIB_savefig(
            fig_name, format=str_format, close_fig=False, **dict_save_fig
        )
    if close_fig:
        plt.close()

def MATPLOTLIB_savefig(title, dpi=200, folder=None, close_fig=True, format=".png"):
    if "." not in format:
        format = "." + format
    plt.savefig(
        folder + To_window_path_compatible_str(title) + format,
        dpi=200,
        bbox_inches="tight",
    )  # prevent x or y labels from cutting off
    if close_fig:
        plt.close()

def PICKLE_Write(path_file, data_object):
    """write binary pickle file of a given data_object"""
    with open(path_file, "wb") as handle:
        pickle.dump(data_object, handle, protocol=pickle.HIGHEST_PROTOCOL)


def PICKLE_Read(path_file):
    """write binary pickle file of a given data_object"""
    with open(path_file, "rb") as handle:
        data_object = pickle.load(handle)
    return data_object


def TIME_GET_timestamp(flag_human_readable=False):
    """Get timestamp of current time in "%Y%m%d_%H%M" format"""
    cur_time = datetime.datetime.now()  # retrieve current time
    return (
        cur_time.strftime("%Y/%m/%d %H:%M")
        if flag_human_readable
        else cur_time.strftime("%Y%m%d_%H%M")
    )


def To_window_path_compatible_str(a_string):
    """
        replace following characters to '_' so that a given string will be compatible for Window file system :
    : (colon)    " (double quote)    / (forward slash)    \ (backslash)    | (vertical bar or pipe)    ? (question mark)    * (asterisk)
        Also, replace new line character into '_'
    """
    return (
        a_string.replace("\n", "_")
        .replace(":", "_")
        .replace('"', "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("|", "_")
        .replace("?", "_")
        .replace("*", "_")
    )

def MATPLOTLIB_savefig(title, dpi=200, folder=None, close_fig=True, format=".png"):
    if "." not in format:
        format = "." + format
    plt.savefig(
        folder + To_window_path_compatible_str(title) + format,
        dpi=200,
        bbox_inches="tight",
    )  # prevent x or y labels from cutting off
    if close_fig:
        plt.close()


def MATPLOTLIB_basic_configuration(
    font_size=None,
    font_size_axes_title=None,
    font_size_axes_label=None,
    font_size_xtick=None,
    font_size_ytick=None,
    font_size_legend=None,
    font_size_figure_title=None,
    x_label=None,
    y_label=None,
    title=None,
    x_scale=None,
    y_scale=None,
    show_grid=True,
    show_legend=False,
    savefig=False,
    y_lim=None,
    x_lim=None,
    save_file_name=None,
    folder=None,
    format=".png",
    show_colorbar=False,
):
    """A basic function for confiquring a matplotlib plot"""
    # set font sizes
    if font_size is not None:
        plt.rc("font", size=20)  # controls default text sizes
    if font_size_axes_title is not None:
        plt.rc("axes", titlesize=20)  # fontsize of the axes title
    if font_size_axes_label is not None:
        plt.rc("axes", labelsize=20)  # fontsize of the x and y labels
    if font_size_xtick is not None:
        plt.rc("xtick", labelsize=20)  # fontsize of the x tick labels
    if font_size_ytick is not None:
        plt.rc("ytick", labelsize=20)  # fontsize of the y tick labels
    if font_size_legend is not None:
        plt.rc("legend", fontsize=20)  # legend fontsize
    if font_size_figure_title is not None:
        plt.rc("figure", titlesize=50)  # fontsize of the figure title
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if title is not None:
        plt.title(title)
    if x_scale is not None:
        plt.xscale(x_scale)
    if y_scale is not None:
        plt.yscale(y_scale)
    if x_lim is not None:
        if isinstance(x_lim, (tuple, list)):
            plt.xlim(left=x_lim[0], right=x_lim[1])
        elif isinstance(x_lim, dict):
            plt.xlim(**x_lim)
    if y_lim is not None:
        if isinstance(y_lim, (tuple, list)):
            plt.ylim(bottom=y_lim[0], top=y_lim[1])
        elif isinstance(y_lim, dict):
            plt.ylim(**y_lim)
    plt.grid(show_grid)
    if show_legend:
        plt.legend()
    if savefig:
        if save_file_name is None:  # if 'save_file_name' is not given
            if title is None:
                title = (
                    "Unnamed Plot_" + TIME_GET_timestamp()
                )  # if title is not given, put a default title to save a plot
            MATPLOTLIB_savefig(title=title, folder=folder, format=format)
        else:
            MATPLOTLIB_savefig(title=save_file_name, folder=folder, format=format)
    if show_colorbar:
        plt.colorbar()

MPL_basic_configuration = MATPLOTLIB_basic_configuration


def TYPE_Convert_NP_Array(data, dtype=None):
    """Default dtype = Float"""
    if dtype is None:
        dtype = float
    if type(data) is pd.Series:
        data = data.values.astype(dtype)
    elif type(data) is list:
        data = np.array(data, dtype=dtype)
    elif type(data) is set:
        data = np.array(list(data), dtype=dtype)
    elif type(data) is not np.ndarray:
        print("ERROR: Invalid data type")
        return -1
    return data

def MPL_1D_Sort_Plot(
    data,
    figsize=(5, 3.5),
    annotate_xy_pos_first_column_label=(0.05, 1.09),
    color_alpha=0.5,
    color_threshold=0,
    line_stype=".-",
    x_label="Sorted Entries",
    title="",
    savefig=False,
    color_above="g",
    color_below="r",
    color_percentile_alpha=0.5,
    color_percentile_thres=None,
    color_percentile_lower="b",
    color_percentile_upper="orange",
    thres_n_points=10000,
    **dict_mpl_basic_configure,
):
    """(1) Convert iterable data like series or list into np.ndarray using 'TYPE_Convert_NP_Array' (2) Sort, (3) Visualize on a plot using green and red colors
    to visualize deviation from the given threshold. if color_percentile_thres is not None, annotate upper and lower percentiles with the color given by color_percentile arguments
    if 'data' is pandas.DataFrame with two columns, first sort values in the second column by the first column, visualize, and annotate unique_entries of the first column on the plot. The NaN values in the first column will be ignored.
    """
    bool_flag_sort_using_two_columns = (
        isinstance(data, pd.DataFrame) and len(data.columns.values) == 2
    )  # 'sort_using_two_columns' if a DataFrame with two columns are given as a 'data'.
    if isinstance(
        data, (pd.DataFrame, pd.Series)
    ):  # set default title and y_label by using pandas.Series name if not given
        data_name = data.name if isinstance(data, pd.Series) else data.columns.values[1]
        if data_name:  # if data_name is not None
            if "y_label" not in dict_mpl_basic_configure:
                dict_mpl_basic_configure["y_label"] = (
                    data.name if isinstance(data, pd.Series) else data.columns.values[1]
                )
            if len(title) == 0:
                title = data_name  # set data_name as default title for the plot
    if (
        not bool_flag_sort_using_two_columns
    ):  # convert data as an numpy array and sort the data if 'series__sort_by_index_first' is set to False or given 'data' is not pandas.Series
        data = TYPE_Convert_NP_Array(data, dtype=float)
        if type(data) is not np.ndarray:
            return -1
    if bool_flag_sort_using_two_columns:
        data = data.dropna()  # remove np.nan values if present
    else:
        arr_mask_isnan = np.isnan(data)
        if arr_mask_isnan.any():
            data = data[~arr_mask_isnan]
    int_ratio_shrinkage = (
        int(len(data) / thres_n_points) if len(data) > thres_n_points else 1
    )  # shrink the number of data values before sorting for efficient plotting. # limit the number of points beling plotting below 'thres_n_points'
    if int_ratio_shrinkage > 1:
        print(
            "Since len( data ) = {} + ({} NaN values) > thres_n_points = {}, Perform random sampling of the data (one value for every {} values) for efficient visualization".format(
                len(data), arr_mask_isnan.sum(), thres_n_points, int_ratio_shrinkage
            )
        )
        data = (
            data.iloc[::int_ratio_shrinkage]
            if bool_flag_sort_using_two_columns
            else data[::int_ratio_shrinkage]
        )
    data = (
        data.sort_values(list(data.columns), inplace=False, ignore_index=True)
        if bool_flag_sort_using_two_columns
        else np.sort(data)
    )  # sort data
    arr_data = (
        data.iloc[:, 1].values if bool_flag_sort_using_two_columns else data
    )  # plot data
    x_range = np.arange(len(arr_data))
    x_axis = np.full_like(x_range, color_threshold)
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.plot(arr_data, line_stype, color="k")
    if color_percentile_thres is not None:  # fill colors in the plot
        index_thres = int(float(len(arr_data)) * color_percentile_thres / 100)
        x_axis_plus_max = np.full_like(x_range, arr_data[-1])
        ax.fill_between(
            x_range[:index_thres],
            x_axis_plus_max[:index_thres],
            x_axis[:index_thres],
            facecolor=color_percentile_lower,
            interpolate=True,
            alpha=color_percentile_alpha,
        )
        ax.fill_between(
            x_range[-index_thres:],
            x_axis_plus_max[-index_thres:],
            x_axis[-index_thres:],
            facecolor=color_percentile_upper,
            interpolate=True,
            alpha=color_percentile_alpha,
        )
    ax.fill_between(
        x_range,
        arr_data,
        x_axis,
        where=arr_data >= x_axis,
        facecolor=color_above,
        interpolate=True,
        alpha=color_alpha,
    )
    ax.fill_between(
        x_range,
        arr_data,
        x_axis,
        where=arr_data <= x_axis,
        facecolor=color_below,
        interpolate=True,
        alpha=color_alpha,
    )
    plt.sca(ax)  # set x_ticks properly after shrinkage
    arr_xticks = plt.xticks()[0][1:-1]
    plt.xticks(arr_xticks, (arr_xticks * int_ratio_shrinkage).astype(int))
    if (
        bool_flag_sort_using_two_columns
    ):  # annotate unique entries in the first columns by which values of second columns were first sorted.
        l_unqiue_entry = sorted(data.iloc[:, 0].unique())
        dict_unique_entry_to_int_representation = dict(
            (unique_entry, int_representation)
            for int_representation, unique_entry in enumerate(l_unqiue_entry)
        )
        dict_int_representation_to_unique_entry = dict(
            (int_representation, unique_entry)
            for int_representation, unique_entry in enumerate(l_unqiue_entry)
        )
        data.iloc[:, 0] = list(
            dict_unique_entry_to_int_representation[entry]
            for entry in data.iloc[:, 0].values
        )
        l_start_of_unique_entry = (
            [0] + list(np.where(np.diff(data.iloc[:, 0].values))[0]) + [len(data)]
        )
        for int_representation, unique_entry in enumerate(l_unqiue_entry):
            x_pos = (
                l_start_of_unique_entry[int_representation]
                + l_start_of_unique_entry[int_representation + 1]
            ) / 2
            ax.annotate(
                unique_entry,
                xy=(x_pos, 1.02),
                xycoords=("data", "axes fraction"),
                ha="center",
            )
        ax.annotate(
            data.columns.values[0],
            xy=annotate_xy_pos_first_column_label,
            xycoords=("axes fraction", "axes fraction"),
            ha="center",
        )
    MATPLOTLIB_basic_configuration(
        x_label=x_label,
        title=TIME_GET_timestamp() + "\n" + title,
        savefig=savefig,
        **dict_mpl_basic_configure,
    )



class Map(object):
    def __init__(self, dict_a2b):
        self.dict_a2b = dict_a2b

    def a2b(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return np.nan

    def a2b_if_mapping_available_else_Map_a2a(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return a


def DICTIONARY_Build_from_arr(arr, order_index_entry=True, index_start=0):
    if order_index_entry:
        return dict(
            (index, entry)
            for entry, index in zip(arr, np.arange(index_start, len(arr) + index_start))
        )
    else:
        return dict(
            (entry, index)
            for entry, index in zip(arr, np.arange(index_start, len(arr) + index_start))
        )

def DICTIONARY_Find_Max(dict_value):  # 2020-07-29 23:55:58
    """find key value with the maximum value in a given dictionary, and return 'key_max', 'value_max'"""
    if len(dict_value) == 0:
        return None, None  # if an empty dictionary is given, return None
    key_max = next(
        dict_value.__iter__()
    )  # retrieve the a key and value from the dictionary as an initial key-value pair with 'max' value
    value_max = dict_value[key_max]
    for key in dict_value:
        value = dict_value[key]
        if value_max < value:
            key_max = key
            value_max = value
    return key_max, value_max

def UUID():
    """return a 128bit universal unique identifier"""
    return uuid4().hex


def OS_Run(
    l_args,
    path_file_stdout=None,
    path_file_stderr=None,
    return_output=True,
    remove_default_output_files=True,
    stdout_binary=False,
):
    """# 2021-03-30 19:41:16
    Run a process and save stdout and stderr as a file.

    'return_output' : return the output as dictionary of strings
    'remove_default_output_files' : remove default stdout and stderr files containing the output of the process when 'path_file_stdout' and 'path_file_stderr' were not given.
    'stdout_binary' : set this flag to True if stdout is binary.
    """
    from subprocess import call as subprocess_call

    uuid_process = UUID()  # set uuid of the process
    # define default stdout and stdin files and set approproate flags
    flag_path_file_stdout_was_given = path_file_stdout is not None
    flag_path_file_stderr_was_given = path_file_stderr is not None

    # default stdout/stderr files will be written to the current working directory
    path_cwd = os.getcwd()
    if not flag_path_file_stdout_was_given:
        path_file_stdout = f"{path_cwd}/{uuid_process}.out.txt"
    if not flag_path_file_stderr_was_given:
        path_file_stderr = f"{path_cwd}/{uuid_process}.err.txt"

    with open(
        path_file_stdout, "w+b" if stdout_binary else "w+"
    ) as fout:  # excute and read std output and std errors of a process
        with open(path_file_stderr, "w+") as ferr:
            out = subprocess_call(l_args, stdout=fout, stderr=ferr)
            fout.seek(0)
            stdout = fout.read()
            ferr.seek(0)
            stderr = ferr.read()
    # remove default output files
    if not flag_path_file_stdout_was_given:
        os.remove(path_file_stdout)
    if not flag_path_file_stderr_was_given:
        os.remove(path_file_stderr)
    return {"stdout": stdout, "stderr": stderr} if return_output else None


def LIST_Split(
    l=None,
    n_split=0,
    return_slice=False,
    flag_contiguous_chunk=False,
    arr_weight_for_load_balancing=None,
    return_split_arr_weight=False,
):
    """# 2022-05-26 10:14:31
    split a list into 'n_split' number of chunks. if 'return_slice' is True, return slice() instances instead of actually spliting the given list-like object.
    performs load balancing based on given list of weights (the 'arr_weight_for_load_balancing' argument)

    'flag_contiguous_chunk' : split the list in a continguous manner so that each chunk contains a region of a original list
    'arr_weight_for_load_balancing' : only valid when 'flag_contiguous_chunk' is True. 'arr_weight_for_load_balancing' should contains the list of weights for each element for load balancing
    'return_split_arr_weight' : return split arr_weights, too (only valid if 'flag_contiguous_chunk' == True and valid arr_weight is given through the 'arr_weight_for_load_balancing' element)
    """
    # retrieve slice
    if flag_contiguous_chunk:
        if (
            arr_weight_for_load_balancing is None
        ):  # process equal number of entries for each chunk
            int_num_entries_for_each_chunk = int(np.ceil(len(l) / n_split))
            l_slice = list(
                slice(
                    index_split * int_num_entries_for_each_chunk,
                    (index_split + 1) * int_num_entries_for_each_chunk,
                )
                for index_split in np.arange(n_split)
            )
        else:  # if an array of weights are given, use the weights to balance the load for each chunk
            # convert dtype of the array to increase the resolution and prevent error due to small resolution of np.float32 # 2022-05-26 10:07:10 by ahs2202
            arr_weight_for_load_balancing = np.array(
                arr_weight_for_load_balancing, dtype=np.float64
            )
            # calculate total weaight for each chunk
            int_total_weight = np.sum(arr_weight_for_load_balancing)
            int_total_weight_for_each_chunk = np.ceil(int_total_weight / n_split)

            # collect the start positions of each chunk
            index_chunk = 0
            l_index_start_of_chunk = [0]
            for index, accumulated_weight in enumerate(
                np.cumsum(arr_weight_for_load_balancing)
            ):
                if (
                    int_total_weight_for_each_chunk * (index_chunk + 1)
                    < accumulated_weight
                ):  # if the accumulated bytes is larger than the 'int_total_weight_for_each_chunk' times the number of chunk written, record a chunk boundary.
                    l_index_start_of_chunk.append(
                        index
                    )  # mark the current position as the start of the chunk (and thus the start of the next chunk)
                    index_chunk += 1  # update the index of the chunk
            if (
                len(l_index_start_of_chunk) > n_split
            ):  # when a possible overflow/errors from too low resolution was detected, correct the boundary
                l_index_start_of_chunk[n_split] = len(arr_weight_for_load_balancing)
                l_pos_start_chunk = l_index_start_of_chunk[: n_split + 1]
            else:
                l_pos_start_chunk = l_index_start_of_chunk + [
                    len(arr_weight_for_load_balancing)
                ]
            l_slice = list(
                slice(
                    l_pos_start_chunk[index_split], l_pos_start_chunk[index_split + 1]
                )
                for index_split in np.arange(n_split)
            )
    else:
        l_slice = list(
            slice(index_split, None, n_split) for index_split in np.arange(n_split)
        )
    if return_slice:
        return l_slice  # if 'return_slice' is True, return slice() instances instead of actually spliting the given list-like object
    else:
        if (
            flag_contiguous_chunk
            and arr_weight_for_load_balancing is not None
            and return_split_arr_weight
        ):  # return split input list and the weights
            return list(l[a_slice] for a_slice in l_slice), list(
                arr_weight_for_load_balancing[a_slice] for a_slice in l_slice
            )
        else:
            return list(l[a_slice] for a_slice in l_slice)


def PD_Threshold(df, AND_operation=True, verbose = False, **dict_thresholds):
    """
    Select rows of a given DataFrame or indices of Series based on a given threshold for each given column or the given series.
    Add 'b' or 'B' at the end of column_label to select rows below the threshold, or add 'a' or 'A' to select rows above the threshold.
    If 'AND_operation' is true, filter generated from the given threshold will be combined with AND operation before filtering rows of a given dataframe
    2024-11-21 20:00 by IEUM An
    """
    # detect data type
    flag_is_dataframe = type(df) is pd.DataFrame
    if flag_is_dataframe :
        set_columns = set( df.columns.values )
        flag_has_multiindex = isinstance( df.columns, pd.MultiIndex ) # detect dataframe containing 'multiindex' columns
        if flag_has_multiindex :
            n_multiindex_levels = len( df.columns.names )
            def compose_multiindex_col( dict_select, l_multiindex = [ ] ) :
                """
                compose multiindex col recursively
                2024-11-21 by IEUM An
                """
                if len( l_multiindex ) + 1 == n_multiindex_levels :
                    return dict( ( tuple( l_multiindex + [ k ] ), dict_select[ k ] ) for k in dict_select )
                else :
                    dict_select_multiindex_col = dict( )
                    for k in dict_select :
                        v = dict_select[ k ]
                        dict_select_multiindex_col.update( compose_multiindex_col( v, l_multiindex = l_multiindex + [ k ] ) )
                    return dict_select_multiindex_col
            dict_thresholds = compose_multiindex_col( dict_thresholds, l_multiindex = [ ] ) # convert 'dict_thresholds' (nested dictionaries) to tuple-based dictionary for the select operation
    else :
        assert( type(df) is pd.Series ) # other allowed data type is pandas Series
        set_columns = set([""]) # allow search values for pandas Series without the need to refer to a column
        
    mask_filter = (
        np.ones(len(df), dtype=bool) if AND_operation else np.zeros(len(df), dtype=bool)
    )
    for col_direction, threshold in dict_thresholds.items():
        if flag_has_multiindex :
            # when input data is a multiindex dataframe
            # extract the filtering direction from the last entry of the input tuple.
            last_col_direction = col_direction[ -1 ]
            col, direction = last_col_direction[:-1], last_col_direction[-1] 
            col = tuple( list( col_direction[ : -1 ] ) + [ col ] )
        else :
            col, direction = col_direction[:-1], col_direction[-1]
        if col not in set_columns:
            if verbose :
                print("'{}' column_label does not exist in the given DataFrame".format(col))
            continue
        data = df[col].values if type(df) is pd.DataFrame else df.values
        if direction.lower() == "a":
            current_mask = data > threshold
        elif direction.lower() == "b":
            current_mask = data < threshold
        else:
            if verbose :
                print(
                    "'{}' direction is not either 'a' or 'b' and thus invalid".format(
                        direction
                    )
                )
            continue
        mask_filter = (
            current_mask & mask_filter if AND_operation else current_mask | mask_filter
        )
    return df[mask_filter]
    
def PD_Select(df, deselect = False, verbose = False, ** dict_select ) :
    """
    Select and filter rows of df according to the given dict_select. If 'deselect' is set to True, deselect rows according to the given dict_select  Usage example : PANDAS_Select( df_meta_imid_ubi, dict(  Data_Type = [ 'Proteome', 'Ubi_Profiling' ], Value_Type = 'log2fc' ) )
    2024-11-21 19:50 by IEUM An
    """
    # detect data type
    flag_is_dataframe = type(df) is pd.DataFrame
    if flag_is_dataframe :
        set_columns = set( df.columns.values )
        flag_has_multiindex = isinstance( df.columns, pd.MultiIndex ) # detect dataframe containing 'multiindex' columns
        if flag_has_multiindex :
            n_multiindex_levels = len( df.columns.names )
            def compose_multiindex_col( dict_select, l_multiindex = [ ] ) :
                """
                compose multiindex col recursively
                2024-11-21 by IEUM An
                """
                if len( l_multiindex ) + 1 == n_multiindex_levels :
                    return dict( ( tuple( l_multiindex + [ k ] ), dict_select[ k ] ) for k in dict_select )
                else :
                    dict_select_multiindex_col = dict( )
                    for k in dict_select :
                        v = dict_select[ k ]
                        dict_select_multiindex_col.update( compose_multiindex_col( v, l_multiindex = l_multiindex + [ k ] ) )
                    return dict_select_multiindex_col
            dict_select = compose_multiindex_col( dict_select, l_multiindex = [ ] ) # convert 'dict_select' (nested dictionaries) to tuple-based dictionary for the select operation
    else :
        assert( type(df) is pd.Series ) # other allowed data type is pandas Series
        
    # detect query type
    for col, query in dict_select.items():
        if not flag_is_dataframe :
            data_values = (
                df.index.values if col == "index" else df.values
            )  # select values or indices of a given pd.Series
        else :
            if col not in set_columns and col != "index":
                if verbose :
                    print("'{}' does not exist in columns of a given DataFrame".format(col))
                continue
            data_values = df.index.values if col == "index" else df[col].values
        if isinstance(
            query, (list, tuple, np.ndarray, set)
        ):  # if data to be selected is iterable
            query = (
                set(query) if isinstance(query, (list, tuple, np.ndarray)) else query
            )  # convert query into set
            df = (
                df[
                    list(
                        False if data_value in query else True
                        for data_value in data_values
                    )
                ]
                if deselect
                else df[
                    list(
                        True if data_value in query else False
                        for data_value in data_values
                    )
                ]
            )
        else:
            df = df[data_values != query] if deselect else df[data_values == query]
    return df


def DF_Deduplicate_without_loading_in_memory(
    path_file_dataframe: str,
    path_file_dataframe_deduplicated: str,
    l_col_for_identifying_duplicates: Union[int, str, list[str]],
    flag_header_is_present: bool = True,
    str_delimiter: str = "\t",
    flag_collect_the_number_of_processed_lines: bool = False,
):
    """
    # 2023-01-06 22:43:43
    (Assumes the given dataframe is gzipped.)
    similar to pandas.DataFrame.drop_duplicates, except that the dataframe will not be loaded into the memory. duplicates are identified and redundant records will be dropped with keep = 'first' setting


    inputs:
    'path_file_dataframe' directory of the dataframe to remove duplicates
    'path_file_dataframe_deduplicated' : an output file directory containing unique records
    'l_col_for_identifying_duplicates' : list of column names (should be all string types) if 'flag_header_is_present' is True else list of column indices (should be all integer types). a single column name/index can be also given.

    returns: { 'int_num_lines' : int_num_lines, 'dict_t_val_count' : dict_t_val_count }
    int_num_lines: collect the number of processed lines
    dict_t_val_count: the counts of each unique combination of values.
    """
    flag_deduplicate_using_a_single_column = isinstance(
        l_col_for_identifying_duplicates, (str, int)
    )  # retrieve a flag indicating the de-duplication will be performed using a single column
    if flag_deduplicate_using_a_single_column:
        l_col_for_identifying_duplicates = [
            l_col_for_identifying_duplicates
        ]  # wrap a single column name/index in a list

    """ open an output file """
    newfile = gzip.open(path_file_dataframe_deduplicated, "wb")
    with gzip.open(path_file_dataframe, "rb") as file:
        """retrieve list of indices of columns for identifying redundant records"""
        if flag_header_is_present:
            # read header
            line = file.readline()
            newfile.write(line)  # write header
            l_col = line.decode().strip().split(str_delimiter)  # parse header
            l_int_index_col_for_identifying_duplicates = list(
                l_col.index(col) for col in l_col_for_identifying_duplicates
            )  # retrieve indices of the columns that should be accessed
        else:
            # if header is not present, 'l_col_for_identifying_duplicates' will be used as list of indices of columns for extracting values by which a unique record will be identified
            l_int_index_col_for_identifying_duplicates = (
                l_col_for_identifying_duplicates
            )

        """ check whether each record is redundant and write only unique record to the output file """
        int_num_lines = 0
        dict_t_val_count = (
            dict()
        )  # a dict that will record the number of occurrences of a unique set of values
        while True:
            line = file.readline()  # read line
            if len(line) == 0:
                break
            int_num_lines += 1
            l_val = line.decode().strip().split(str_delimiter)  # parse the current line
            t_val = (
                l_val[l_int_index_col_for_identifying_duplicates[0]]
                if flag_deduplicate_using_a_single_column
                else tuple(
                    l_val[int_index]
                    for int_index in l_int_index_col_for_identifying_duplicates
                )
            )  # retrieve a single value / a tuple of values containined in the line
            if t_val not in dict_t_val_count:
                newfile.write(
                    line
                )  # write line if the line contain a unqiue set of values
                dict_t_val_count[t_val] = 0  # initialize the count
            dict_t_val_count[t_val] += 1  # increase the count
    newfile.close()
    return {"int_num_lines": int_num_lines, "dict_t_val_count": dict_t_val_count}


def PD_Binary_Flag_Select(
    df, name_col_binary_flag, int_bit_flag_position, flag_select=True
):
    """
    'flag_select' : if True, select rows containing the binary flags at the given position 'int_bit_flag_position'
    """
    # handle empty dataframe
    if len(df) == 0:
        return df
    int_bit_flag = 1 << int_bit_flag_position
    return df[
        list(
            (f & int_bit_flag) > 0 if flag_select else (f & int_bit_flag) == 0
            for f in df[name_col_binary_flag].values
        )
    ]


def check_gtf_or_gff3_format( path_file_input : str ) -> bool :
    """
    based on the input file path, check GTF or GFF3 format of the input annotation file

    return:
    flag_gtf_format # a boolean flag indicating whether the file is in GTF format (in which case the flag will be True) or GFF3 format (in which case the flag will be False)
    
    2024-12-03 17:10 by IEUM An
    """
    path_file_input = os.path.abspath( path_file_input ) # retrieve absolute path
    
    name = path_file_input.rsplit( '/', 1 )[ -1 ] # retrieve file name
    # if the file contains gzipped content according to the file name (file extension), drop '.gz' from the file name.
    if name.rsplit( '.', 1 )[ -1 ].lower( ) in { 'gz' } :
        name = name.rsplit( '.', 1 )[ 0 ]
    
    # check whether the file contains GTF annotation according to the file name (file extension)
    flag_gtf_format = name.rsplit( '.', 1 )[ -1 ].lower( ) in { 'gtf' }

    return flag_gtf_format


def GTF_Parse_Attribute(attr):
    """
    parse attribute string of a gtf file
    # 2021-02-06 18:51:47
    2024-12-03 by IEUM An, support for integer/float data type was added 
    """
    dict_data = dict()
    l_e = attr.split('; ')
    if len( l_e ) > 0 :
        last_e = l_e[ -1 ]
        if len( last_e ) > 0 and last_e[ -1 ] == ';' :
            l_e[ -1 ] = last_e[ : -1 ] # discard the ';' character at the end of the last element
    for e in l_e :
        e = e.strip()
        # skip empty element
        if len(e) == 0:
            continue
        # check dtype
        if '"' in e :
            # string dtype
            str_key, str_value = e.split(' "', 1) # there should be 1 occurrence of ' "'
            # remove the trailing '"'
            if str_value[-1] == '"':
                str_value = str_value[:-1]
            value = str_value # use the value as-is
        else :
            # integer or float dtype
            str_key, str_value = e.split(' ', 1) # the attribute name should not contain ' ' when integer/float value is contained
            value = float( str_value ) # for simplicity, it will be converted to float first
            if value == int( value ) : # perform the equivalent test, and convert the value to the integer
                value = int( value )
        dict_data[str_key] = value
    return dict_data


def GFF3_Parse_Attribute(attr):
    """
    # 2021-08-16 16:35:15
    Parse attribute string of GFF3 formated gene annotation file
    """
    return dict(e.split("=", 1) for e in attr.split(";") if "=" in e)


def GTF_Read(
    path_gtf,
    parse_attr=True,
    remove_chr_from_seqname=True,
    flag_verbose=False,
    flag_gtf_format=True, # deprecated
    flag_gtf_gzipped=False, # deprecated
    **dict_filter_gtf,
):
    """
    Load gzipped or plain text GTF files into pandas DataFrame. the file's gzipped-status can be explicitly given by 'flag_gtf_gzipped' argument.
    'path_gtf' : directory to the gtf file or a dataframe containing GTF records to parse attributes
    'parse_attr' : parse gtf attribute if set to True
    'dict_filter_gtf' : keyworded arguments for 'PD_Select', which will be used to filter df_gtf before parsing attributes
    remove_chr_from_seqname=True, # by default, remove 'chr' prefix from the sequence names.

    # 2022-02-08 08:55:37
    # 2024-12-03 21:57 by IEUM An 
    """
    # read GTF/GFF annotation and sort the entries
    try:
        df = (
            pd.read_csv(
                path_gtf,
                sep="\t",
                header=None,
                low_memory=False,
                comment="#",
                skip_blank_lines=True,
            )
            if isinstance(path_gtf, (str))
            else path_gtf
        )  # if 'path_gtf' is a string, read the given gtf file from disk using the given directory # ignore comments
        df.columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
    except:
        # return empty GTF when an error occurs during reading a GTF file
        if flag_verbose:
            print(
                "error reading GTF file. Might be an empty GTF file, returning empty dataframe"
            )
        df = pd.DataFrame(
            [[""] * 9],
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
            ],
        )
        df = df.iloc[:0]
        return df
    df = df.sort_values(["seqname", "start"]).reset_index(drop=True)

    # detect file format
    flag_gtf_format = check_gtf_or_gff3_format( path_gtf )

    # remove 'chr' prefix 
    if remove_chr_from_seqname:
        df["seqname"] = list(
            seqname if seqname[:3] != "chr" else seqname[3:]
            for seqname in df.seqname.values
        )

    # filter entries using GTF/GFF columns (excluding attribute column) 
    if len(dict_filter_gtf) > 0:
        df = PD_Select(df, **dict_filter_gtf)
        df.reset_index(drop=True, inplace=True)
    
    # parse attributes
    if parse_attr:
        return df.join(
            pd.DataFrame(
                list(GTF_Parse_Attribute(attr) for attr in df.attribute.values)
                if flag_gtf_format
                else list(GFF3_Parse_Attribute(attr) for attr in df.attribute.values)
            )
        )

    # return the result
    return df


def GTF_Write(
    df_gtf, path_file, flag_update_attribute=True, flag_filetype_is_gff3=False
):
    """
    write gtf file as an unzipped tsv file
    'flag_update_attribute' : ignore the 'attribute' column present in the given dataframe 'df_gtf', and compose a new column based on all the non-essential columns of the dataframe.
    'flag_filetype_is_gff3' : a flag indicating the filetype of the output file. According to the output filetype, columns containing attributes will be encoded into the values of the attribute column before writing the file.

    # 2021-08-24 21:02:08 - 2024-11-25 14:30 by IEUM-An
    """
    if flag_update_attribute:
        l_col_essential = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

        l_col_attribute = list(
            col for col in df_gtf.columns.values if col not in l_col_essential
        )  # all non-essential columns will be considered as the columns

        empty_attribute_value = '' # define empty attribute value
        if len( l_col_attribute ) == 0 :
            # when there is no columns to include in the attribute
            l_attribute_new = list( empty_attribute_value for _ in range( len( df_gtf ) ) )
        else :
            # compose the attribute column using the given columns 
            l_attribute_new = list()
            for arr in df_gtf[l_col_attribute].values:
                str_attribute = empty_attribute_value # initialize as 'empty_attribute_value'
                for name, val in zip(l_col_attribute, arr):
                    if isinstance(val, float) and np.isnan(val):
                        continue
                    str_attribute += (
                        f"{name}={val};" if flag_filetype_is_gff3 else f'{name} "{val}"; '
                    )  # encode attributes according to the gff3 file format
                str_attribute = str_attribute.strip()
                l_attribute_new.append(str_attribute)
        # add the attribute column
        df_gtf["attribute"] = l_attribute_new  # update attributes
    df_gtf[
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
        ]
    ].to_csv(path_file, index=False, header=None, sep="\t", quoting=csv.QUOTE_NONE)

def GTF_Interval_Tree(
    path_file_gtf, feature=["gene"], value="gene_name", drop_duplicated_intervals=False
):
    """# 2022-05-20 22:48:35
    Return an interval tree containing intervals retrieved from the given gtf file.

    'path_file_gtf' : directory to the gtf file or dataframe iteslf
    'feature' : list of features in the gtf to retrieve intervals from the gtf file
    'value' : list of columne names (including GTF attributes) to include as a list of values for each interval in the returned interval tree.
    'drop_duplicated_intervals' : drop duplicated intervals
    """
    # read GTF file
    if isinstance(
        path_file_gtf, (str)
    ):  # if 'path_file_gtf' is a string object, used the value to read a GTF file.
        df_gtf = GTF_Read(path_file_gtf, parse_attr=True)
    else:  # assumes 'path_file_gtf' is a dataframe containing GTF records if it is not a string object
        df_gtf = path_file_gtf
    # retrieve gtf records of given list of features if valid query is given
    if feature is not None:
        df_gtf = PD_Select(df_gtf, feature=feature)
    if len(df_gtf) == 0:  # return an empty dictionary if df_gtf is an empty dataframe
        return dict()
    df_gtf.dropna(
        subset=[value] if isinstance(value, str) else value, inplace=True
    )  # value should be valid
    # remove duplicated intervals
    if drop_duplicated_intervals:
        df_gtf.drop_duplicates(subset=["seqname", "start", "end"], inplace=True)
    dict_it = dict()
    for arr_interval, arr_value in zip(
        df_gtf[["seqname", "start", "end"]].values, df_gtf[value].values
    ):
        seqname, start, end = arr_interval
        if seqname not in dict_it:
            dict_it[seqname] = intervaltree.IntervalTree()
        # add the interval with a given list of value
        dict_it[seqname].addi(
            start - 1,
            end,
            tuple(arr_value) if isinstance(arr_value, np.ndarray) else arr_value,
        )  # 1-based coordinate system to 1-based coordinate
    return dict_it
    

def FASTA_Read(
    path_file_fasta,
    print_message=False,
    remove_space_in_header=False,
    return_dataframe=False,
    parse_uniprot_header=False,
    header_split_at_space=False,
):  # 2020-08-21 14:38:09
    """# 2021-12-28 20:19:44
    for a file-like object of file of 'path_file_fasta' directory, parse the content into a dictionary. 'remove_space_in_header' option remove space in the header line (required for ABySS output)
    'parse_uniprot_header': if set to True and 'return_dataframe' is set to True, parse uniprot sequence header into [ 'accession', 'description', 'uniprot_source', 'uniprot_acc', 'uniprot_name' ]
    'header_split_at_space' : split a header at the first space, and use the string before the first space as a header string
    """
    dict_header_to_seq = dict()
    dict_duplicated_header_count = dict()
    bool_flag_input_gzipped = False
    if hasattr(path_file_fasta, "readline"):
        file = path_file_fasta  # if file-like object was given instead of path_file_fasta, use path_file_fasta as file
    else:
        bool_flag_input_gzipped = path_file_fasta[-3:] == ".gz"
        file = (
            gzip.open(path_file_fasta, "rb")
            if bool_flag_input_gzipped
            else open(path_file_fasta, "r")
        )  # open file of path_file_fasta depending on the detected gzipped status
    line = (file.readline().decode() if bool_flag_input_gzipped else file.readline())[
        :-1
    ]
    while True:
        str_header = line
        l_seq = list()
        while True:
            line = (
                file.readline().decode() if bool_flag_input_gzipped else file.readline()
            )[:-1]
            if not line or line[:1] == ">":
                break
            l_seq.append(line)
        header, seq = str_header[1:], "".join(l_seq)
        if remove_space_in_header:
            header = header.replace(" ", "_")
        elif (
            header_split_at_space
        ):  # if 'header_split_at_space' is true : split a header at the first space, and use the string before the first space as a header string
            header = header.split(" ", 1)[0]
        if header in dict_duplicated_header_count:
            dict_duplicated_header_count[
                header
            ] += 1  # handle sequences with duplicated headers
        else:
            dict_duplicated_header_count[header] = 1
        if dict_duplicated_header_count[header] > 1:
            dict_header_to_seq[
                header
                + "_dup_{index}".format(index=dict_duplicated_header_count[header])
            ] = "".join(
                l_seq
            )  # if current fastq header in already exists, add '_dup_{index}' to the header to make it unique
        else:
            dict_header_to_seq[header] = "".join(l_seq)
        if not line:
            break
    file.close()  # close the file
    if len(dict_header_to_seq) == 1 and tuple(dict_header_to_seq) == tuple([""]):
        dict_header_to_seq = (
            dict()
        )  # if 'dict_header_to_seq' is empty, return an empty dictionary as 'dict_header_to_seq'
    if print_message:
        print(pd.Series(dict_duplicated_header_count))
    if return_dataframe:  # return parsed fasta file as a dataframe
        df = pd.Series(dict_header_to_seq).reset_index()
        df.columns = ["header", "seq"]
        df["length"] = df.seq.apply(len)
        if (
            parse_uniprot_header
        ):  # if set to True and 'return_dataframe' is set to True, parse uniprot sequence header into [ 'accession', 'description', 'uniprot_source', 'uniprot_acc', 'uniprot_name' ]
            l_l_value = list()
            for header in df.header.values:
                accession, description = header.split(" ", 1)
                uniprot_source, uniprot_acc, uniprot_name = accession.split("|")
                l_l_value.append(
                    [accession, description, uniprot_source, uniprot_acc, uniprot_name]
                )
            df = df.join(
                pd.DataFrame(
                    l_l_value,
                    columns=[
                        "accession",
                        "description",
                        "uniprot_source",
                        "uniprot_acc",
                        "uniprot_name",
                    ],
                )
            )
        return df
    else:
        return dict_header_to_seq


def FASTA_Write(
    path_file_fasta,
    dict_fasta=None,
    l_id=None,
    l_seq=None,
    overwrite_existing_file=False,
    int_num_characters_for_each_line=60,
):
    """# 2022-01-12 17:16:43 (improved performance when writing very large fasta file, such as human genome sequence)
    'int_num_characters_for_each_line' : insert the newline character for every 'int_num_characters_for_each_line' number of characters.
    write fasta file at the given directory with dict_fastq (key = fasta_header, value = seq) or given list of id (fasta_header) and seq
    write gzipped fasta file if 'path_file_fasta' ends with '.gz'
    """

    """ if 'dict_fasta' was given, compose 'dict_fasta' from the lists of 'l_id' and 'l_seq'. """
    if dict_fasta is None:
        dict_fasta = dict()
        if isinstance(l_id, str):  # when a single pair of sequence and id were given
            dict_fasta[l_id] = l_seq
        else:  # when a lists of sequences and ids were given
            for str_id, str_seq in zip(l_id, l_seq):
                dict_fasta[str_id] = str_seq

    """ open file """
    # detect gzipped status of the output file
    flag_file_gzipped = path_file_fasta.rsplit(".", 1)[1] == "gz"
    if os.path.exists(path_file_fasta) and not overwrite_existing_file:
        print("the file already exists")
        return -1
    with gzip.open(path_file_fasta, "wb") if flag_file_gzipped else open(
        path_file_fasta, "w"
    ) as newfile:
        for str_id in dict_fasta:
            seq = dict_fasta[str_id]
            len_seq = len(seq)

            """ write a header line """
            line_header = ">" + str_id + "\n"
            newfile.write(line_header.encode() if flag_file_gzipped else line_header)

            pos_seq = 0
            while True:
                if (len_seq - pos_seq) > int_num_characters_for_each_line:
                    line = (
                        seq[pos_seq : pos_seq + int_num_characters_for_each_line] + "\n"
                    )
                    newfile.write(line.encode() if flag_file_gzipped else line)
                    pos_seq += int_num_characters_for_each_line
                else:
                    """if less than 'int_num_characters_for_each_line' number of character remains, write all remaining characters and exit the loop"""
                    line = seq[pos_seq:len_seq] + "\n"
                    newfile.write(line.encode() if flag_file_gzipped else line)
                    pos_seq += len_seq - pos_seq
                    break


def DF_Build_Index_Using_Dictionary(
    df, l_col_for_index, str_delimiter=None, function_transform=None
):  # 2020-08-06 17:12:59
    """# 2021-09-07 19:34:16
    return a dictionary with key = index or multi-index (when list of columns is givne through the 'l_col_for_index' argument) for a given 'l_col_for_index' and value = list of integer index for df.values (not df.index.values)
    Using Python dictionary and numpy array can be upto ~2000 times faster than using pandas.DataFrame.loc[ ] method for accessing multi-indexed rows.

    'l_col_for_index' : (1) a string that is equal to the name of the column to be indexed or (2) a list of names of columns for composing a multi-index.
    'str_delimiter' : when given, values in the SINGLE column referred to by the 'l_col_for_index' argument (should be a string and not a list) will be split using the character(s) of 'str_delimiter', and indexed separately. When the value of a row contains multiple unique entries as determined by 'str_delimiter', the row will be indexed multiple times with different keys
    'function_transform' : a function for transforming values for building the dictionary-based index a using single column with non-None 'str_delimiter' value (default: None)
    """
    dict_index = dict()
    if isinstance(
        l_col_for_index, str
    ):  # when only single col_name was given for index
        if (
            str_delimiter is None
        ):  # when each value of the column will be used to index a row
            for int_index, index in enumerate(df[l_col_for_index].values):
                index = (
                    index if function_transform is None else function_transform(index)
                )
                if index in dict_index:
                    dict_index[index].append(int_index)
                else:
                    dict_index[index] = [int_index]
        elif isinstance(
            str_delimiter, str
        ):  # when rows will be indexed by unique entries contained in each value (as determined by 'str_delimiter' characters) of the given column
            for int_index, val in enumerate(df[l_col_for_index].values):
                flag_val_is_nan = isinstance(
                    val, float
                )  # check whether the current value is np.nan value
                if flag_val_is_nan or str_delimiter not in val:
                    index = (
                        val
                        if flag_val_is_nan or function_transform is None
                        else function_transform(val)
                    )
                    if index in dict_index:
                        dict_index[index].append(int_index)
                    else:
                        dict_index[index] = [int_index]
                else:
                    for index in val.split(str_delimiter):
                        index = (
                            index
                            if function_transform is None
                            else function_transform(index)
                        )
                        if index in dict_index:
                            dict_index[index].append(int_index)
                        else:
                            dict_index[index] = [int_index]

    else:  # when a list of col_names was given for index
        for int_index, arr_index in enumerate(df[l_col_for_index].values):
            t_index = tuple(arr_index)
            if t_index in dict_index:
                dict_index[t_index].append(int_index)
            else:
                dict_index[t_index] = [int_index]
    return dict_index


def GTF_Build_Mask(
    dict_seqname_to_len_seq,
    df_gtf=None,
    str_feature=None,
    remove_chr_from_seqname=True,
    path_folder_output=None,
):
    """# 2021-10-10 01:27:23
    build a bitarray mask of entries in the gtf file (0 = none, 1 = at least one entry exists).
    if 'path_folder_output' is given, save the generated mask to the given folder (an empty directory is recommended)
    if 'path_folder_output' is given but 'df_gtf' is not given, load the previously generated mask from the 'path_folder_output'

    'dict_seqname_to_len_seq' : dictionary containing sequence length information of the genome
    'df_gtf' : directory to gtf file or a dataframe containing gtf. If none is given, load previously generated masks from 'path_folder_output'
    'path_folder_output' : directory to save or load masks
    """
    from bitarray import bitarray

    # retrieve the absolute path of the 'path_folder_output'
    if path_folder_output is not None:
        path_folder_output = os.path.abspath(path_folder_output)
        path_folder_output += "/"
    if df_gtf is None:
        if path_folder_output is None:
            print("required inputs were not given")
            return -1
        else:
            """load mask from the 'path_folder_output'"""
            dict_seqname_to_ba = dict()
            for seqname, path_file in GLOB_Retrive_Strings_in_Wildcards(
                f"{path_folder_output}*.bin"
            ).values:
                ba = bitarray()
                with open(path_file, "rb") as file:
                    ba.fromfile(file)
                dict_seqname_to_ba[seqname] = ba[
                    : dict_seqname_to_len_seq[seqname]
                ]  # drop additional '0' added to the end of the binary array
            return dict_seqname_to_ba
    # create output folder if it does not exist
    if not os.path.exists(path_folder_output):
        os.makedirs(path_folder_output, exist_ok=True)
    # remove 'chr' characters from seqnames in the 'dict_seqname_to_len_seq'
    if remove_chr_from_seqname:
        dict_seqname_to_len_seq = dict(
            (
                seqname if seqname[:3] != "chr" else seqname[3:],
                dict_seqname_to_len_seq[seqname],
            )
            for seqname in dict_seqname_to_len_seq
        )

    """ initialize bitarrys to zeros (the length bitarrays are same as the reference sequences) """
    dict_seqname_to_ba = dict()
    for seqname in dict_seqname_to_len_seq:
        ba = bitarray(dict_seqname_to_len_seq[seqname])
        ba.setall(0)
        dict_seqname_to_ba[seqname] = ba

    """ read gtf and build mask """
    if isinstance(df_gtf, str):
        df_gtf = GTF_Read(df_gtf)
    # handle an empty GTF file
    if len(df_gtf) == 0:
        """save empty masks as files"""
        if path_folder_output is not None:
            for seqname in dict_seqname_to_ba:
                with open(f"{path_folder_output}{seqname}.bin", "wb") as file:
                    dict_seqname_to_ba[seqname].tofile(file)
        return dict_seqname_to_ba
    if (
        str_feature is not None
    ):  # select only specific features form the gtf if query is given
        df_gtf = PD_Select(df_gtf, feature=str_feature)
    if remove_chr_from_seqname:
        df_gtf["seqname"] = list(
            seqname if seqname[:3] != "chr" else seqname[3:]
            for seqname in df_gtf.seqname.values
        )
    df_gtf.start -= 1  # 1-based coordinate -> 0-based coordinate
    for seqname, start, end in df_gtf[
        ["seqname", "start", "end"]
    ].values:  # 0-based coordinate
        # only consider sequences in 'dict_seqname_to_len_seq'
        if seqname not in dict_seqname_to_ba:
            continue
        dict_seqname_to_ba[seqname][start:end] = 1
    """ save masks as files """
    if path_folder_output is not None:
        for seqname in dict_seqname_to_ba:
            with open(f"{path_folder_output}{seqname}.bin", "wb") as file:
                dict_seqname_to_ba[seqname].tofile(file)
    return dict_seqname_to_ba


# functions related to the command line interface
def Parse_Printed_Table(str_output):
    """# 2022-08-08 11:39:01
    Parse printed table by identifying the position of columns and inferring datatypes using pandas module
    """
    l_line = str_output.split("\n")
    # survey the number of ' ' space characters in each line
    int_max_characters_in_each_line = max(
        len(line) for line in l_line
    )  # retrieve the maximum number of characters a line contains
    arr_space_counter = np.zeros(int_max_characters_in_each_line, dtype=int)
    for line in l_line:
        # add padding containing ' ' to the end of the sentence for accurate detection of columns
        if len(line) < int_max_characters_in_each_line:
            line += " " * (int_max_characters_in_each_line - len(line))
        for i, c in enumerate(line):
            if c == " ":
                arr_space_counter[i] += 1

    arr_pos_col_space = np.where(arr_space_counter == max(arr_space_counter))[
        0
    ]  # retrieve positions of the columns of ' ' space characters

    arr_pos_of_columns_marking_the_boundary = (
        [-1, arr_pos_col_space[0]]
        + list(arr_pos_col_space[1:][np.diff(arr_pos_col_space) != 1])
        + [int_max_characters_in_each_line + 1]
    )  # retrieve positions of the columns marking the boundaries

    # collect values
    l_l = []
    for line in l_line:
        l = []
        for i in range(len(arr_pos_of_columns_marking_the_boundary) - 1):
            col_before, col_after = arr_pos_of_columns_marking_the_boundary[i : i + 2]
            l.append(line[col_before + 1 : col_after].strip())
        l_l.append("\t".join(l))
    df = pd.read_csv(StringIO("\n".join(l_l)), sep="\t")
    return df


def PIP_List_Packages():
    """# 2022-08-08 11:32:04
    list installed packages
    """
    return (
        Parse_Printed_Table(os.popen("pip list").read().strip())
        .drop(index=[0])
        .set_index("Package")
    )


# slice functions
def Slice_to_Range(sl, length):
    """# 2022-06-28 21:47:51
    iterate indices from the given slice
    """
    assert isinstance(sl, slice)  # make sure sl is slice
    # convert slice to integer indices
    for i in range(*sl.indices(length)):
        yield i


def COUNTER(l_values, dict_counter=None, ignore_float=True):  # 2020-07-29 23:49:51
    """
    Count values in l_values and return a dictionary containing count values. if 'dict_counter' is given, countinue counting by using the 'dict_counter'. if 'ignore_float' is True, ignore float values, including np.nan

    l_values # list of values OR another 'dict_counter'

    # 2024-11-18 17:00 by IEUM An
    """
    if dict_counter is None:
        dict_counter = dict()
    if isinstance( l_values, dict ) :
        # when 'l_values' is another 'dict_counter'
        dict_counter_incoming = l_values
        # update 'dict_counter' with 'dict_counter_incoming'
        for k in dict_counter_incoming :
            if k in dict_counter :
                dict_counter[ k ] += dict_counter_incoming[ k ]
            else :
                dict_counter[ k ] = dict_counter_incoming[ k ]
    else :
        # when 'l_values' is list of values
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


def LIST_COUNT(
    iterable,
    return_series=True,
    duplicate_filter=2,
    dropna=True,
    sort_series_by_values=True,
    convert_tuple_to_string=False,
):
    """
    # 20210224
    return a dictionary where key = each element in a given list, value = counts of the element in the list. if 'duplicate_filter' is not None, return entries that are duplicated 'duplicate_filter' times or more.
    """
    if dropna and isinstance(iterable, pd.Series):
        iterable = (
            iterable.dropna()
        )  # if dropna is set to 'True', dropn NaN values before counting
    if isinstance(next(iterable.__iter__()), (np.ndarray, list)):
        iterable = list(
            map(tuple, iterable)
        )  # if value is non-hashable list of numpy array, convert a value to a hashable format, tuple
    dict_counted = COUNTER(iterable)
    if (
        convert_tuple_to_string
    ):  # if 'convert_tuple_to_string' is True and values in a given list are tuples, convert tuples into string
        dict_counted__tuple_converted_to_string = dict()
        for key in dict_counted:
            value = dict_counted[key]
            if isinstance(key, (tuple)):
                dict_counted__tuple_converted_to_string[
                    ("{}, " * len(key))[:-2].format(*key)
                ] = value  # convert tuple into string concatanated with ', '
            else:
                dict_counted__tuple_converted_to_string[key] = value
        dict_counted = dict_counted__tuple_converted_to_string
    if return_series:
        s_counted = pd.Series(dict_counted)
        if duplicate_filter is not None:
            s_counted = s_counted[s_counted >= duplicate_filter]
        if sort_series_by_values:
            s_counted = s_counted.sort_values(ascending=False)
        return s_counted
    else:
        return dict_counted


def DICTIONARY_Find_keys_with_max_value(dict_value):
    """# 2021-11-24 20:44:07
    find a list of key values with the maximum value in a given dictionary, and return 'l_key_max', 'value_max'
    """
    value_max = None  # initialize max value
    l_key_max = []  # list of key with max_values
    if len(dict_value) != 0:  # if the dictionary is not empty
        for key in dict_value:
            value = dict_value[key]
            if value_max is None:
                value_max = value
                l_key_max.append(key)
            elif value_max > value:
                continue
            elif value_max < value:
                l_key_max = [key]
                value_max = value
            elif (
                value_max == value
            ):  # if the another key contains the current max value, add the key to the list of keys with max values
                l_key_max.append(key)
    return l_key_max, value_max


def DICTIONARY_Build_from_arr(arr, order_index_entry=True, index_start=0):
    if order_index_entry:
        return dict(
            (index, entry)
            for entry, index in zip(arr, np.arange(index_start, len(arr) + index_start))
        )
    else:
        return dict(
            (entry, index)
            for entry, index in zip(arr, np.arange(index_start, len(arr) + index_start))
        )


def GLOB_Retrive_Strings_in_Wildcards(
    str_glob,
    l_path_match=None,
    return_dataframe=True,
    retrieve_file_size=False,
    retrieve_last_modified_time=False,
    time_offset_in_seconds=3600 * 9,
):  # 2020-11-16 18:20:52
    """# 2022-01-09 23:25:48
    retrieve strings in '*' wildcards in list of matched directories for the given string containing '*' wildcards. return strings in wildcards as a nested lists. Consecutive wildcards should not be used ('**' should not be used in the given string)
    'retrieve_file_size': if 'return_dataframe' is True, return file sizes in bytes by using os.stat( path_match ).st_size
    'retrieve_last_modified_time': return the last modified time with pandas datetime datatype
    'time_offset_in_seconds': offset in seconds to Coordinated Universal Time (UTC)"""
    l_path_match = (
        glob.glob(str_glob) if l_path_match is None else l_path_match
    )  # retrive matched directories using glob.glob if 'l_path_match' is not given
    l_intervening_str = str_glob.split(
        "*"
    )  # retrive intervening strings in a glob string containing '*' wildcards
    l_l_str_in_wildcard = list()
    for (
        path_match
    ) in l_path_match:  # retrive strings in wildcards for each matched directory
        path_match_subset = path_match.split(l_intervening_str[0], 1)[1]
        l_str_in_wildcard = list()
        for intervening_str in l_intervening_str[1:]:
            if len(intervening_str) > 0:
                str_in_wildcard, path_match_subset = path_match_subset.split(
                    intervening_str, 1
                )
            else:
                str_in_wildcard, path_match_subset = (
                    path_match_subset,
                    "",
                )  # for the wildcard at the end of the given string, put remaining string into 'str_in_wildcard' and empties 'path_match_subset'
            l_str_in_wildcard.append(str_in_wildcard)
        l_l_str_in_wildcard.append(l_str_in_wildcard)
    if (
        return_dataframe
    ):  # return dataframe containing strings in wildcards and matched directory
        df = pd.DataFrame(
            l_l_str_in_wildcard,
            columns=list(
                "wildcard_" + str(index) for index in range(str_glob.count("*"))
            ),
        )
        df["path"] = l_path_match
        if retrieve_file_size:
            df["size_in_bytes"] = list(
                os.stat(path_match).st_size for path_match in l_path_match
            )
            df["size_in_gigabytes"] = df["size_in_bytes"] / 2**30
        if retrieve_last_modified_time:
            df["time_last_modified"] = list(
                datetime.datetime.utcfromtimestamp(
                    os.path.getmtime(path_file) + time_offset_in_seconds
                ).strftime("%Y-%m-%d %H:%M:%S")
                for path_file in df.path.values
            )
            df.time_last_modified = pd.to_datetime(
                df.time_last_modified
            )  # convert to datetime datatype
        return df
    else:
        return l_l_str_in_wildcard


def Multiprocessing_Batch_Generator_and_Workers(
    gen_batch,
    process_batch,
    post_process_batch=None,
    int_num_threads : int =15,
    int_num_seconds_to_wait_before_identifying_completed_processes_for_a_loop : float = 0.2,
    flag_wait_for_a_response_from_worker_after_sending_termination_signal : bool = True, # wait until all worker exists before resuming works in the main process
):
    """
    'Multiprocessing_Batch_Generator_and_Workers' : multiprocessing using batch generator and workers.
    all worker process will be started using the default ('fork' in UNIX) method.
    perform batch-based multiprocessing using the three components, (1) gen_batch, (2) process_batch, (3) post_process_batch. (3) will be run in the main process, while (1) and (2) will be offloaded to worker processes.
    the 'batch' and result returned by 'process_batch' will be communicated to batch processing workers through pipes.

    'gen_batch' : a generator object returning batches
    'process_batch( pipe_receiver, pipe_sender )' : a function that can process batch from 'pipe_receiver'. should terminate itself when None is received. 'pipe_sender' argument is to deliver the result to the main process, and should be used at the end of code to notify the main process that the work has been completed. sending 'None' through 'pipe_sender' will terminate the block and the main process will be unblocked (however, the works will be continued to be distributed and performed by the child processes).
    'post_process_batch( result )' : a function that can process return value from 'process_batch' function in the main process. operations that are not thread/process-safe can be done here, as these works will be serialized in the main thread.
    'int_num_threads' : the number of threads(actually processes) including the main process. For example, when 'int_num_threads' is 3, 2 worker processes will be used. one thread is reserved for batch generation.
    'int_num_seconds_to_wait_before_identifying_completed_processes_for_a_loop' : number of seconds to wait for each loop before checking which running processes has been completed
    flag_wait_for_a_response_from_worker_after_sending_termination_signal : bool = True, # wait until all worker exists before resuming works in the main process.
    # 2024-10-28 00:32
    """

    def __batch_generating_worker(
        gen_batch,
        l_pipe_sender_input,
        l_pipe_receiver_output,
        pipe_sender_output_to_main_process,
    ):
        """# 2023-10-06 01:52:24 
        define a worker for generating batch and distributing batches across the workers, receives results across the workers, and send result back to the main process
        """
        q_batch = collections.deque()  # initialize queue of batchs
        int_num_batch_processing_workers = len(l_pipe_sender_input)
        flag_batch_generation_completed = False  # flag indicating whether generating batchs for the current input sam file was completed
        arr_num_batch_being_processed = np.zeros(
            int_num_batch_processing_workers, dtype=int
        )  # retrieve the number of batches currently being processed in each worker. if this number becomes 0, assumes the worker is available
        while True:
            """retrieve batch (if available)"""
            if not flag_batch_generation_completed and len( q_batch ) < int_num_batch_processing_workers * int_max_num_batches_in_a_queue_for_each_worker : # if batch generation has not been completed, or the number of batches that have been generated are not large, continue to generate batches.
                try:
                    batch = next(gen_batch)  # retrieve the next barcode
                    q_batch.appendleft(batch)  # append batch
                except StopIteration:
                    flag_batch_generation_completed = True
            else:
                # if batch generation has been completed, all batches have been distributed, and processed, exit the loop
                if flag_batch_generation_completed and len(q_batch) == 0 and arr_num_batch_being_processed.sum() == 0:
                    break
                # sleep for a while
                time.sleep(
                    int_num_seconds_to_wait_before_identifying_completed_processes_for_a_loop
                )  # sleep

            """ collect completed works """
            for index_worker in range(int_num_batch_processing_workers):
                while l_pipe_receiver_output[
                    index_worker
                ].poll():  # until results are available
                    pipe_sender_output_to_main_process.send(
                        l_pipe_receiver_output[index_worker].recv()
                    )  # retrieve result, and send the result back to the main process
                    arr_num_batch_being_processed[
                        index_worker
                    ] -= 1  # update the number of batches being processed by the worker

            """ if workers are available and there are remaining works to be distributed, distribute works """
            index_worker = 0  # initialize the index of the worker
            while (
                len(q_batch) > 0
                and (
                    arr_num_batch_being_processed
                    < int_max_num_batches_in_a_queue_for_each_worker
                ).sum()
                > 0
            ):  # if there is remaining batch to be distributed or at least at least one worker should be available
                if (
                    arr_num_batch_being_processed[index_worker]
                    <= arr_num_batch_being_processed.mean()
                    and arr_num_batch_being_processed[index_worker]
                    < int_max_num_batches_in_a_queue_for_each_worker
                ):  # the load of the current worker should be below the threshold # if the load for the current worker is below the pool average, assign the work to the process (load-balancing)
                    l_pipe_sender_input[index_worker].send(q_batch.pop())
                    arr_num_batch_being_processed[index_worker] += 1
                index_worker = (
                    1 + index_worker
                ) % int_num_batch_processing_workers  # retrieve index_worker of the next worker

        # notify batch-processing workers that all workers are completed
        for pipe_s in l_pipe_sender_input:
            pipe_s.send(None)
        # notify the main process that all batches have been processed
        pipe_sender_output_to_main_process.send(None)
        return

    int_num_batch_processing_workers = max(
        1, int_num_threads - 2
    )  # retrieve the number of workers for processing batches # minimum number of worker is 1
    # compose pipes
    l_pipes_input = list(mp.Pipe() for i in range(int_num_batch_processing_workers))
    l_pipes_output = list(mp.Pipe() for i in range(int_num_batch_processing_workers))
    pipe_sender_output_to_main_process, pipe_receiver_output_to_main_process = mp.Pipe()
    # compose workers
    l_batch_processing_workers = list(
        mp.Process(
            target=process_batch, args=(l_pipes_input[i][1], l_pipes_output[i][0])
        )
        for i in range(int_num_batch_processing_workers)
    )  # compose a list of batch processing workers
    p_batch_generating_worker = mp.Process(
        target=__batch_generating_worker,
        args=(
            gen_batch,
            list(s for s, r in l_pipes_input),
            list(r for s, r in l_pipes_output),
            pipe_sender_output_to_main_process,
        ),
    )
    # start workers
    for p in l_batch_processing_workers:
        p.start()
    p_batch_generating_worker.start()

    # post-process batches
    while True:
        res = pipe_receiver_output_to_main_process.recv()
        if res is None:
            break
        if post_process_batch is not None:
            post_process_batch(
                res
            )  # process the result returned by the 'process_batch' function in the 'MAIN PROCESS', serializing potentially not thread/process-safe operations in the main thread.
    
    # if 'flag_wait_for_a_response_from_worker_after_sending_termination_signal' is True, wait until a response is received from the worker
    if flag_wait_for_a_response_from_worker_after_sending_termination_signal :
        for s, r in l_pipes_output : # pipe receiving responses from batch workers
            r.recv( )

    # allow the worker processes to join
    for p in l_batch_processing_workers:
        p.join()
    p_batch_generating_worker.join() 
    return # exit


def Multiprocessing(
    arr,
    Function,
    n_threads=12,
    path_temp="/tmp/",
    Function_PostProcessing=None,
    global_arguments=[],
    col_split_load=None,
):
    """# 2022-02-23 10:55:34
    split a given iterable (array, dataframe) containing inputs for a large number of jobs given by 'arr' into 'n_threads' number of temporary files, and folks 'n_threads' number of processes running a function given by 'Function' by givning a directory of each temporary file as an argument. if arr is DataFrame, the temporary file will be split DataFrame (tsv format) with column names, and if arr is 1d or 2d array, the temporary file will be tsv file without header
    By default, given inputs will be randomly distributed into multiple files. In order to prevent separation of a set of inputs sharing common input variable(s), use 'col_split_load' to group such inputs together.

    'Function_PostProcessing' : if given, Run the function before removing temporary files at the given temp folder. uuid of the current session and directory of the temporary folder are given as arguments to the function.
    'global_arguments' : a sort of environment variables (read only) given to each process as a list of additional arguments in addition to the directory of the input file. should be used to use local variables inside main( ) function if this function is called inside the main( ) function.
                         'global_arguments' will be passed to 'Function_PostProcessing', too.
    'col_split_load' : a name of column or a list of column names (or integer index of column or list of integer indices of columns if 'arr' is not a dataframe) for grouping given inputs when spliting the inputs into 'n_threads' number of dataframes. Each unique tuple in the column(s) will be present in only one of split dataframes.
    'n_threads' : if 'n_threads' is 1, does not use multiprocessing module, but simply run the function with the given input. This behavior is for enabiling using functions running Multiprocessing in another function using Multiprocessing, since multiprocessing.Pool module does not allow nested pooling.
    """
    from multiprocessing import Pool

    if isinstance(
        arr, (list)
    ):  # if a list is given, convert the list into a numpy array
        arr = np.array(arr)
    str_uuid = UUID()  # create a identifier for making temporary files
    l_path_file = (
        list()
    )  # split inputs given by 'arr' into 'n_threads' number of temporary files
    if (
        col_split_load is not None
    ):  # (only valid when 'arr' is dataframe) a name of column for spliting a given dataframe into 'n_threads' number of dataframes. Each unique value in the column will be present in only one split dataframe.
        if isinstance(
            arr, pd.DataFrame
        ):  # if arr is DataFrame, the temporary file will be split DataFrame (tsv format) with column names
            if isinstance(
                col_split_load, (str)
            ):  # if only single column name is given, put it in a list
                col_split_load = [col_split_load]

            # randomly distribute distinct tuples into 'n_threads' number of lists
            dict_index = DF_Build_Index_Using_Dictionary(
                arr, col_split_load
            )  # retrieve index according to the tuple contained by 'col_split_load'
            l_t = list(dict_index)
            n_t = len(l_t)  # retrieve number of tuples
            if (
                n_t < n_threads
            ):  # if the given number of thread is larger than the existing number of tuples, set the number of tuples as the number of threads
                n_threads = n_t
            np.random.shuffle(l_t)
            l_l_t = list(l_t[i::n_threads] for i in range(n_threads))

            arr_df = arr.values
            l_col = arr.columns.values

            for index_chunk in range(n_threads):
                l_t_for_the_chunk = l_l_t[index_chunk]
                # retrieve integer indices of the original array for composing array of the curreht chunk
                l_index = []
                for t in l_t_for_the_chunk:
                    l_index.extend(dict_index[t])
                path_file_temp = (
                    path_temp + str_uuid + "_" + str(index_chunk) + ".tsv.gz"
                )
                pd.DataFrame(arr_df[np.sort(l_index)], columns=l_col).to_csv(
                    path_file_temp, sep="\t", index=False
                )  # split a given dataframe containing inputs with groupping with a given list of 'col_split_load' columns
                l_path_file.append(path_file_temp)
        else:
            print(
                "'col_split_load' option is only valid when the given 'arr' is dataframe, exiting"
            )
            return -1
    else:
        # if number of entries is larger than the number of threads, reduce the n_threads
        if len(arr) < n_threads:
            n_threads = len(arr)
        if isinstance(
            arr, pd.DataFrame
        ):  # if arr is DataFrame, the temporary file will be split DataFrame (tsv format) with column names
            for index_chunk in range(n_threads):
                path_file_temp = (
                    path_temp + str_uuid + "_" + str(index_chunk) + ".tsv.gz"
                )
                arr.iloc[index_chunk::n_threads].to_csv(
                    path_file_temp, sep="\t", index=False
                )
                l_path_file.append(path_file_temp)
        else:  # if arr is 1d or 2d array, the temporary file will be tsv file without header
            l_chunk = LIST_Split(arr, n_threads)
            for index, arr in enumerate(
                l_chunk
            ):  # save temporary files containing inputs
                path_file_temp = path_temp + str_uuid + "_" + str(index) + ".tsv"
                if len(arr.shape) == 1:
                    df = pd.DataFrame(arr.reshape(arr.shape[0], 1))
                elif len(arr.shape) == 2:
                    df = pd.DataFrame(arr)
                else:
                    print("invalid inputs: input array should be 1D or 2D")
                    return -1
                df.to_csv(path_file_temp, sep="\t", header=None, index=False)
                l_path_file.append(path_file_temp)

    if n_threads > 1:
        with Pool(n_threads) as p:
            l = p.starmap(
                Function,
                list([path_file] + list(global_arguments) for path_file in l_path_file),
            )  # use multiple process to run the given function
    else:
        """if n_threads == 1, does not use multiprocessing module"""
        l = [Function(l_path_file[0], *list(global_arguments))]

    if Function_PostProcessing is not None:
        Function_PostProcessing(str_uuid, path_temp, *global_arguments)

    for path_file in glob.glob(path_temp + str_uuid + "*"):
        os.remove(path_file)  # remove temporary files

    return l  # return mapped results


class Offload_Works:
    """# 2023-08-27 18:26:25 
    a class for offloading works in a separate server process without blocking the main process. similar to async. methods, but using processes instead of threads.

    int_max_num_workers : Union[ int, None ] = None # maximum number of worker processes. if the maximum number of works are offloaded, submitting a new work will fail. By default, there will be no limit of the number of worker processes
    """

    def __init__(
        self,
        int_max_num_workers: Union[int, None] = None,
    ):
        """# 2023-01-07 12:13:22"""
        # handle default arguments
        if not isinstance(int_max_num_workers, int) or int_max_num_workers <= 0:
            int_max_num_workers = None  # by default, there will be no limit of the number of worker processes
        self._int_max_num_workers = int_max_num_workers

        # initialize
        self._dict_worker = dict()  # a dictionary that store the list of active workers
        self._dict_res = (
            dict()
        )  # a dictionary that store the completed results of submitted works until it is retrieved.

    @property
    def is_worker_available(self) :
        """ # 2023-08-12 15:42:56 
        return a binary flag indicating whether a worker is available
        """
        return self.int_num_active_workers != self.int_max_num_workers
        
    @property
    def int_num_active_workers(self):
        """# 2023-01-07 12:30:24
        return the number of active worker processes that are actively doing the works
        """
        self.collect_results( flag_return_results = False )  # collect completed results
        return len(self._dict_worker)

    @property
    def int_num_completed_results(self):
        """# 2023-01-07 12:30:31
        return the number of completed results
        """
        self.collect_results( flag_return_results = False )  # collect completed results
        return len(self._dict_res)

    @property
    def int_max_num_workers(self):
        """# 2023-01-07 12:17:17
        'int_max_num_workers' is read only property
        """
        return self._int_max_num_workers

    def collect_results(self, flag_return_results : bool = False):
        """# 2023-08-27 18:18:40 
        collect completed results of the submitted works to the current object and store the results internally.

        flag_return_results : bool = False # If True, return the dictionary containing the all the completed results. The returned results will be flushed from the internal data container. If False, return the number of collected results
        """
        int_num_collected_results = 0  # count the number of collected works
        for str_uuid_work in list(self._dict_worker):
            worker = self._dict_worker[str_uuid_work]
            if worker["pipe_receiver"].poll():
                self._dict_res[str_uuid_work] = worker[ "get_result" ]( ) # collect and store the result
                worker["process"].join()  # dismiss the worker process
                del self._dict_worker[str_uuid_work]  # delete the record of the worker
                int_num_collected_results += 1
            del worker
        if flag_return_results :
            dict_result = self._dict_res # retrieve all completed results
            self._dict_res = dict( ) # initialize the internal container for saving results
            return dict_result   
        else :
            return int_num_collected_results

    def submit_work(self, func, args: tuple = (), kwargs: dict = dict(), associated_data = None ):
        """# 2023-08-12 16:59:57 
        submit a work

        return a 'str_uuid_work' that identify the submitted work. the id will be used to check the status of the work and retrieve the result of the work
        associated_data = None # a data object associated with the work. the data will be returned with the result returned by the function as a dictionary using the following format:
            { 'result' : res, 'associated_data' : associated_data }
            If None is given, simply the result returned by the function will be returned.
        """
        # raise error when already the maximum number of workers are working
        if (
            self.int_max_num_workers is not None
            and self.int_num_active_workers >= self.int_max_num_workers
        ):
            raise RuntimeError(
                f"maximum number of workers are working. an addtional work cannot be received."
            )

        import multiprocessing as mp
        import uuid

        # initialize a worker process
        worker = dict()
        # add 'associated_data'
        if associated_data is not None :
            worker[ 'associated_data' ] = associated_data
        
        # create a pipe for communication
        pipe_sender, pipe_receiver = mp.Pipe()
        worker["pipe_receiver"] = pipe_receiver

        def __worker(pipe_sender):
            res = func(*args, **kwargs)  # run the work
            pipe_sender.send(res)  # send the result

        worker["process"] = mp.Process(target=__worker, args=(pipe_sender,))
        worker["process"].start()  # start the work
        
        # define a function to get result
        def __get_result( ) :
            res = worker[ "pipe_receiver" ].recv()
            return res if associated_data is None else { 'result' : res, 'associated_data' : associated_data }
        worker["get_result"] = __get_result

        str_uuid_work = uuid.uuid4().hex  # assign uuid to the current work
        self._dict_worker[str_uuid_work] = worker
        return str_uuid_work  # return the id of the work

    def check_status(self, str_uuid_work: str):
        """# 2023-01-07 12:19:45
        check the status of a submitted work using the 'str_uuid_work' of the work that has been given when the work was submitted.

        checking the result of unfinished work will not block

        return True if the work has been completed. return False if the work has not been completed.
        """
        self.collect_results( flag_return_results = False )  # collect the completed results
        return str_uuid_work in self._dict_res

    def wait(self, str_uuid_work: str):
        """# 2023-01-07 13:24:14
        wait until the given work has been completed.
        """
        self.collect_results( flag_return_results = False )  # collect the completed results

        # if the work currently active, wait until the work is completed and collect the result
        if str_uuid_work in self._dict_worker:
            worker = self._dict_worker[str_uuid_work]  # retrieve the worker
            self._dict_res[str_uuid_work] =  worker[ "get_result" ]( ) # wait until the result become available and collect and store the result
            worker["process"].join()  # dismiss the worker process
            del self._dict_worker[str_uuid_work]  # delete the record of the worker

    def retrieve_result(self, str_uuid_work: str):
        """# 2023-01-07 12:20:02
        retrieving the result of unfinished work will block until the result become available.
        once a result is retrieved, all the information about the submitted work will be deleted from the object.

        if the work cannot be identified, None will be returned.
        """
        self.wait(str_uuid_work)  # wait until the work is completed

        # if a result is available, return the result and delete from the internal storage
        if str_uuid_work in self._dict_res:
            return self._dict_res.pop(str_uuid_work)

    def wait_all(self, flag_return_results : bool = False):
        """# 2023-08-12 15:54:11 
        wait all works to be completed.
        flag_return_results : bool = False # wait for all submitted works, and return results as a dictionary
        """
        self.collect_results( flag_return_results = False )  # collect the completed results

        for str_uuid_work in list(self._dict_worker):  # for each uncompleted work
            self.wait(str_uuid_work)  # wait until the work is completed
            
        if flag_return_results : # return the results
            dict_result = self._dict_res # retrieve all completed results
            self._dict_res = dict( ) # initialize the internal container for saving results
            return dict_result

def Workers( 
    func, 
    * args,
    int_num_workers_for_Workers : int = 8, 
    ** kwargs,
) :
    """
    func, # a function to run across workers
    * args, # arguments to the function
    int_num_workers_for_Workers : int = 8, # number of workers for the current function ('Workers')
    ** kwargs, # keyworded arguments to the function
    
    For the arguments or keyworded arguments that are unique to individual runs, Python List should be used to list the arguments for each run. If the original argument is Python List, the argument should be wrapped inside another Python List, repeated N times where N is the total number of runs.
    
    For example, the following usage is possible,
    >>> def func( d, e, a = 0, b = 1, c = 3 ) :
    >>>     print( f"{d + e + a + b + c = }" )
    >>> Workers( func, [ 10, 20, 30, 40, 50 ], 30, int_num_workers_for_Workers = 8, a = [ 1, 2, 3, 4, 5, 6, 7, 8 ], b = 10, c = [ 7, 8, 9, 10, 11, 12, 13 ] )
    d + e + a + b + c = 58d + e + a + b + c = 70d + e + a + b + c = 94d + e + a + b + c = 106d + e + a + b + c = 82
    
    # 2024-01-03 17:33:49 
    """
    def _gen_arguments( args_workers, kwargs_workers ) :
        """
        generate arguments with given args and kwargs
        # 2024-01-03 16:23:42 
        """
        ''' survey the number of independent runs '''
        l_int_num_runs = [ ] 
        for arg_workers in args_workers : # survey 'args_workers'
            if isinstance( arg_workers, list ) : 
                l_int_num_runs.append( len( arg_workers ) )
        for k in kwargs_workers : # survey 'kwargs_workers'
            kwarg_workers = kwargs_workers[ k ] # retrieve the argument
            if isinstance( kwarg_workers, list ) : 
                l_int_num_runs.append( len( kwarg_workers ) )
        int_num_runs = np.min( l_int_num_runs ) # determin the number of runs (minimum of the length of arguments)
        
        ''' generate arguments for individual runs '''
        for idx_run in range( int_num_runs ) : # for 'idx_run'
            ''' compose arguments for a run '''
            args_run, kwargs_run = [ ], dict( )
            for arg_workers in args_workers : # compose 'args_run'
                args_run.append( arg_workers[ idx_run ] if isinstance( arg_workers, list ) else arg_workers )
                
            for k in list( kwargs_workers ) : # modify 'kwargs_workers'
                kwarg_workers = kwargs_workers[ k ] # retrieve the argument
                kwargs_run[ k ] = kwarg_workers[ idx_run ] if isinstance( kwarg_workers, list ) else kwarg_workers
            yield args_run, kwargs_run # yield generated arguments and key-worded arguments

    def _run_func( p_i, p_o ) :
        """
        run the given function with given arguments
        # 2024-01-03 16:23:42 
        """
        while True :
            ins = p_i.recv( )
            if ins is None :
                break
            args, kwargs = ins # parse the input

            # run function
            func( * args, ** kwargs, ) 
            p_o.send( 'completed' ) # notify the work has been completed
        p_o.send( 'completed' ) # notify all works have been completed

    Multiprocessing_Batch_Generator_and_Workers(
        gen_batch = _gen_arguments( args, kwargs ), # generate arguments
        process_batch = _run_func,
        int_num_threads = int_num_workers_for_Workers + 2,
        flag_wait_for_a_response_from_worker_after_sending_termination_signal = True,
    )


def OS_FILE_Combine_Files_in_order(
    l_path_file,
    path_newfile,
    overwrite_existing_file=False,
    delete_input_files=False,
    header=None,
    remove_n_lines=0,
    flag_use_header_from_first_file=False,
    flag_bgzip_output=False,
    int_byte_chuck=100000,
):  # 2020-07-20 11:47:29
    """# 2021-10-28 10:58:38
    combine contents of files in l_path_file and write at path_newfile. if header is given, append header (string type with \n at the end) at the front of the file. if 'remove_n_lines' > 0, remove n lines from each files.
    gzipped files are also supported. However, when input files have mixed gzipped status (i.e. first file is gzipped, while second file is a plain text file), it will cause a TypeError.
    Mixed Bgzip and gzipped statues are allowed.

    'flag_use_header_from_first_file' : copy header from the first file to the new file
    'flag_bgzip_output' : a flag indicating whether the output should be written using bgzip module (Biopython's bgzf module). When the output file's file extension is '.bgz', the output file will be block-gzipped file regardless of the 'flag_bgzip_output' flag's status.
    'int_byte_chuck' : the size of the chuck to be read/write when the input file and output file has different datatype (input file = binary file, output file = plain text, or vice versa)
    """
    # check the validity of inputs
    if os.path.exists(path_newfile) and not overwrite_existing_file:
        print("[OS_FILE_Combine_Files_in_order][ERROR] output file already exists")
        return -1
    if len(l_path_file) == 0:
        print("[OS_FILE_Combine_Files_in_order][ERROR] given list of files is empty")
        return -1

    # set flags of input/output filetypes
    str_output_file_extension = path_newfile.rsplit(".", 1)[
        1
    ]  # retrieve the file extension of the output files
    bool_flag_gzipped_output = str_output_file_extension in [
        "gz",
        "bgz",
    ]  # set boolean flag for gzipped output file
    bool_flag_bgzipped_output = (
        str_output_file_extension == "bgz" or flag_bgzip_output
    )  # set boolean flag for block-gzipped output file # if 'flag_bgzip_output' is True, write a block-gzipped file even though the output file extension is '.gz'
    bool_flag_gzipped_input = l_path_file[0].rsplit(".", 1)[1] in [
        "gz",
        "bgz",
    ]  # set boolean flag for gzipped input file

    """ open an output file """
    if bool_flag_gzipped_output:
        if bool_flag_bgzipped_output:
            import Bio.bgzf

            newfile = Bio.bgzf.BgzfWriter(path_newfile, "wb")  # open bgzip file
        else:
            newfile = gzip.open(path_newfile, "wb")  # open simple gzip file
    else:
        newfile = open(path_newfile, "w")  # open plain text file

    """ write a header to the output file """
    if (
        flag_use_header_from_first_file
    ):  # if a flag indicating copying header from the first file to the new file is set, open the first file and read the header line
        path_file = l_path_file[0]
        with (
            gzip.open(path_file, "rb")
            if bool_flag_gzipped_input
            else open(path_file, "r")
        ) as file:
            header = file.readline()
    if header:
        header = (
            header.decode() if not isinstance(header, (str)) else header
        )  # convert header byte string to string if header is not a string type
        newfile.write(
            (header.encode() if bool_flag_gzipped_output else header)
        )  # write a header line to the output file

    """ copy input files to the output file """
    for path_file in l_path_file:
        with (
            gzip.open(path_file, "rb")
            if bool_flag_gzipped_input
            else open(path_file, "r")
        ) as file:
            if remove_n_lines:
                for index in range(remove_n_lines):
                    file.readline()
            if not (
                bool_flag_gzipped_output ^ bool_flag_gzipped_input
            ):  # if output & input are both binary or plain text (should be same datatype), simply copy byte to byte
                shutil.copyfileobj(file, newfile)
            else:  # if input is binary while output is plain text or vice versa
                while True:
                    content = file.read(int_byte_chuck)
                    if (
                        len(content) == 0
                    ):  # if all content has been read, terminate transferring procedure
                        break
                    newfile.write(
                        content.encode()
                        if bool_flag_gzipped_output
                        else content.decode()
                    )
    newfile.close()
    if delete_input_files:
        for path_file in l_path_file:
            os.remove(path_file)


def INTERVAL_Overlap(
    interval_1,
    interval_2,
    flag_sort_to_retrieve_start_and_end=False,
    flag_0_based_coordinate_system=True,
):  # 2020-08-06 20:44:47
    """Fast, basic function for retrieving overlap between two intervals.
    return number of overlapped length between two intervals (each interval is a tuple or list containing start and end position).
    'flag_0_based_coordinate_system': if interval contains float numbers, set 'flag_0_based_coordinate_system' to True.
    'flag_sort_to_retrieve_start_and_end': if interval is always (start, end), set 'flag_sort_to_retrieve_start_and_end' to False to improve performance. (about 200ns faster)
    """
    if flag_sort_to_retrieve_start_and_end:
        start_1, end_1 = sorted(interval_1)
        start_2, end_2 = sorted(interval_2)
    else:
        start_1, end_1 = interval_1
        start_2, end_2 = interval_2
    if not flag_0_based_coordinate_system:
        start_1, start_2 = start_1 - 1, start_2 - 1
    if (end_1 <= start_2) | (end_2 <= start_1):
        return 0
    else:
        return min(end_1, end_2) - max(start_1, start_2)


"""
classes and functions for sharing data across multiple forked processes
"""
from multiprocessing.managers import BaseManager


class ManagerReadOnly(BaseManager):
    pass


class HostedDict:
    """# 2023-01-08 23:00:20
    A class intended for sharing large read-only object using multiprocessing Manager
    Hosted Dictionary (Read-Only)
    """

    # constructor
    def __init__(self, path_file_pickle):
        # read the pickle object
        self._path_file_pickle = path_file_pickle
        with open(path_file_pickle, "rb") as handle:
            self._data = pickle.load(handle)

    def get_keys(self):
        """# 2023-01-08 23:05:40
        return the list of keys
        """
        return list(self._data)

    def getitem(self, key):
        """# 2023-01-08 23:02:10
        get item of a given key
        """
        return self._data[key] if key in self._data else None

    def subset(self, keys):
        """# 2023-01-08 23:00:59
        return a subset of a dictionary for a given list of keys
        """
        dict_subset = dict()
        for key in keys:
            if key in self._data:
                dict_subset[key] = self._data[key]
        return dict_subset


class HostedDictIntervalTree:
    """# 2023-01-08 23:00:20
    A class intended for sharing large read-only object using multiprocessing Manager
    Hosted Dictionary of Interval Trees (Read-Only)
    """

    def __init__(self, path_file_pickle):
        # read the pickle object
        self._path_file_pickle = path_file_pickle
        with open(path_file_pickle, "rb") as handle:
            self._data = pickle.load(handle)

    def search_query(self, seqname: str, query):
        """# 2023-01-08 23:03:44
        perform interval search for a given query
        """
        if seqname not in self._data:
            return []
        return [[st, en, values] for st, en, values in self._data[seqname][query]]  #

    def search_queries(self, seqname: str, queries):
        """# 2023-01-08 23:01:05
        perform interval search for a given queries

        queries # should be iterable
        """
        n_queries = len(queries)
        if seqname not in self._data:
            return [[] for _ in range(n_queries)]
        return [
            [[st, en, values] for st, en, values in self._data[seqname][query]]
            for query in queries
        ]


# register
ManagerReadOnly.register("HostedDict", HostedDict)
ManagerReadOnly.register("HostedDictIntervalTree", HostedDictIntervalTree)


def FASTQ_Iterate(path_file, return_only_at_index=None):
    """# 2020-12-09 22:22:34
    iterate through a given fastq file.
    'return_only_at_index' : return value only at the given index. For example, for when 'return_only_at_index' == 1, return sequence only.
    """
    if return_only_at_index is not None:
        return_only_at_index = (
            return_only_at_index % 4
        )  # 'return_only_at_index' value should be a value between 0 and 3
    bool_flag_file_gzipped = (
        ".gz" in path_file[-3:]
    )  # set a flag indicating whether a file has been gzipped.
    with gzip.open(path_file, "rb") if bool_flag_file_gzipped else open(
        path_file
    ) as file:
        while True:
            record = (
                [file.readline().decode()[:-1] for index in range(4)]
                if bool_flag_file_gzipped
                else [file.readline()[:-1] for index in range(4)]
            )
            if len(record[0]) == 0:
                break
            if return_only_at_index is not None:
                yield record[return_only_at_index]
            else:
                yield record

def FASTQ_Read(
    path_file,
    return_only_at_index=None,
    flag_add_qname=True,
    int_num_reads: Union[int, None] = None,
):  # 2020-08-18 22:31:31
    """# 2021-08-25 07:06:50
    read a given fastq file into list of sequences or a dataframe (gzipped fastq file supported). 'return_only_at_index' is a value between 0 and 3 (0 = readname, 1 = seq, ...)
    'flag_add_qname' : add a column containing qname in the bam file (space-split read name without '@' character at the start of the read name)
    int_num_reads : Union[ int, None ] = None # the number of reads to include in the output, starting from the start. By default, include all reads.
    """
    if return_only_at_index is not None:
        return_only_at_index = (
            return_only_at_index % 4
        )  # 'return_only_at_index' value should be a value between 0 and 3
    bool_flag_file_gzipped = (
        ".gz" in path_file[-3:]
    )  # set a flag indicating whether a file has been gzipped.
    l_seq = list()
    l_l_values = list()
    """ read fastq file """
    file = gzip.open(path_file, "rb") if bool_flag_file_gzipped else open(path_file)
    int_read_count = 0  # count the number of reads parsed.
    while True:
        record = (
            [file.readline().decode()[:-1] for index in range(4)]
            if bool_flag_file_gzipped
            else [file.readline()[:-1] for index in range(4)]
        )
        if len(record[0]) == 0:
            break
        if return_only_at_index is not None:
            l_seq.append(record[return_only_at_index])
        else:
            l_l_values.append([record[0], record[1], record[3]])
        int_read_count += 1  # increase the read count
        if (
            int_num_reads is not None and int_read_count >= int_num_reads
        ):  # if the number of parsed reads equal to the desired number of reads in the output, stop reading the file
            break
    file.close()
    if return_only_at_index is not None:
        return l_seq
    else:
        df_fq = pd.DataFrame(l_l_values, columns=["readname", "seq", "quality"])
        if flag_add_qname:
            df_fq["qname"] = list(
                e.split(" ", 1)[0][1:] for e in df_fq.readname.values
            )  # retrieve qname
        return df_fq

def Series_Subset(s, set_index):
    """# 2021-08-01 19:10:01
    subset index of a given series using a given set of index 'index'"""
    if not isinstance(set_index, set):
        set_index = set(set_index)
    return s[list(True if e in set_index else False for e in s.index.values)]

class Map(object):
    def __init__(self, dict_a2b):
        self.dict_a2b = dict_a2b

    def a2b(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return np.nan

    def a2b_if_mapping_available_else_Map_a2a(self, a):
        if a in self.dict_a2b:
            return self.dict_a2b[a]
        else:
            return a