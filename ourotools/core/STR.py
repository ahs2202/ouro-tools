import numpy as np
import regex

from bitarray import bitarray

from . import BA

# In[ ]:


setting__n_characters__Insert_characters_every_n_characters = 20


def Insert_characters_every_n_characters(
    a_string, n_characters=None, insert_characters="\n"
):
    """Insert new line characters in the sequences every 'n_characters' characters"""
    n_characters = (
        setting__n_characters__Insert_characters_every_n_characters
        if n_characters is None
        else n_characters
    )
    return insert_characters.join(
        list(
            a_string[n_characters * index_line : n_characters * (index_line + 1)]
            for index_line in np.arange(int((len(a_string) - 1) / n_characters) + 1)
        )
    )  # Since there should not be additional no new additional character at the end for the sequences with 60 amino acid, ( len( seq ) - 1 ) / 60 will be used


def Replace_a_character_at_an_index(text, index, replacement):
    return text[:index] + replacement + text[index + 1 :]


# In[ ]:


setting__n_characters__replace_a_character_every_n_characters = 11


def Replace_a_character_every_n_characters(
    a_string, n_characters=None, replacement_characters=""
):
    """Insert new line characters in the sequences every 'n_characters' characters"""
    n_characters = (
        setting__n_characters__replace_a_character_every_n_characters
        if n_characters is None
        else n_characters
    )
    return replacement_characters.join(
        list(
            a_string[n_characters * index_line : n_characters * (index_line + 1) - 1]
            for index_line in np.arange(int((len(a_string) - 1) / n_characters) + 1)
        )
    )  # Since there should not be additional no new additional character at the end for the sequences with 60 amino acid, ( len( seq ) - 1 ) / 60 will be used


# In[ ]:


def Find_all(a_str, sub):
    """Find all positions of a substring 'sub' in a given string 'a_str'"""
    start = 0
    l_positions = list()
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            break
        l_positions.append(start)
        start += len(sub)  # use start += 1 to find overlapping matches
    return l_positions


# In[ ]:


def Generate_substring_indexes(substring, string):
    """Generate indices of where substring begins in string
    >>> list(find_substring('me', "The cat says meow, meow"))
    [13, 19]"""
    last_found = -1  # Begin at -1 so the next position to search from is 0
    while True:
        # Find next index of substring, by starting after its last known position
        last_found = string.find(substring, last_found + 1)
        if last_found == -1:
            break  # All occurrences have been found
        yield last_found


# In[ ]:


def Locate_substring_last_index(substring, string):
    """return an index for the last substring in a string"""
    last_found = -1  # Begin at -1 so the next position to search from is 0
    while True:
        # Find next index of substring, by starting after its last known position
        index = string.find(substring, last_found + 1)
        if index == -1:
            return last_found
        else:
            last_found = index


# In[4]:


def To_python_compatible_str(
    a_string,
    dict_replacement={
        "%": "_Percent_",
        "+": "_Plus_",
        "-": "_Minus_",
        "&": "_and_",
        "=": "_Equal_",
    },
):
    """convert a string into python-compatible string."""
    l_incompatible_character = [
        " ",
        "?",
        "(",
        ")",
        "&",
        "%",
        "/",
        ",",
        ":",
        ".",
        "-",
        "+",
        "[",
        "]",
        "#",
        "=",
        "\n",
        '"',
        "\\",
        "|",
        "?",
        "*",
    ]
    for incompatible_character in l_incompatible_character:
        a_string = a_string.replace(
            incompatible_character, dict_replacement.get(incompatible_character, "_")
        )
    return a_string


# In[ ]:


def To_path_compatible_str(a_string):
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
        .replace(" ", "_")
        .replace(";", "_")
    )


# In[ ]:


def Read_Field_Locations(
    str_field_locations,
    character_indicating_fields="-",
    character_indicating_non_fields=" ",
):
    """return start and end position of masks annotating field locations of a given strings"""
    int_position = 0
    l_l_annotation_start_end = list()
    while True:
        str_field_locations_start = str_field_locations[int_position:]
        if "-" in str_field_locations_start:
            int_annotation_start = (
                str_field_locations_start.find(character_indicating_fields)
                + int_position
            )  # retrive field annotation start position
        else:
            break  # finish the loop if there is no annotation left
        str_field_locations_end = str_field_locations[int_annotation_start:]
        int_annotation_end = (
            str_field_locations_end.find(character_indicating_non_fields)
            + int_annotation_start
            if " " in str_field_locations_end
            else len(str_field_locations)
        )  # retrive field annotation end position
        l_l_annotation_start_end.append([int_annotation_start, int_annotation_end])
        int_position = int_annotation_end
    return l_l_annotation_start_end


# In[2]:


def Search_Subsequence(
    str_seq, str_subseq, error_rate=0.1, index_start_window=None, index_end_window=None
):
    """# 2021-02-04 18:12:04
    search a sequence with a given subsequence and return the integer index after the subsequence (for example, when subsequence presents from 10-19, return 20)
    'error_rate' : maximum allowed error_rate (including substitutions, insertions and deletions)
    'index_start_window' :
    """
    # set default values
    str_matched_subsequence = None
    index_start_subsequence = -1  # 0-based coordinates
    index_end_subsequence = -1  # 0-based coordinates
    int_num_errors = -1
    # search all positions of a given sequence by default
    len_seq = len(str_seq)
    if index_start_window is None:
        index_start_window = 0
    if index_end_window is None:
        index_end_window = len_seq
    # handle negative indices (indices from the end of the sequence)
    if index_start_window < 0:
        index_start_window = len_seq + index_start_window
    if index_end_window < 0:
        index_end_window = len_seq + index_end_window
    str_seq = str_seq[index_start_window:index_end_window]
    if str_subseq in str_seq:
        str_matched_subsequence = str_subseq  # return the index after the subsequence
        index_start_subsequence = index_start_window + str_seq.find(str_subseq)
        index_end_subsequence = index_start_subsequence + len(str_subseq)
        int_num_errors = 0
    else:
        int_max_num_allowed_errors = int(
            len(str_subseq) * error_rate
        )  # set maximum number of errors
        for int_num_allowed_errors in range(1, int_max_num_allowed_errors + 1):
            l_match = regex.findall(
                f"({str_subseq}){{e<={int_num_allowed_errors}}}", str_seq
            )
            if len(l_match) > 0:
                str_matched_subsequence = l_match[0]
                index_start_subsequence = index_start_window + str_seq.find(
                    str_matched_subsequence
                )
                index_end_subsequence = index_start_subsequence + len(
                    str_matched_subsequence
                )
                int_num_errors = int_num_allowed_errors
                break

    return {
        "index_start_subsequence": index_start_subsequence,
        "index_end_subsequence": index_end_subsequence,
        "matched_subsequence": str_matched_subsequence,
        "num_errors": int_num_errors,
    }  # return dictionary containing the index after the subsequence


def Search_Subsequences(
    str_seq,
    str_subseq,
    error_rate=0.1,
    index_start_window=None,
    index_end_window=None,
    int_min_num_errors=0,
):
    """# 2022-06-15 23:48:25
    perform an exhuastive search of a sequence with a given subsequence and return the integer index after the subsequence (for example, when subsequence presents from 10-19, return 20)
    ignore overlapping matches. if one region was matched with the subsequence, the subsequent search will be done in the sequences excluding the region

    'error_rate' : maximum allowed error_rate (including substitutions, insertions and deletions)
    'index_start_window', 'index_end_window' : start and end position of the region to search the subsequence

    'int_min_num_errors' : an internal argument. the minimum number of errors for the search. used in the recursive calling of this function internally.
    """
    l_sr = []  # initialize search results

    def _initialize_search_result_():
        # set default values
        return {
            "index_start_subsequence": None,
            "index_end_subsequence": -1,
            "matched_subsequence": -1,
            "num_errors": -1,
        }  # 0-based coordinates

    def _add_search_results_of_a_subsequence_(str_subseq, int_num_errors):
        len_subseq = len(str_subseq)
        for index_found in Find_all(str_seq_of_a_given_window, str_subseq):
            ba[
                index_found : index_found + len_subseq
            ] = 1  # set the mask to mark the region matched with the subsequence
            sr = _initialize_search_result_()
            sr[
                "matched_subsequence"
            ] = str_subseq  # return the index after the subsequence
            sr["index_start_subsequence"] = index_start_window + index_found
            sr["index_end_subsequence"] = sr["index_start_subsequence"] + len_subseq
            sr["num_errors"] = int_num_errors
            l_sr.append(sr)

    """ indexing """
    # search all positions of a given sequence by default
    len_seq = len(str_seq)
    if index_start_window is None:
        index_start_window = 0
    if index_end_window is None:
        index_end_window = len_seq
    # handle negative indices (indices from the end of the sequence)
    if index_start_window < 0:
        index_start_window = len_seq + index_start_window
    if index_end_window < 0:
        index_end_window = len_seq + index_end_window
    # retrieve a sequence of the given window
    str_seq_of_a_given_window = str_seq[index_start_window:index_end_window]

    """ handle premature exit condition """
    if (
        len(str_seq_of_a_given_window) < len(str_subseq) - int_min_num_errors
    ):  # if the searched region is smaller than the query sequence (considering errors), exit without performing the search
        return l_sr

    """ initialize bitarray mask """
    # 'ba' : an internal argument. a bitarray object of 'str_seq_of_a_given_window', marking the regions were matched subsequences were found (marked by 1). these region will be ignored
    ba = bitarray(len(str_seq_of_a_given_window))
    ba.setall(0)

    """ performing search for 'int_num_allowed_errors' """
    int_max_num_allowed_errors = int(
        len(str_subseq) * error_rate
    )  # retrieve maximum number of allowed errors
    if int_min_num_errors == 0:
        _add_search_results_of_a_subsequence_(str_subseq, 0)
    elif int_min_num_errors <= int_max_num_allowed_errors:
        for str_matched_subsequence in regex.findall(
            f"({str_subseq}){{e<={int_min_num_errors}}}", str_seq_of_a_given_window
        ):
            _add_search_results_of_a_subsequence_(
                str_matched_subsequence, int_min_num_errors
            )

    if (
        int_min_num_errors < int_max_num_allowed_errors
    ):  # recursive calls only when terminating condition (int_min_num_errors == int_max_num_allowed_errors) is not met
        for start_next_window, end_next_window in BA.Find_Segment(
            ba, background=1
        ):  # ignore matched regions in the subsequent searches
            l_sr.extend(
                Search_Subsequences(
                    str_seq,
                    str_subseq,
                    error_rate=error_rate,
                    index_start_window=index_start_window + start_next_window,
                    index_end_window=index_start_window + end_next_window,
                    int_min_num_errors=int_min_num_errors + 1,
                )
            )  # perform recursive search for the larger number of errors and collect the result
    return l_sr  # return the search result


# In[ ]:


def Find_stretch_of_a_character(a_str, a_char, int_len_threshold=1):
    """
    Find locations of stretches of a character
    'int_len_threshold' : a number of characters in a stretch to be classified as a 'stretch'
    return a list of [ index_start, index_end, length_of_stretch ] of identified stretches
    """
    # set default values
    flag_stretch_started = False
    index_start = None
    index = 0
    l_l = list()
    for current_char in a_str:
        if flag_stretch_started:
            if current_char != a_char:  # when stretch has been started and has ended
                flag_stretch_started = False
                len_stretch = index - index_start
                if len_stretch >= int_len_threshold:
                    l_l.append([index_start, index, len_stretch])
                index_start = None
        else:
            if current_char == a_char:  # when stretch has just started
                flag_stretch_started = True
                index_start = index
        index += 1
    if index_start is not None:  # if stretch of a character ends with the given string
        len_stretch = index - index_start
        if len_stretch >= int_len_threshold:
            l_l.append([index_start, index, len_stretch])
    return l_l
