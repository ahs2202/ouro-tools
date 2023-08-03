from bitarray import bitarray
import numpy as np


def Find(ba, val=1):
    """deprecated
    # 2022-06-12 21:09:22
    search the given value in the bitarray iteratively and return the start position of the occurrences
    """
    len_ba = len(ba)
    l_int_pos_occurrence = []
    int_pos_start = 0
    while int_pos_start < len_ba:
        int_pos_occurrence = ba.find(val, int_pos_start + 1)
        if int_pos_occurrence < 0:
            break
        else:
            l_int_pos_occurrence.append(int_pos_occurrence)
            int_pos_start = int_pos_occurrence
    return l_int_pos_occurrence


def find(ba, val=1):
    """# 2022-07-03 21:41:36
    generator that returns the location of 'val' from the start of the bitarray.
    Use this function for memory-efficient iteration of position of active entries in a bitarray (np.int64 takes 64 times more memory than an entry in bitarray, which is 1 bit).
    """
    len_ba = len(ba)
    if len(ba) == 0:  # handle empty bitarray input
        return
    int_pos_start = 0
    if ba[0] == val:  # handle the case when the first bit is equal to 'val'
        yield int_pos_start
    while int_pos_start < len_ba:
        int_pos_occurrence = ba.find(val, int_pos_start + 1)
        if int_pos_occurrence < 0:
            break
        else:
            yield int_pos_occurrence
            int_pos_start = int_pos_occurrence


# def Find_Segment( ba, background = 1, flag_use_bitwise_operation = True ) :
#     ''' # 2022-06-12 19:36:38
#     find segment of a given bitarray for the given background 'background'. for example, when background = 1, find segment of 0

#     'flag_use_bitwise_operation' : use bitwise operation. useful when the length of bitarray is very long
#     '''
#     len_ba = len( ba )
#     if flag_use_bitwise_operation :
#         ''' the implementation using bitwise operation '''
#         mask = ( ba ^ ba >> 1 ) # find boundary of value change
#         l_int_pos = Find( mask, 1 )
#         flag_start_is_a_segment = ba[ 0 ] != background
#         if flag_start_is_a_segment :
#             l_int_pos = [ 0 ] + l_int_pos
#         if len( l_int_pos ) % 2 != 0 :
#             l_int_pos += [ len_ba ]
#         return np.array( l_int_pos ).reshape( ( int( len( l_int_pos ) / 2 ), 2 ) )
#     else :
#         ''' the implementation using simple iteration '''
#         int_pos = 0
#         int_seg_start = None
#         l_seg = [ ]
#         while int_pos < len_ba :
#             if int_seg_start is None and ba[ int_pos ] != background :
#                 int_seg_start = int_pos
#             elif int_seg_start is not None and ba[ int_pos ] == background :
#                 l_seg.append( [ int_seg_start, int_pos ] )
#                 int_seg_start = None
#             int_pos += 1
#         if int_seg_start is not None :
#             l_seg.append( [ int_seg_start, len_ba ] )
#         return l_seg


def Find_Segment(ba, background=1):
    """# 2022-06-12 19:36:38
    find segment of a given bitarray for the given background 'background'. for example, when background = 1, find segment of 0

    # 2022-06-22 18:36:45
    updated to a new implementation that does not require generating a copy of the bitarray.
    """

    def toggle_bit(bit):
        return (bit + 1) % 2  # toggle the bit

    # initialize
    l_int_pos = []
    state_of_interest = toggle_bit(background)  # look for a non-background state
    int_pos = 0

    while True:
        int_pos = ba.find(state_of_interest, int_pos)
        if int_pos < 0:
            break
        l_int_pos.append(int_pos)
        state_of_interest = toggle_bit(
            state_of_interest
        )  # toggle the state of interest
    if len(l_int_pos) % 2 != 0:
        l_int_pos.append(len(ba))
    return np.array(l_int_pos).reshape(
        (int(len(l_int_pos) / 2), 2)
    )  # reshape a list of int_pos to list of segments


def Retrieve_Integer_Indices(ba, background=0):
    """deprecated (slow)
    # 2022-06-23 08:31:45
    returns a list of integer indices of non-background bits in a given bitarray 'ba'

    'background' : a bit representing background. either 0 or 1
    """
    # compose a list of integer indices of active rows after applying filter
    l = []
    for st, en in Find_Segment(
        ba, background=background
    ):  # retrieve active segments from bitarray filter
        l.extend(range(st, en))  # retrieve integer indices of the active rows
    return l


def to_array(ba):
    """# 2022-06-24 23:54:12
    return a boolean numpy array of the given bitarray"""
    return np.frombuffer(ba.unpack(), dtype=bool)


def to_bitarray(arr_bool):
    """# 2022-06-27 15:29:56
    convert numpy boolean array to bitarray
    """
    ba = bitarray()
    ba.pack(arr_bool.tobytes())
    return ba


def from_integer_indices_to_bitarray(l_int_index, length):
    """# 2022-07-01 11:21:12
    convert list of integer indices to bitarray

    'length' : length of the output bitarray object
    """
    ba = bitarray(length)
    ba.setall(0)
    for int_index in l_int_index:
        ba[int_index] = True
    return ba


def to_integer_indices(ba):
    """# 2022-07-01 21:38:17
    retrieve integer indices of the active entries of the given bitarray

    'ba' : input bitarray object
    """
    return np.where(to_array(ba))[0]


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


def detect_boolean_mask(ba):
    """# 2022-08-10 23:29:06
    detect boolean mask by looking up to 10 values
    """
    # extract the first row from the ndarray data
    if isinstance(ba, np.ndarray) and len(ba.shape) > 1:
        for i in range(len(ba.shape) - 1):
            ba = ba[0]
    if not hasattr(ba, "__iter__"):  # 'ba' should be iterable to be a boolean mask
        return False
    ba = ba[:10]
    # if the length of array of interest is <= 2 and only consists of 0, 1, consider it as an integer array, not boolean array
    if len(ba) <= 2 and set(COUNTER(ba[:10])).issubset({0, 1}):
        return False
    return set(COUNTER(ba[:10])).issubset({0, 1, True, False})


def convert_mask_to_bitarray(ba):
    """# 2022-07-03 00:54:46
    convert boolean mask to a bitarray
    """
    # handle when a list type has been given (convert it to np.ndarray)
    if isinstance(ba, list):
        ba = np.array(ba, dtype=bool)
    # handle when a numpy ndarray has been given (convert it to bitarray)
    if isinstance(ba, np.ndarray):
        ba = to_bitarray(ba)
    assert isinstance(ba, bitarray)  # check the return value is bitarray object
    return ba


def convert_mask_to_array(ba):
    """# 2022-07-03 00:54:46
    convert boolean mask to a array
    """
    """ handle non-bitarray mask types """
    # handle when a list type has been given (convert it to np.ndarray)
    if isinstance(ba, list):
        ba = np.array(ba, dtype=bool)
    # handle when a numpy ndarray has been given (convert it to np.ndarray)
    if isinstance(ba, bitarray):
        ba = to_array(ba)
    assert isinstance(ba, np.ndarray)  # check the return value is np.ndarray object
    return ba
