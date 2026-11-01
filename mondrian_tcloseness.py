# coding:utf-8
"""
main module of mondrian
"""

# !/usr/bin/env python
# coding=utf-8

import pdb
import time
from datetime import datetime
from utils.utility import cmp_value, value
from functools import cmp_to_key

# warning all these variables should be re-inited, if
# you want to run mondrian with different parameters
__DEBUG = False
QI_LEN = 10
GL_K = 0
RESULT = []
QI_RANGE = []
QI_DICT = []
QI_ORDER = []

'''
TODO (t-closeness Part 0):
t-closeness threshold in [0, 1]; smaller means stricter privacy; update this parameter for experiments.
'''
T_THRESHOLD = 0.5

class Partition(object):

    """
    Class for Group (or EC), which is used to keep records
    self.member: records in group
    self.low: lower point, use index to avoid negative values
    self.high: higher point, use index to avoid negative values
    self.allow: show if partition can be split on this QI
    """

    def __init__(self, data, low, high):
        """
        split_tuple = (index, low, high)
        """
        self.low = list(low)
        self.high = list(high)
        self.member = data[:]
        self.allow = [1] * QI_LEN

    def add_record(self, record, dim):
        """
        add one record to member
        """
        self.member.append(record)

    def add_multiple_record(self, records, dim):
        """
        add multiple records (list) to partition
        """
        for record in records:
            self.add_record(record, dim)

    def __len__(self):
        """
        return number of records
        """
        return len(self.member)
    

def get_normalized_width(partition, index):
    """
    return Normalized width of partition
    similar to NCP
    """
    d_order = QI_ORDER[index]
    width = value(d_order[partition.high[index]]) - value(d_order[partition.low[index]])
    if width == QI_RANGE[index]:
        return 1
    return width * 1.0 / QI_RANGE[index]


def choose_dimension(partition):
    """
    choose dim with largest norm_width from all attributes.
    This function can be upgraded with other distance function.
    """
    max_width = -1
    max_dim = -1
    for dim in range(QI_LEN):
        if partition.allow[dim] == 0:
            continue
        norm_width = get_normalized_width(partition, dim)
        if norm_width > max_width:
            max_width = norm_width
            max_dim = dim
    if max_width > 1:
        pdb.set_trace()
    return max_dim


def frequency_set(partition, dim):
    """
    get the frequency_set of partition on dim
    """
    frequency = {}
    for record in partition.member:
        try:
            frequency[record[dim]] += 1
        except KeyError:
            frequency[record[dim]] = 1
    return frequency


def find_median(partition, dim):
    """
    find the middle of the partition, return split_val
    """
    # use frequency set to get median
    frequency = frequency_set(partition, dim)
    split_val = ''
    next_val = ''
    value_list = list(frequency.keys())
    value_list.sort(key=cmp_to_key(cmp_value))
    total = sum(frequency.values())
    middle = total // 2
    if middle < GL_K or len(value_list) <= 1:
        try:
            return '', '', value_list[0], value_list[-1]
        except IndexError:
            return '', '', '', ''
    index = 0
    split_index = 0
    for i, qi_value in enumerate(value_list):
        index += frequency[qi_value]
        if index >= middle:
            split_val = qi_value
            split_index = i
            break
    else:
        print("Error: cannot find split_val")
    try:
        next_val = value_list[split_index + 1]
    except IndexError:
        # there is a frequency value in partition
        # which can be handle by mid_set
        # e.g.[1, 2, 3, 4, 4, 4, 4]
        next_val = split_val
    return (split_val, next_val, value_list[0], value_list[-1])

def merge_qi_value(x_left, x_right, connect_str='~'):
    '''Connect the interval boundary value as a generalized interval and return the result as a string
    return:
        result:string
    '''
    if isinstance(x_left, (int, float)):
        if x_left == x_right:
            result = '%d' % (x_left)
        else:
            result = '%d%s%d' % (x_left, connect_str, x_right)
    elif isinstance(x_left, str):
        if x_left == x_right:
            result = x_left
        else:
            result = x_left + connect_str + x_right
    elif isinstance(x_left, datetime):
        # Generalize the datetime type value
        begin_date = x_left.strftime("%Y-%m-%d %H:%M:%S")
        end_date = x_right.strftime("%Y-%m-%d %H:%M:%S")
        result = begin_date + connect_str + end_date
    return result

def compute_emd(partition):
    """
    TODO (t-closeness Part 2):
    Compute EMD between the partition's SA distribution and GLOBAL_SA_DIST.
    GLOBAL_SA_DIST is computed in the init function, and its calculation corresponds to (t-closeness Part 1).
    For categorical SA, EMD = 0.5 * sum over all SA values of |p_v - q_v|.
    Return a float in [0, 1].
    """
    sa_index = len(partition.member[0]) - 1 if partition.member else -1
    # Step 1: count SA values inside this partition
    # Step 2: for every value in GLOBAL_SA_DIST, compute |p_v - q_v|
    # Step 3: return 0.5 * sum of absolute differences
    # your code here

def satisfy_t_closeness(partition, t):
    if partition is None or not partition.member:
        return True
    return compute_emd(partition) <= t + 1e-9

def anonymize_strict(partition):
    """
    recursively partition groups until not allowable
    """
    allow_count = sum(partition.allow)
    # only run allow_count times
    if allow_count == 0:
        RESULT.append(partition)
        return
    for index in range(allow_count):

        '''
        TODO (t-closeness Part 3):
        Reject this split if either child violates t-closeness.
        If rejected, mark this dim as not allowed, then continue to try next dim.
        You only need to insert the t-closeness condition check right after the k-anonymity check in the code from Experiment 1.
        Hint: use satisfy_t_closeness(lhs, T_THRESHOLD) and satisfy_t_closeness(rhs, T_THRESHOLD) above.
        '''
        # your code here

        anonymize_strict(lhs)
        anonymize_strict(rhs)
        return
    
    RESULT.append(partition)


def anonymize_relaxed(partition):
    """
    recursively partition groups until not allowable
    """
    if sum(partition.allow) == 0:
        # can not split
        RESULT.append(partition)
        return
    
    '''
    TODO (t-closeness Part 4):
    Reject this split if either child violates t-closeness.
    If rejected, mark this dim as not allowed, then RETRY on the same partition.
    (This function has no for-loop, so there is no `continue` here.
    Follow the same retry pattern as the k-anonymity check above.)
    Hint: use satisfy_t_closeness(lhs, T_THRESHOLD) and satisfy_t_closeness(rhs, T_THRESHOLD).
    '''
    # your code here   

    anonymize_relaxed(lhs)
    anonymize_relaxed(rhs)


def init(data, k, QI_num=-1):
    """
    reset global variables
    """
    global GL_K, RESULT, QI_LEN, QI_DICT, QI_RANGE, QI_ORDER
    if QI_num <= 0:
        QI_LEN = len(data[0]) - 1
    else:
        QI_LEN = QI_num
    GL_K = k
    RESULT = []
    # static values
    QI_DICT = []

    global GLOBAL_SA_DIST
    '''
    TODO (t-closeness Part 1):
    Compute the global distribution of the sensitive attribute.
    The sensitive attribute is the LAST column: index = len(record) - 1.
    Store the result in the global variable GLOBAL_SA_DIST as {value: probability}.
    Example: if 6 records have SA values [A, A, B, C, C, C],
             then GLOBAL_SA_DIST = {'A': 2/6, 'B': 1/6, 'C': 3/6}.
    '''
    #your code here

    QI_ORDER = []
    QI_RANGE = []
    att_values = []
    for i in range(QI_LEN):
        att_values.append(set())
        QI_DICT.append(dict())
    for record in data:
        for i in range(QI_LEN):
            att_values[i].add(record[i])
    for i in range(QI_LEN):
        value_list = list(att_values[i])
        value_list.sort(key=cmp_to_key(cmp_value))
        QI_RANGE.append(value(value_list[-1]) - value(value_list[0]))
        QI_ORDER.append(list(value_list))
        for index, qi_value in enumerate(value_list):
            QI_DICT[i][qi_value] = index


def mondrian(data, k, relax=False, QI_num=-1):
    """
    Main function of mondrian, return result in tuple (result, (ncp, rtime)).
    data: dataset in 2-dimensional array.
    k: k parameter for k-anonymity
    QI_num: Default -1, which exclude the last column. Othewise, [0, 1,..., QI_num - 1]
            will be anonymized, [QI_num,...] will be excluded.
    relax: determine use strict or relaxed mondrian,
    Both mondrians split partition with binary split.
    In strict mondrian, lhs and rhs have not intersection.
    But in relaxed mondrian, lhs may be have intersection with rhs.
    """
    init(data, k, QI_num)
    result = []
    data_size = len(data)
    low = [0] * QI_LEN
    high = [(len(t) - 1) for t in QI_ORDER]
    whole_partition = Partition(data, low, high)
    # begin mondrian
    start_time = time.time()
    if relax:
        # relax model
        anonymize_relaxed(whole_partition)
    else:
        # strict model
        anonymize_strict(whole_partition)
    rtime = float(time.time() - start_time)
    # generalization result and
    # evaluation information loss
    ncp = 0.0
    dp = 0.0
    for partition in RESULT:
        rncp = 0.0
        for index in range(QI_LEN):
            rncp += get_normalized_width(partition, index)
        rncp *= len(partition)
        ncp += rncp
        dp += len(partition) ** 2
        for record in partition.member[:]:
            for index in range(QI_LEN):
                record[index] = merge_qi_value(QI_ORDER[index][partition.low[index]],
                                QI_ORDER[index][partition.high[index]])
            result.append(record)
    # If you want to get NCP values instead of percentage
    # please remove next three lines
    ncp /= QI_LEN
    ncp /= data_size
    ncp *= 100

    emds = [compute_emd(p) for p in RESULT]
    actual_t = max(emds) if emds else 0.0
    violations = sum(1 for e in emds if e > T_THRESHOLD + 1e-9)
    print("Target t = %.4f \nActual max EMD = %.4f \nNumber of t-closeness violations = %d \nNumber of equivalence classes = %d"
          % (T_THRESHOLD, actual_t, violations, len(RESULT)))
    if violations > 0:
        print("  [Warning] Some partitions exceed T_THRESHOLD. "
              "Check the t-closeness check in TODO 3 / 4b.")

    if __DEBUG:
        from decimal import Decimal
        print("Discernability Penalty=%.2E" % Decimal(str(dp)))
        print("size of partitions=%d" % len(RESULT))
        print("K = %d" % k)
        print("GCP = %.2f %%" % ncp)
    return (result, (ncp, rtime))
