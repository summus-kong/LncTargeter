#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import logging
from datetime import timedelta
from src import utils as ut
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from src import parse_sequence as ps

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def get_bwtResult(file):
    """
    Read bowtie result file as data frame
    :param file: bowtie result file
    :return: logical
    """
    # check input
    if not ut.check_files_exist(file):
        logging.error('Please check if bowtie align program done')
        raise FileNotFoundError('Please check if bowtie align program done')

    bwt_names = ['query_id', 'strand', 'ref_chr', 'first_position',
                 'sequence', 'quality', 'mapper_numbers', 'mismatch']
    logging.info('Read Bowtie align result file')
    bwt_result = pd.read_csv(file, header=None, sep='\t', names=bwt_names, dtype={'ref_chr': object})
    bwt_result['end_position'] = bwt_result.apply(lambda x: x['first_position'] + len(x['sequence']), axis=1)

    return bwt_result


def get_fs_result(fs):
    """
    Read in the TXT file from which the sub sequence was extracted
    :param fs: TXT file of sub sequence
    :return: data frame
    """
    # check input
    if not ut.check_files_exist(fs):
        logging.error('Please check if extracting sub sequence program done')
        raise FileNotFoundError('Please check if extracting sub sequence program done')

    logging.info('Read sub sequence TXT file')
    fs_usecols = ['id', 'sequence_fs', 'cm_value', 'RNA-RNA_deltaG', 'start_base', 'end_base', 'sub_GC',
                  'sequence_length']
    fs_df = pd.read_table(fs, sep='\t', usecols=fs_usecols)

    return fs_df


def merge_res(bwt_df, fs_df):
    """
    Merge the result file produced by extracting sub sequence and the bowtie align result file
    :param bwt_df: bowtie result data frame
    :param fs_df: sub sequence data frame
    :return: merged data frame
    """
    # check input
    if bwt_df.empty or fs_df.empty:
        logging.error('Please check input arguments of merge_res function')
        raise ValueError('Please check input arguments of merge_res function')

    logging.info('Merge bowtie align file and subsequences file for XGBoost model')
    merge_df = pd.merge(bwt_df, fs_df, left_on='query_id', right_on='id', how='inner')
    del merge_df['id']

    return merge_df


def base_sub(string, p, c):
    """
    Replaces a single character in a string
    :param string: string
    :param p: character to be converted
    :param c: A character used for conversion
    :return: converted string
    """
    # check input
    if not string or not p or not c:
        logging.error('Please check input arguments of base_sub function')
        raise ValueError('Please check input arguments of base_sub function')

    new = []
    for s in string:
        new.append(s)
    new[p] = c
    return ''.join(new)


def multi_sub(string, p, c):
    """
    Replace multiple characters in a string
    :param string: string
    :param p: characters to be converted
    :param c: characters used for conversion
    :return: converted string
    example: multi_sub['asddf', [2, 4], ['P', 'F']
    """
    # check input
    if not string or not p or not c:
        logging.error('Please check input arguments of multi_sub function')
        raise ValueError('Please check input arguments of multi_sub function')

    new = []
    for s in string:
        new.append(s)
    for index, point in enumerate(p):
        new[point] = c[index]

    return ''.join(new)


def intervaloverlap(list1, list2, overlap=0):
    """
    Determine whether the two intervals overlap
    :param list1: list one contain start and end
    :param list2: list two contain start and end
    :param overlap: overlapping length
    :return: logical
    """
    # check input
    if not list1 or not list2:
        logging.error('Please check input arguments of intervaloverlap function')
        raise ValueError('Please check input arguments of intervaloverlap function')

    l1 = np.max([list1[0], list2[0]])
    l2 = np.min([list1[1], list2[1]])
    if l2 - l1 > overlap:
        return True
    else:
        return False


def filter_bwt(bwt_df):
    """
    Screen the results of the bwt alignment. If multiple subsequences are aligned to the same position of
    the same chirp site, the subsequence with the least mismatches and the longest will be selected;
    Align to different positions of the same chirp region, then keep them;
    After screening at the same chirp region, screen on lncRNA based on the region
    :param bwt_df: bowtie result data frame
    :return: filtered bowtie result data frame
    """
    # check input
    if not bwt_df:
        logging.error('Please check input arguments of filter_bwt function')
        raise FileNotFoundError('Please check input arguments of filter_bwt function')

    logging.info('=== Screen bowtie align result ===')
    start_time = time.time()
    bwt_df = bwt_df.sort_values(by=['ref_chr', 'strand', 'first_position'],
                                ignore_index=True,
                                ascending=True)
    df_idx = []
    idx = 0
    while idx < bwt_df.shape[0]:
        tem_idx = []  # data frame index
        tem_seq = []  # base sequence length
        tem_mismatch = []  # number of mismatch
        a = [bwt_df.iloc[idx, 3], bwt_df.iloc[idx, 3] + len(bwt_df.iloc[idx, 4])]  # interval a
        out_lncsite = [bwt_df.iloc[idx, 10], bwt_df.iloc[idx, 11]]

        for inner_idx in bwt_df.index[idx:]:
            if bwt_df.iloc[inner_idx, 2] == bwt_df.iloc[idx, 2] and \
                bwt_df.iloc[inner_idx, 1] == bwt_df.iloc[idx, 1]:
                b = [bwt_df.iloc[inner_idx, 3],
                     bwt_df.iloc[inner_idx, 3] + len(bwt_df.iloc[inner_idx, 4])]  # interval b
                inner_lncsite = [bwt_df.iloc[inner_idx, 10], bwt_df.iloc[inner_idx, 11]]
                if intervaloverlap(a, b, 0) and intervaloverlap(out_lncsite, inner_lncsite, 1):
                    tem_idx.append(inner_idx)
                    tem_seq.append(len(bwt_df.iloc[inner_idx, 4]))
                    tem_mismatch.append(str(bwt_df.iloc[inner_idx, 7]).count(';') + 1)
                if not intervaloverlap(a, b, 0) and intervaloverlap(out_lncsite, inner_lncsite, 1):
                    break

        if tem_idx and tem_seq:
            # The element subscript that stores the longest sequence
            lidx = [i for i, x in enumerate(tem_seq) if x == max(tem_seq)]
            minmismatch = 4
            # filter the longest sequence with the least mismatches
            for tmp in lidx:
                if tem_mismatch[tmp] < minmismatch:
                    minmismatch = tem_mismatch[tmp]
            df_idx.append(tem_idx[tem_mismatch.index(minmismatch)])

        idx += len(tem_idx)

    df_idx = np.sort(pd.unique(df_idx))
    tem_bwt_df = bwt_df.iloc[df_idx].reset_index(drop=True)
    end_time = time.time()
    timediff = round(end_time - start_time)

    logging.info('Process screen bowtie align file done')
    logging.info("Time taken:" + str(timedelta(seconds=timediff)))

    return tem_bwt_df


def get_sequence(bwt):
    """
    Prepare a file to calculate the deltaG
    :param bwt: bowtie align data frame
    :return: bowtie align data frame
    """
    # check input
    if not bwt:
        logging.error('Please check input arguments of get_sequence function')
        raise ValueError('Please check input arguments of get_sequence function')

    start_time = time.time()
    logging.info('=== Prepare file to calculate deltaG ===')
    comp_sequence = []
    length = []
    strand = bwt['strand'].values
    sequence = bwt['sequence'].values
    mismatch = bwt['mismatch'].values
    for index, content in enumerate(sequence):
        tem_len = len(content)
        if pd.isna(mismatch[index]):
            tem_seq = str(Seq(content).complement())
        elif str(mismatch[index]).count(',') == 0:
            position = int(str(mismatch[index]).split(':')[0])
            base = str(mismatch[index]).split(':')[1].split('>')[0]
            tem_seq = str(Seq(content).complement())
            if strand[index] == '+':
                tem_seq = base_sub(tem_seq, position, base)
            else:
                tem_seq = base_sub(tem_seq, -(position + 1), base)
        else:
            poses = []
            bases = []
            for i in str(mismatch[index]).split(','):
                poses.append(int(i.split(':')[0]))
                bases.append(i.split(':')[1].split('>')[0])
            tem_seq = str(Seq(content).complement())
            if strand[index] == '+':
                tem_seq = multi_sub(tem_seq, poses, bases)
            else:
                tem_seq = multi_sub(tem_seq, [-(x + 1) for x in poses], bases)
        comp_sequence.append(tem_seq)
        length.append(tem_len)
    bwt['comp_sequence'] = comp_sequence
    bwt['length(bp)'] = length
    bwt['mismatch'].fillna('NA', inplace=True)
    end_time = time.time()
    timediff = round(end_time - start_time)

    logging.info('Process done')
    logging.info("Time taken:" + str(timedelta(seconds=timediff)))

    return bwt
