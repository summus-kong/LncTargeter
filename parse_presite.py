#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.14
@Author: Weibo Hou
@Description: Extract presite information
"""

from src import parse_sequence as ps
import pandas as pd
import numpy as np
import logging
from src import utils as ut
import src

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )

NR='dfile'
d2 = src.get_data(NR)
NNRD = np.load(d2, allow_pickle=True)
NNRD = NNRD.item()


def parser_site(file):
    """
    Analysis of loci related characteristics
    :param file: presite file
    :return: data frame
    """
    logging.info('Prepare and calculate binding region characteristics')
    # check input file
    if not ut.check_files_exist(file):
        logging.error('Please check input')
        raise FileNotFoundError('Please check input arguments')

    names = ['chr', 'start', 'end', 'sequence_site']
    logging.info('Open the tab-separated intermediate file containing chromosome position information and binding '
                 'region sequences')
    assume_site = pd.read_csv(file, header=None, sep='\t', names=names)
    return assume_site


def calculate_deltaG(df):
    """
    Calculate sequence length,deltaG,cm,GC content
    :param df:
    :return:
    """
    # check input
    if df.empty:
        logging.error('Please provide valid parameter')
        raise ValueError

    fe_list = []
    length_list = []
    cm_list = []
    gc_list = []

    sequence = [x.upper() for x in df['sequence_site'].values]
    window = 2
    logging.info('Calculate binding region characteristics: eg GC content, complexity...')
    for seq in sequence:
        sf = ps.SequenceInformation
        fe = 2 * NNRD['init']
        length = len(seq)
        cm = sf.topological_entropy(seq)
        gc = round((seq.count('G') + seq.count('C')) / length * 100, 2)
        if length < window:
            fe_list.append(fe)
        else:
            start = 0
            end = 2
            while end <= length:
                subseq = str(seq[start:end])
                if subseq in NNRD.keys():
                    fe = fe + NNRD[subseq]
                start += 1
                end += 1
            fe_list.append(fe)
        length_list.append(length)
        cm_list.append(cm)
        gc_list.append(gc)

    df['RD_site'] = fe_list
    df['length_site'] = length_list
    df['cm_site'] = cm_list
    df['gc_site'] = gc_list

    return df


def merge_site_fs(fs, site_df):
    # names = ['chr', 'start', 'end', 'lncid', 'sequence_fs', 'cm_fs', 'RR_fs', 'start_fs', 'length_fs', 'gc_fs']
    # fs_df = pd.read_csv(fs, header=None, sep='\t', names=names)
    # fs is data frame
    # check input
    if fs.empty and site_df.empty:
        logging.error('Please provide valid parameter')
        raise ValueError

    logging.info('Merge DNA binding region associated file and subsequences file for XGBoost model')
    merge_df = pd.merge(fs, site_df, on=['chr', 'start', 'end'], how='inner')
    return merge_df


