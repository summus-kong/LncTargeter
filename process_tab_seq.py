#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import pandas as pd
from src import parse_sequence as ps
from src import utils as ut

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def read_bed_file(infile, inames, dtypes=None):
    """
    Read tab-separate file
    :param dtypes:
    :param infile: filename
    :param inames: colnames list
    :return: pandas data frame
    """
    # check input
    if not ut.check_files_exist(infile):
        logging.error('Please check input for intersect_gene function')
        raise FileNotFoundError
    if not isinstance(inames, list):
        logging.error('Please provide colnames list')
        raise ValueError
    # names = ['chr', 'start', 'end', 'sequence', 'length']
    logging.info(f'Reading file: {infile}')
    data = pd.read_csv(infile, header=None, sep='\t', names=inames, low_memory=False, dtype=dtypes)

    return data


def process_bed(indf, out_cols):
    """
    Calculate GC content and cm value
    :param out_cols: colnames
    :param indf: input pandas data frame
    :return: pandas data frame
    """
    # check input
    if indf.empty or not isinstance(out_cols, list):
        logging.error('Please provide valid parameter')
        raise ValueError

    gcs = []
    cms = []
    sequence = [x.upper() for x in indf['sequence'].values]
    logging.info('Calculate GC content, complexity...')
    for seq in sequence:
        sf = ps.SequenceInformation
        length = len(seq)
        cm = sf.topological_entropy(seq)
        if length > 0:
            gc = round((seq.count('G') + seq.count('C')) / length * 100, 2)
        else:
            gc = 0

        gcs.append(gc)
        cms.append(cm)
    col1, col2 = out_cols
    indf[col1] = gcs
    indf[col2] = cms
    del indf['sequence']

    return indf


def merge_table(t1, t2, keys):
    """
    Merge data frame
    :param t1: dt 1
    :param t2: dt 2
    :param keys: key words
    :return: merge df
    """
    # check input
    if t1.empty or t2.empty or not keys:
        logging.error('Please provide valid parameter')
        raise ValueError

    logging.info('Merge temporary file')
    merge_df = pd.merge(t1, t2, on=keys, how='inner')

    return merge_df
