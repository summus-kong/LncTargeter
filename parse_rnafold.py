#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Date: 2021.10.12
@Author: Weibo Hou
@Description: Calculate and extract base sequence information
"""

import os
from src import utils as ut
import logging
import time
import src
from datetime import timedelta
import pandas as pd
import numpy as np
import random
import string
from src import parse_sequence as ps

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )
NN='rfile'
d1 = src.get_data(NN)
NNRR = np.load(d1, allow_pickle=True)
NNRR = NNRR.item()


def openfile(file):
    """
    Open RNAfold result file
    :param file: RNA fold result file
    :return: return a list divided by >
    """
    logging.info('Open RNAfold result file')
    # check input file
    if not ut.check_files_exist(file):
        logging.error('Failed to open RNAfold result file')
        raise ValueError('Please check RNAfold server')

    if ut.check_files_exist(file):
        with open(file, 'r') as f:
            lines = f.readlines()
            pline = [i.rstrip() for i in lines]
            flines = [pline[i: i + 3] for i in range(0, len(pline), 3)]

        return flines


def energy(seq, align):
    """
    Calculate RNA-RNA free energy
    :param seq: input sequence
    :param align: input alignment file
    :return: RNA-RNA free energy
    """
    # check input
    if not seq or not align:
        logging.error('Please check input arguments of energy function')
        raise ValueError('Please check input arguments')

    fe = NNRR['init']
    for base_idx in range(len(seq) - 1):
        incr = NNRR[(seq[base_idx: base_idx + 2], align[base_idx: base_idx + 2])]
        fe += incr

    return round(fe, 3)


def getRRenergy(rnafold, mincm=0.5, minlength=17, maxlength=50, myseed=1):
    """
    Calculate free energy by RNAfold result
    :param myseed: seed value
    :param maxlength: The maximum length of a sub sequence
    :param mincm: minimum complexity value for filtering sub sequence
    :param minlength: minimum sequence length for filtering sub sequence
    :param rnafold: RNAfold file
    :return: dict sequence id and sequence information
    """
    # check input
    if not rnafold:
        logging.error('Please check input arguments for parse RNAfold result file')
        raise ValueError('Please check input arguments')

    start_time = time.time()
    seqs = {}
    for info in rnafold:
        seqid = info[0]
        sequence = info[1].upper()
        seqalign = info[2].replace('(', ')')
        for step in range(minlength, maxlength + 1):
            start = 0
            for i in range(len(sequence) + 1 - step):
                seq = sequence[start: step + i]
                sf = ps.SequenceInformation()
                cm = sf.topological_entropy(seq)
                fe = energy(seq, seqalign[start: step + i])
                gc = round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)
                if cm >= mincm:
                    np.random.seed(myseed)
                    random_str = ''.join(random.sample(string.ascii_letters + string.digits, 5))
                    seqs[seqid.replace(" ", "").lstrip('>') + '_' + random_str] = [seq.replace('U', 'T'), cm, fe,
                                                                                   len(seq), start, step + i, gc]
                start += 1

    end_time = time.time()
    timediff = round(end_time - start_time)
    logging.info('Parse RNA secondary structure result file done')
    logging.info("Time taken of parse RNA secondary structure result file:" + str(timedelta(seconds=timediff)))

    return seqs


def out_csv(seqs, csv):
    """
    Output to csv which include sequence id， sequence cm value， sequence RNA-RNA free energy and length
    :param seqs: getRRenergy function result
    :param csv: tab-separated CSV file
    :return: logical
    """
    logging.info(f'Output extracted subsequences and related feature file: {csv}')
    flag = False
    # check input
    if not seqs:
        logging.error('Please check input arguments for out_csv function')
        raise ValueError('Please check input arguments')
    # check output path
    out_path = ut.get_file_dir(csv)
    if not ut.check_paths_exist(out_path):
        logging.info(f'Creating output path: {out_path}')
        if ut.mkdir(out_path):
            logging.info('Created path done')

    df = pd.DataFrame.from_dict(seqs, orient='index').reset_index()
    index = ['id', 'sequence_fs', 'cm_value', 'RNA-RNA_deltaG', 'sequence_length', 'start_base', 'end_base', 'sub_GC']
    df.columns = index
    df.to_csv(csv, sep='\t', index=False)

    if ut.check_files_exist(csv):
        flag = True

    return flag


def to_fa(seqs, fa):
    """
    Output to fasta file
    :param seqs: getRRenergy function result
    :param fa: output fasta file
    :return: logical
    """
    logging.info(f'Output extracted subsequences fasta file: {fa}')
    flag = False
    # check input
    if not seqs:
        logging.error('Please check input arguments for out_csv function')
        raise ValueError('Please check input arguments')
    # check output path
    out_path = ut.get_file_dir(fa)
    if not ut.check_paths_exist(out_path):
        logging.info(f'Creating output path: {out_path}')
        if ut.mkdir(out_path):
            logging.info('Created path done')

    with open(fa, 'w') as file:
        for i, j in seqs.items():
            f = '>' + i + '\n' + j[0] + '\n'
            file.write(f)

    if ut.check_files_exist(fa):
        flag = True

    return flag
