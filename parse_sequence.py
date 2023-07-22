#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Date: 2021.10.11
@Author: Weibo Hou
@Description: Calculate and extract base sequence information
"""
import math
import time
from datetime import timedelta
import logging
import pandas as pd
from src import utils as ut
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


class SequenceInformation:

    @staticmethod
    def get_sequence(fa):
        """
        Get the potential target sequence id and sequence
        :param fa: input fasta file
        :return: dict of sequence id and sequence
        """
        sids = []
        seqs = []
        # check input fasta file
        if not ut.check_files_exist(fa):
            logging.error('Please provide a valid input fasta file')
            raise ValueError('Please check input fasta file')

        flag = 0
        if ut.check_files_exist(fa):
            flag = 1
        if flag:
            logging.info('Parse fasta sequence to dict for further analysis')
            for seq in SeqIO.parse(fa, "fasta"):
                tmp = str(seq.seq).upper().replace('U', 'T')
                sids.append(seq.id)
                seqs.append(tmp)
        seqdict = dict(zip(sids, seqs))

        return seqdict

    def get_GCcontent(self, seq):
        """
        Calculate sequence GC content
        :param seq: sequence object
        :return: GC content
        """
        if not seq:
            logging.error('Please check get_GCcontent function')
            raise ValueError('Please check get_GCcontent arguments')

        gc_content = round(GC(seq), 2)

        return gc_content

    @staticmethod
    def topological_entropy(seq):
        """
        compute topological entropy
        :param seq: string base sequence
        :return: complex value
        """
        if not seq:
            raise ValueError('Please check input sequence')

        length = len(seq)
        temp_n = math.log(length, 4)
        up_n = math.ceil(temp_n)
        down_n = math.floor(temp_n)
        if length >= math.pow(4, down_n) + down_n - 1:
            n = down_n
        else:
            n = up_n
        store = []
        sts = []
        cnt = 0
        if n > 0:
            while cnt < length - n:
                i = 0
                while i < n:
                    i_i = cnt + i
                    store.append(seq[i_i])
                    i = i + 1
                st = ''.join(store)
                sts.append(st)
                store.clear()
                cnt = cnt + 1
            stt = set(sts)
            kind_n = math.log(len(stt), 4)
            entropy = round(kind_n / n, 2)
        else:
            entropy = 0

        return entropy

    @staticmethod
    def get_entropys(seqs):
        """
        get comentropy
        :param seqs:
        :return:
        """
        # check input
        if not seqs:
            logging.error('Please provide valid input sequences')
            raise ValueError('Please check input sequences')
        start_time = time.time()
        logging.info('=== Calculate sequence entropy ===')
        entropys = []
        for seq in seqs.values():
            en = SequenceInformation.topological_entropy(seq)
            entropys.append(en)
        end_time = time.time()
        timediff = end_time - start_time
        logging.info('Calculate sequence entropy done')
        logging.info("Time taken:" + str(timedelta(seconds=timediff)))

        return entropys
