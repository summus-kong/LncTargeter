#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.14
@Author: Weibo Hou
@Description: Parse Argument
"""

import argparse
from src import utils as ut


def parseArguments():
    """
    Parse command-line arguments passed to the pipeline.
    :return:
    """
    parser = argparse.ArgumentParser(prog='LncTarget',
                                     add_help=False,
                                     description='Predicting LncRNA binding sites and their potential regulatory genes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required Arguments',
                                         'Parameters must be supplied, otherwise throw an exception.')
    required.add_argument('-i', '--input_fasta',
                          type=str,
                          metavar='file',
                          required=True,
                          help='LncRNA sequence file used to predict DNA binding regions(Fasta file format).')
    required.add_argument('-r', '--reference',
                          type=str,
                          metavar='file',
                          required=True,
                          help='Reference genome file, the default file format is Fasta.')
    required.add_argument('-b', '--bt_index_base',
                          type=str,
                          metavar='basename',
                          required=True,
                          help='If bowtie index is not exist, write bt data to files with this dir/basename.')
    required.add_argument('-o', '--output',
                          type=str,
                          metavar='file',
                          required=True,
                          help='The output potential binding regions file, which is tab-separated. The file is '
                               'revised to include the chromosome location information of the binding sites, '
                               'feature information and the probability that the site is an lncRNA binding site.')
    optional = parser.add_argument_group('Optional Arguments',
                                         'Specify additional non-essential parameters.')
    optional.add_argument('-g', '--genome',
                          type=str,
                          metavar='<str>',
                          required=False,
                          default='mm10',
                          choices=['mm10', 'hg38'],
                          help='specify reference genome')
    optional.add_argument('-p', '--threads',
                          type=int,
                          metavar='<int>',
                          required=False,
                          default=1,
                          help='# of threads')
    optional.add_argument('-c', '--min_complexity',
                          required=False,
                          type=float,
                          metavar='<float>',
                          default=0.5,
                          help='Minimum complexity value for screening Sub sequence from input lncRNA fasta sequence. '
                               'Through many experiments, it is found that when the value is equal to 0.5, the model '
                               'has a good effect, so it is recommended to use the default value.')
    optional.add_argument('-l', '--min_length',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=17,
                          help='Minimum sequence length for extract lncRNA fasta sub sequence. For better performance, '
                               'we recommend using the default value')
    optional.add_argument('-m', '--max_length',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=50,
                          help='Maximum sequence length for extract lncRNA fasta sub sequence. It is found that when '
                               'the value is equal to 50, the performance of the model is more impressive.')
    optional.add_argument('-d', '--extend_distance',
                          required=False,
                          default=200,
                          type=int,
                          metavar='<int>',
                          help='Distance upstream and downstream of the potential lncRNA binding regions. In order to '
                               'better analyze the LncRNA DNA binding regions, expand the upstream and downstream of '
                               'the start and end of the regions as its potential binding regions for further analysis.')
    optional.add_argument('-k',
                          type=int,
                          metavar='<int>',
                          required=False,
                          default=1,
                          help='report up to <int> good alignments per read')
    optional.add_argument('--seed',
                          type=int,
                          metavar='<int>',
                          required=False,
                          default=1,
                          help='Seed for random number generator.')
    optional.add_argument('--mismatch',
                          type=int,
                          metavar='<int>',
                          required=False,
                          default=3,
                          help='Report end-to-end hits w/ <=v mismatches.')
    optional.add_argument('--overlap',
                          type=int,
                          metavar='<int>',
                          required=False,
                          default=0,
                          help='Maximum distance between features allowed for features to be merged. Default is 0. '
                               'That is, overlapping and/or book-ended features are merged.')
    optional.add_argument('--target_upstream',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=50000,
                          help="Distance upstream of the potential lncRNA binding regions for prediction of target genes.")
    optional.add_argument('--target_downstream',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=50000,
                          help="Distance downstream of the potential lncRNA binding regions for prediction of target genes.")
    optional.add_argument('--tss_upstream',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=500,
                          help="Distance upstream of the gene TSS.")
    optional.add_argument('--tss_downstream',
                          required=False,
                          type=int,
                          metavar='<int>',
                          default=500,
                          help="Distance downstream of the gene TSS.")
    optional.add_argument('--min_score',
                          required=False,
                          type=float,
                          metavar='<float>',
                          default=0.5,
                          help="Minimum score for screening lncRNA DNA target regions.")
    optional.add_argument('--gene_score',
                          required=False,
                          type=float,
                          metavar='<float>',
                          default=0.5,
                          help='Minimum score for target genes.')
    optional.add_argument('--core_threshold',
                          required=False,
                          type=float,
                          metavar='<float>',
                          default=0.5,
                          help='Threshold for screening core sequence')
    optional.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit.")
    optional.add_argument('-v', '--version',
                          action='version',
                          version=ut.get_version(),
                          help='Print version information and exit.')

    args = parser.parse_args()

    return args
