#!usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Date: 2021.11.30
@Author: Weibo Hou
@Description: Filter binding sites
"""
import logging
import sys
import subprocess
import pandas as pd
from src import utils as ut

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def search_peak(bases_counts, threshold, rnafold):
    """
    Search peaks and return a list contains element subscript
    :param rnafold: opened rnafold file
    :param threshold: line y = threshold
    :param bases_counts: base counts
    :return: a list contains element subscript
    """
    peak_pos = []
    inter_pos = []
    for i in bases_counts.index:
        if bases_counts.index[0] < i < bases_counts.index[-1]:
            if bases_counts[i - 1] < bases_counts[i] > bases_counts[i + 1]:
                peak_pos.append(i)
            if bases_counts[i] <= threshold <= bases_counts[i + 1] or bases_counts[i + 1] <= threshold <= bases_counts[i]:
                inter_pos.append(i)
    # print(inter_pos)
    out_start = []
    out_end = []
    star_seq = []
    for index, content in enumerate(inter_pos):
        if index < len(inter_pos) - 1:
            seq_pos = [inter_pos[index], inter_pos[index + 1]]
            for i in peak_pos:
                if seq_pos[0] <= i <= seq_pos[1]:
                    out_start.append(seq_pos[0])
                    out_end.append(seq_pos[1])
                    for info in rnafold:
                        sequence = info[1].upper()
                        star_seq.append(sequence[seq_pos[0]: seq_pos[1]])

    out_dict = {'chr': 'chr1', 'start': out_start, 'end': out_end, 'core_sequence': star_seq}

    return out_dict


def identity_star_seq(total_merge, outbed, rnafold, cut_value):
    """
    Screen feature sequence
    :param rnafold: opened rnafold file
    :param outbed: output bed file
    :param cut_value: value for screen feature sequence
    :param total_merge: data frame
    :return: run state
    """
    cols = ['start_fs', 'end_fs']
    lnc_sites = total_merge[cols]
    bases = []
    start = lnc_sites['start_fs'].values
    end = lnc_sites['end_fs'].values
    for index, content in enumerate(start):
        base = range(content, end[index])
        bases.extend(base)

    # sbc single_base_counts
    sbc = pd.value_counts(bases).sort_index()
    idx = sbc.index
    # scale to [0, 1]
    counts_max = max(sbc)
    counts_min = min(sbc)
    scale_counts = pd.Series([(i - counts_min) / (counts_max - counts_min) for i in sbc], index=idx)
    seqs_out = search_peak(scale_counts, threshold=cut_value, rnafold=rnafold)
    pd.DataFrame(seqs_out).drop_duplicates().to_csv(outbed, sep='\t', header=False, index=False)
    if not ut.check_files_exist(outbed):
        return False

    return True


def extract_sites(feature_sequence, result_xgb, tem_path):
    """
    Extract sequence related binding sites
    :param tem_path: output tem path
    :param feature_sequence: sequence bed file
    :param result_xgb: run xgboost result
    :return: run state
    """
    # check input
    if not ut.check_files_exist(feature_sequence) or not ut.check_paths_exist(tem_path):
        logging.error('Please provide a valid input')
        raise ValueError('Please check input reference fasta file to extract sites')

    usecols = ['start_fs', 'end_fs', 'chr', 'start', 'end', 'score']
    result_xgb = result_xgb[usecols]
    result_xgb.insert(0, 'chr_lnc', 'chr1')
    tem_sites = f'{tem_path}/predicted_sites_unfiltered.bed'
    result_xgb.to_csv(tem_sites, sep='\t', header=False, index=False)
    tem_bed_res = f'{tem_path}/bedtools_feature_sequence.bed'
    cmd = f'bedtools intersect -a {tem_sites} -b {feature_sequence} -wa -wb > {tem_bed_res}'
    try:
        # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
        res = ut.run_biotools(ut.covert_cmd(cmd))
        ut.print_cc(res)
        exitcode = res.returncode
        # print(exitcode)
        if exitcode == 0:
            colnames = ['chr_lnc', 'chr_fe', 'start_fe', 'end_fe', 'core_sequence']
            colnames[1:1] = usecols
            presites = pd.read_csv(tem_bed_res, names=colnames, sep='\t')
            # print(presites.head())
            presites = presites[['chr', 'start', 'end', 'core_sequence', 'score']]
            unmerge_sites = f'{tem_path}/unmerge_binding_sites.bed'
            presites.to_csv(unmerge_sites, sep='\t', index=False, header=False)
            return True
        else:
            logging.info(f"Following error occurred executing above command (return code={exitcode}):")
            print("STDOUT:\n" + res.stdout)
            print("STDERR:\n" + res.stderr)
            return False
    except subprocess.CalledProcessError as e:
        logging.error("CalledProcessError exception occurred.\n" + str(e))
        return False
    except:
        logging.error("Fatal error occurred during execution.\n" + str(sys.exc_info()[0]))
        return False
