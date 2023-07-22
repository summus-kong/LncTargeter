#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.14
@Author: Weibo Hou
@Description: Predict lncRNA binding site
"""

import os
import pathlib
import pandas as pd
import src
import xgboost as xgb
from src import parse_presite as pp
from src import parse_rnafold as pr
from src import utils as ut
from src import run_biotools as rb
from src import parse_bowtie_result as pbr
from src import xgb_model as model
from src import process_tab_seq as pts
import time
from datetime import timedelta
import logging
from src import parse_args as pa
from src import filter_sites as fl


def mkdir_tem(path='./'):
    """
    Create temporary folder for storing temporary files
    :param path: temporary folder path
    :return: logical
    """
    tem_path = path + 'tem'
    if not ut.check_paths_exist(tem_path):
        logging.warning('Temporary path is not exist and creats it now')
        if not ut.mkdir(tem_path):
            raise OSError('Failed to create Temporary folder')

    if not ut.check_files_exist(tem_path):
        return False

    return True


def out_rnafold_path(inpath):
    """
    Get RNAfold output file path
    :param inpath: temporary file path
    :return: path
    """
    inpath = pathlib.Path(inpath)
    outfile = 'rnafold_result.str'
    res = str(inpath / outfile)

    return res


def out_subseq_path(inpath):
    """
    Get sub sequence temporary file path
    :param inpath: temporary file path
    :return: path
    """
    inpath = pathlib.Path(inpath)
    outfa = 'subsequence.fa'
    outcsv = 'subsequence_information.txt'
    fa_res = str(inpath / outfa)
    csv_res = str(inpath / outcsv)

    return fa_res, csv_res


def run_predict_bs():
    """
    Perform the prediction of lncRNA potential DNA binding sites and target genes
    :return: logical run state
    """
    args = pa.parseArguments()
    infa = args.input_fasta
    out = args.output
    ref = args.reference
    index = args.bt_index_base
    threads = args.threads
    mincm = args.min_complexity
    maxlen = args.max_length
    minlen = args.min_length
    myseed = args.seed
    v = args.mismatch
    k = args.k
    dis = args.extend_distance
    target_up_int = args.target_upstream
    target_down_int = args.target_downstream
    tss_up = args.tss_upstream
    tss_down = args.tss_downstream
    min_score = args.min_score
    gene_score = args.gene_score
    ct = args.core_threshold
    d = args.overlap
    genome = args.genome

    print(' ')
    print('Commands parameters: \n')
    options = vars(args)
    for key, value in options.items():
        print(f'    {key}: {value}')
    print(' ')

    outpath = ut.get_file_dir(out)
    outtem = str(outpath / 'tem')
    if not ut.check_paths_exist(outtem):
        ut.mkdir(outtem)

    logging.info('========== Prediction of putative LncRNA DNA target regions ==========')
    start_time = time.time()
    print(' ')
    # step 1
    # run ranfold
    rnafold_res = out_rnafold_path(outtem)
    # flag for rnafold state
    flag_rnafold = 0
    if rb.rnafold(infile=infa, outfile=rnafold_res):
        flag_rnafold = 1
    if not flag_rnafold:
        logging.error('Please check RNAfold server')
        return False

    print(' ')
    # step 2
    # extract sub sequence from lncRNA
    logging.info('=== 2. Extract sub sequence and calculate free energy ===')
    flag_csv = 0
    flag_fa = 0
    if not ut.check_files_exist(rnafold_res):
        logging.error('Please check RNAfold result')
        return False
    sub_res = out_subseq_path(outtem)
    sub_fa_res, sub_csv_res = sub_res
    rnafold = pr.openfile(rnafold_res)
    seqs = pr.getRRenergy(rnafold, mincm=mincm, minlength=minlen, maxlength=maxlen, myseed=myseed)
    pr.out_csv(seqs, sub_csv_res)
    pr.to_fa(seqs, sub_fa_res)
    if ut.check_files_exist(sub_fa_res) and ut.check_files_exist(sub_csv_res):
        flag_fa = 1
        flag_csv = 1

    if not flag_csv and not flag_fa:
        logging.error('Please check extracted sub sequence')
        return False

    print(' ')
    # step 3
    # run bowtie_build
    flag_bb = 0
    if flag_fa and flag_csv:
        if rb.build_index(index=index, genome=ref):
            flag_bb = 1
    if not flag_bb:
        logging.error('Please check Bowtie build server')
        return False

    print(' ')
    # step 4
    # run bowtie
    flag_ba = 0
    outbw_path = outtem + '/bowtie/lncRNA_bowtie_align.bwt'
    if flag_bb:
        extend_args = {'--threads': threads, '-v': v, '--seed': myseed, '-k': k}
        if rb.bowtie_align(**extend_args, index_path=index, fs_seq=sub_fa_res, out_bwt=outbw_path):
            flag_ba = 1
    if not flag_ba:
        logging.error('Please check Bowtie align server')
        return False

    print(' ')
    # step 5
    logging.info('=== 5. Prepare input data and execute XGBoost model ===')
    # fasta index
    flag_fi = 0
    if flag_ba:
        if rb.faidx(ref=ref):
            flag_fi = 1
    if not flag_fi:
        logging.error('Please check samtools faidx server')
        return False

    # extends DNA binding sites and get fasta sequence
    flag_gf = 0
    outfa_path = outtem + '/tem_target_pos_tab.fa'
    if flag_fi:
        if rb.getfasta(fi=ref, bwt=outbw_path, fo=outfa_path, inpath=outtem, distince=dis):
            flag_gf = 1
    if not flag_gf:
        logging.error('Please check extract potential lncRNA binging region fasta sequence server')
        return False

    # parse presite
    flag_pp = 0
    if flag_gf:
        site = pp.parser_site(outfa_path)
        pro_site = pp.calculate_deltaG(site)
        # read bowtie aligned file
        bwt = pbr.get_bwtResult(outbw_path)
        fs = pbr.get_fs_result(sub_csv_res)
        merge_bwt = pbr.merge_res(bwt, fs)

        del_cols = ['strand', 'sequence', 'quality', 'mapper_numbers', 'mismatch']
        merge_bwt.drop(del_cols, axis=1, inplace=True)
        reorder_col = ['ref_chr', 'first_position', 'end_position', 'query_id', 'sequence_fs', 'cm_value',
                       'RNA-RNA_deltaG', 'start_base', 'end_base', 'sequence_length', 'sub_GC']
        merge_bwt = merge_bwt[reorder_col]
        rename_col = ['chr', 'start', 'end', 'lncid', 'sequence_fs', 'cm_fs', 'RR_fs', 'start_fs', 'end_fs', 'length_fs', 'gc_fs']
        merge_bwt.columns = rename_col

        # produce data frame for xgboost input
        total_merge = pp.merge_site_fs(merge_bwt, pro_site)
        if not total_merge.empty:
            logging.info('Prepare input data for XGBoost model done')
            flag_pp = 1

        if not flag_pp:
            logging.error('Please check prepare XGBoost input file server')
            return False

    # run XGBoost
    # np.set_printoptions(suppress=True)
    # pd.set_option('display.float_format', lambda x: '%.6f' % x)
    flag_xgb = 0
    if flag_pp:
        use_cols = ['cm_fs', 'gc_fs', 'RR_fs', 'start_fs', 'length_fs', 'RD_site', 'cm_site', 'gc_site']
        input_xgb = xgb.DMatrix(total_merge[use_cols])
        site_model = src.get_model('LncRNABinder.json')
        xgb_prob = model.run_xgboost(input_xgb, site_model)
        total_merge['score'] = xgb_prob
        if not total_merge.empty:
            flag_xgb = 1

        if not flag_xgb:
            logging.error('Please check XGBoost server')
            return False

    # output predicted result
    flag_out = 0
    for_gene = f'{outtem}/for_gene_sites.bed'
    if flag_xgb:
        # output feature bed file
        bedpath = f'{outtem}/feature_sequence.bed'
        if not fl.identity_star_seq(total_merge, bedpath, rnafold=rnafold, cut_value=ct):
            logging.warning('Please check star sequence procedure')
            return False
        total_merge = total_merge[total_merge.score > min_score]
        total_merge = total_merge.loc[~total_merge['chr'].isin(['X', 'Y'])]
        total_merge.to_csv(for_gene, sep='\t', index=False)
        # print(total_merge.head())
        if fl.extract_sites(bedpath, total_merge, outtem):
            unmerged = f'{outtem}/unmerge_binding_sites.bed'
            if not rb.merge(inbed=unmerged, out=out, tem_path=outtem, d=d):
                logging.warning('Please check merge procedure')
        else:
            logging.info('Please check extract procedure')
            return False
        # cols_out = ['chr', 'start', 'end', 'lncid', 'scores']
    if ut.check_files_exist(out) and ut.check_files_exist(for_gene):
        flag_out = 1

    if not flag_out:
        logging.error('Please check output server')
        return False

    print(' ')
    logging.info('========== Prediction of LncRNA target genes ==========')
    # step 1
    logging.info('=== 1. Produce temporary files for prediction of target genes ===')
    flag_gene_pre = 0
    # Please prepare file pcg.bed (all protein coding genes)
    if flag_out:
        tsbed = src.get_mouse('pcg.bed')
        if genome == 'hg38':
            tsbed = src.get_human('pcg.bed')
        ing = rb.intersect_gene(for_gene, outtem, tsbed, min_score=min_score, target_up=target_up_int,
                                target_down=target_down_int, tss_up=tss_up, tss_down=tss_down)  # out需要改
        if ing:
            flag_gene_pre = 1
        if not ing:
            logging.error('Please check intersect_gene server')
            return False

    # get fasta
    flag_gene_fa = 0
    if flag_gene_pre:
        getfa = rb.getfastaTgene(ref, outtem)
        if getfa:
            flag_gene_fa = 1
        if not getfa:
            logging.error('Please check getfastaTgene server')
            return False

    # merge data
    flag_gene_merge = 0
    # read tem total
    if flag_gene_fa:
        regions = f'{outtem}/tem_process.txt'
        out_linker = f"{outtem}/linker_sequence_tab.txt"
        out_tss = f"{outtem}/tss_sequence_tab.txt"
        re_cols = ['chr_ex', 'start_ex', 'end_ex', 'cm_site', 'gc_site', 'chr_region', 'start_region',
                   'end_region', 'chr_gene', 'start_gene', 'end_gene', 'gene_symbol', 'gene_name', 'overlap',
                   'chr_linker', 'start_linker', 'end_linker']
        tem_total = pts.read_bed_file(regions, re_cols)
        linker_cols = ['chr_linker', 'start_linker', 'end_linker', 'sequence', 'linker_length']
        tss_cols = ['chr_gene', 'start_gene', 'end_gene', 'gene_symbol', 'sequence']
        gene_pos_cols = ['chr_gene', 'gene_body_start', 'gene_body_end', 'gene_symbol']
        chr_length_cols = ['chr_gene', 'chr_num', 'chr_length']
        linkers = pts.read_bed_file(out_linker, linker_cols)
        tss = pts.read_bed_file(out_tss, tss_cols)
        gene_bed = src.get_mouse('gene_pos.bed')
        chr_bed = src.get_mouse('chr_length')
        if genome == 'hg38':
            gene_bed = src.get_human('gene_pos.bed')
            chr_bed = src.get_human('chr_length')
        true_gene = pts.read_bed_file(gene_bed, gene_pos_cols)
        chr_length = pts.read_bed_file(chr_bed, chr_length_cols)
        # process linker and tss
        linkers_processed = pts.process_bed(linkers, ['linker_gc', 'linker_cm'])
        tss_processed = pts.process_bed(tss, ['tss_gc', 'tss_cm'])
        # merge
        m1 = pd.merge(tem_total, tss_processed, on=['gene_symbol', 'chr_gene'], how='inner')
        m2 = pd.merge(m1, linkers_processed, on=['chr_linker', 'start_linker', 'end_linker'], how='inner')
        m3 = pd.merge(m2, true_gene, on=['chr_gene', 'gene_symbol'], how='inner')
        m4 = pd.merge(m3, chr_length, on='chr_gene', how='inner')

        chr_number = max(m4['chr_num'])
        m4['chr_trans'] = m4['chr_num'] / chr_number
        m4['rs_trans'] = m4['start_region'] / m4['chr_length']
        m4['re_trans'] = m4['end_region'] / m4['chr_length']
        m4['gbs_trans'] = m4['gene_body_start'] / m4['chr_length']
        m4['gbe_trans'] = m4['gene_body_end'] / m4['chr_length']
        # print(m4.columns)
        # clean
        cols_further = ['chr_region', 'start_region', 'end_region', 'gene_symbol', 'gene_name', 'cm_site', 'gc_site',
                        'linker_length', 'linker_gc', 'linker_cm', 'tss_gc', 'tss_cm', 'chr_trans', 'rs_trans',
                        're_trans', 'gbs_trans', 'gbe_trans', 'chr_gene', 'gene_body_start', 'gene_body_end']
        m4 = m4[cols_further]
        if not m4.empty:
            flag_gene_merge = 1

        if not flag_gene_merge:
            logging.error('Please check merge_table server')
            return False

    # step 2
    logging.info('=== 2. Run XGBoost for prediction of target genes ===')
    flag_gene_xgb = 0
    cols_for_model = ['cm_site', 'gc_site', 'linker_length', 'linker_gc', 'linker_cm', 'tss_gc', 'tss_cm',
                      'chr_trans', 'rs_trans', 're_trans', 'gbs_trans', 'gbe_trans']
    input_xgb_gene = xgb.DMatrix(m4[cols_for_model])
    gene_path = src.get_model('LncRNABinder_Tgene_XGBoost.json')
    xgb_prob_gene = model.run_xgboost(input_xgb_gene, gene_path)
    m4['score'] = xgb_prob_gene
    if not m4.empty:
        flag_gene_xgb = 1

    if not flag_gene_xgb:
        logging.error('Please check XGBoost server')
        return False
    # output predicted result
    flag_out = 0
    out_gene_res = f"{outpath}/lncRNA_target_gene.txt"
    if flag_xgb:
        out_predicted_genes = m4[['chr_gene', 'gene_body_start', 'gene_body_end', 'gene_symbol', 'gene_name', 'score']]
        out_predicted_genes.columns = ['chr', 'start', 'end', 'gene_symbol', 'gene_name', 'score']
        out_predicted_genes = out_predicted_genes[out_predicted_genes.score > gene_score]
        out_predicted_genes = out_predicted_genes.groupby(['chr', 'start', 'end', 'gene_symbol', 'gene_name'])['score'].mean().reset_index()
        out_predicted_genes.to_csv(out_gene_res, sep='\t', index=False)
    if ut.check_files_exist(out_gene_res):
        flag_out = 1
        if ut.check_paths_exist(outtem):
            rmdir = f'rm -r {outtem}'
            ut.run_biotools(rmdir)

    if not flag_out:
        logging.error('Please check output server')
        return False

    end_time = time.time()
    timediff = round(end_time - start_time)
    logging.info('Process prediction of putative LncRNA DNA target regions and genes done')
    logging.info("Time taken of prediction of putative LncRNA DNA target regions and genes: " + str(timedelta(seconds=timediff)))
    logging.info('END')

    return True
