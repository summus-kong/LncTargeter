#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.9
@Author: Weibo Hou
@Description: Run bioinformatics tools using python
"""

import os
import time
import sys
from datetime import timedelta
import subprocess
import logging
import pandas as pd
from src import utils as ut

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def rnafold(*args, infile, outfile, **kw):
    """
    Run RNAfold from ViennaRNA
    :param infile: input fasta file
    :param outfile: output RNAfold file
    :return: return the state of RNAfold
    """
    logging.info('=== 1. Run RNAfold to predict RNA secondary structure ===')
    # check input fasta file
    if not ut.check_files_exist(infile):
        logging.error('Please provide a valid input fasta file to build RNA secondary structure')
        raise FileNotFoundError
    # check output path
    out_path = ut.get_file_dir(outfile)
    if not ut.check_paths_exist(out_path):
        logging.warning('Output path is not exist and is being created')
        if not ut.mkdir(out_path):
            raise OSError('Failed to create output directory')
    options = ut.parse_biotools_options(*args, **kw)
    cmd = ['RNAfold', '<', infile, '--noPS']
    cmd.extend(options)
    cmd_suffix = ['>', outfile]
    cmd.extend(cmd_suffix)

    time_start = time.time()
    logging.info('Execute commands for predicting RNA secondary structure')
    try:
        # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
        res = ut.run_biotools(ut.covert_cmd(cmd))
        end_time = time.time()
        timediff = round(end_time - time_start)

        ut.print_cc(res)
        logging.info('Process RNAfold done')
        logging.info("Time taken of execute RNAfold commands:" + str(timedelta(seconds=timediff)))

        exitcode = res.returncode
        if exitcode == 0:
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


def build_index(*args, index=None, genome=None, verbose=False, **kw):
    """
    Build bowtie index
    :param verbose: output stdout and stderror
    :param *args: tuple positional arguements
    :param **kw: dict keyword arguments
    :param index: path where the index will be created
    :param genome: input genome fasta file
    :return: return the state of bowtie-build
    """

    logging.info('=== 3. Checking and build Bowtie index files ===')
    # if index already exist, then exit
    if index:
        if ut.check_bowtieindex(index):
            logging.info('Bowtie index already exist')
            logging.info(f'Bowtie index path: {index}')
            return True

    # check input genome file
    if not genome or not ut.check_files_exist(genome):
        logging.error('Please provide a valid input genome fasta file to build Bowtie index')
        raise FileNotFoundError("Please check input to bowtie build index")

    index_dir = ut.get_file_dir(index)
    # create index folder
    if not ut.check_paths_exist(index_dir):
        if not ut.mkdir(index_dir):
            raise OSError('Error creating bowtie index. Failed to create index directory')

    # parse arguments
    outlog = ['>', str(index_dir / 'bowtie_build.log'), '2>&1']
    # arguments list
    internal_arguments = ut.parse_biotools_options(*args, **kw)
    # command
    cmd = ['bowtie-build']
    cmd.extend(internal_arguments)
    cmd.append(genome)
    cmd.append(index)
    cmd.extend(outlog)
    # get current time
    time_start = time.time()
    logging.info('Execute Bowtie build commands')
    logging.info(f'Bowtie index path: {index_dir}')
    logging.info(f'Bowtie build log file path: {outlog[1]}')
    try:
        # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
        res = ut.run_biotools(ut.covert_cmd(cmd))
        end_time = time.time()
        timediff = round(end_time - time_start)

        if verbose:
            ut.print_cc(res)
        logging.info('Build bowtie index done')
        logging.info("Time taken of bowtie-build commands:" + str(timedelta(seconds=timediff)))

        exitcode = res.returncode
        if exitcode == 0:
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


def bowtie_align(*args, index_path, fs_seq, out_bwt, **kw):
    """
    Run Bowtie align program
    :param out_bwt: bowtie align result output path
    :param args: tuple positional arguments
    :param index_path: reference genome index path
    :param fs_seq: Sub Sequence fasta file
    :param kw: dict keywords arguments
    :return: command running status
    """
    logging.info('=== 4. Align subsequences to reference genome ===')
    # check index
    if not ut.check_bowtieindex(index_path):
        logging.error('Bowtie index is not exist and exit')
        raise FileNotFoundError("Please check input to bowtie align index")
    # check fs_sequence
    if not ut.check_files_exist(fs_seq):
        logging.error('Please provide a valid input Sub Sequence fasta file to Bowtie align program')
        raise FileNotFoundError("Please check input fasta to bowtie align")
    if not ut.check_paths_exist(ut.get_file_dir(out_bwt)):
        if not ut.mkdir(ut.get_file_dir(out_bwt)):
            raise OSError(f'Failed to create {ut.get_file_dir(out_bwt)} directory')

    flag = 1
    if flag:
        start_time = time.time()
        logging.info('Execute Bowtie align commands')
        options = ut.parse_biotools_options(*args, **kw)
        try:
            out_path = ut.get_file_dir(out_bwt)
            if not ut.check_paths_exist(out_path):
                logging.info(f'Creating output path: {out_path}')
                if ut.mkdir(out_path):
                    logging.info('Created path done')

            outlog = ['>', str(out_path / 'bowtie_build.log'), '2>&1']
            internal = ['--al', str(out_path / 'Reads_aligned'),
                        '--un', str(out_path / 'Reads_unaligned'),
                        index_path, fs_seq, out_bwt]
            cmd = ['bowtie', '-f', '--nofw']
            cmd.extend(internal)
            cmd.extend(options)
            cmd.extend(outlog)

            # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
            res = ut.run_biotools(ut.covert_cmd(cmd))
            end_time = time.time()
            timediff = round(end_time - start_time)

            ut.print_cc(res)
            logging.info('Execute bowtie align commands done')
            logging.info("Time taken of bowtie commands:" + str(timedelta(seconds=timediff)))

            exitcode = res.returncode
            if exitcode == 0:
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


def getfasta(*args, fi, bwt, fo, inpath, distince=200, **kw):
    """
    Extract fasta sequence using bedtools
    :param inpath: tem directory
    :param distince: distance upstream and downstream of the potential lncRNA binding regions
    :param args: tuple positional arguments
    :param fi: reference fasta file
    :param bwt: input bowtie aligned file
    :param fo: output fasta file
    :param kw: dict keywords arguments
    :return: logical
    """
    # check input files
    if not ut.check_files_exist(fi) or not ut.check_files_exist(bwt):
        logging.error('Please provide a valid input')
        raise FileNotFoundError("Please check input reference fasta file to bedtools getfasta")
    if not ut.check_paths_exist(inpath):
        if not ut.mkdir(inpath):
            raise OSError(f'Failed to create {inpath} directory')
    if not ut.check_paths_exist(ut.get_file_dir(fo)):
        if not ut.mkdir(ut.get_file_dir(fo)):
            raise OSError(f'Failed to create {ut.get_file_dir(fo)} directory')

    flag = 1
    if flag:
        logging.info('Extract putative sequences and character of the DNA binding region')
        options = ut.parse_biotools_options(*args, **kw)
        try:
            out_path = ut.get_file_dir(fo)
            if not ut.check_paths_exist(out_path):
                logging.info(f'Creating output path: {out_path}')
                if ut.mkdir(out_path):
                    logging.info('Created path done')
            awk = """awk 'BEGIN{OFS="\t"} """
            awk = awk + """{if($4-dis<0) print $3,0,$4+length($4)+dis,$3"_"$4"_"$4+length($5);else print $3,$4-dis,$4+length($4)+dis,$3"_"$4"_"$4+length($5)}' """
            txt = f'{inpath}/tem_target_pos.bed'
            awk = f"{awk} dis={distince} {bwt} > {txt}"
            awk = str(awk)
            logging.info('Prepare temporary bed file for getfasta')
            logging.info(f'Tem bed file path: {txt}')
            # logging.info(f'Codes: {repr(awk)}')
            awk_run = ut.run_biotools(awk)
            fa = inpath + '/tem_target_pos.fa'
            if not ut.check_files_exist(txt):
                logging.error('Please check getfasta awk server')
                return False
            if ut.check_files_exist(txt):
                internal = ['-fi', fi,
                            '-bed', txt,
                            '-fo', fa]
                cmd = ['bedtools', 'getfasta', '-nameOnly', '-tab']
                cmd.extend(internal)
                cmd.extend(options)

                logging.info('Prepare tab-separated binding regions fasta file')
                logging.info(f'Tab-separated file path: {fa}')
                # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
                res = ut.run_biotools(ut.covert_cmd(cmd))

                ut.print_cc(res)
                logging.info('Process DNA binding regions done')

                if ut.check_files_exist(fa):
                    logging.info('Convert and process the previous file')
                    parse = 'cat '
                    parse = parse + fa
                    parse = parse + """ | tr -s '_' '\t'"""
                    parse = parse + ' > ' + fo
                    parse = str(parse)
                    logging.info('Clean and covert tab-separated binding regions fasta file')
                    logging.info(f'Cleaned file path: {fo}')
                    # logging.info(f'Codes: {repr(parse)}')
                    run_parse = ut.run_biotools(parse)

                    exitcode = run_parse.returncode
                    if exitcode == 0:
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


def faidx(*args, ref, **kw):
    """
    Build fasta index .fai
    :param args:
    :param ref: reference fasta file
    :param kw:
    :return: state
    """
    logging.info('Index regions from reference genome fasta file')
    # check input files
    if not ut.check_files_exist(ref):
        logging.error('Please provide a valid input')
        raise FileNotFoundError

    # check .fai file
    fai = ref + '.fai'
    flag = 0
    if not ut.check_files_exist(fai):
        flag = 1
    if ut.check_files_exist(fai):
        logging.info('Reference genome file .fai index already exist')
        return True
    if flag:
        logging.info('Build reference genome fasta file index .fai')
        options = ut.parse_biotools_options(*args, **kw)
        cmd = ['samtools', 'faidx', ref]
        cmd.extend(options)

        try:
            # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
            res = ut.run_biotools(ut.covert_cmd(cmd))

            ut.print_cc(res)
            logging.info('Build fasta index done')

            exitcode = res.returncode
            if exitcode == 0:
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


def merge(*args, inbed, out, tem_path, d=0, **kw):
    """
    Merge putative LncRNA target regions with scores greater than 0.5 and overlap greater than maximum overlap distance
    :param d: Maximum distance between features allowed for features to be merged. Default is 0. That is, overlapping and/or book-ended features are merged
    :param tem_path:
    :param out:
    :param args:
    :param inbed: LncRNA target regions bed file
    :param kw:
    :return:
    """
    print('Additional: ')
    logging.info('Merge overlapped putative LncRNA DNA target regions for removing redundant')

    # check input
    if not ut.check_files_exist(inbed):
        logging.error('Please check prediction server')
        raise FileNotFoundError
    if not ut.check_paths_exist(ut.get_file_dir(out)):
        if not ut.mkdir(ut.get_file_dir(out)):
            raise OSError(f'Failed to create {ut.get_file_dir(out)} directory')
    if not ut.check_paths_exist(tem_path):
        if not ut.mkdir(tem_path):
            raise OSError(f'Failed to create {tem_path} directory')

    options = ut.parse_biotools_options(*args, **kw)
    tem_bed = f'{tem_path}/merge_temsite.bed'
    cmd = ['sort -k1,1 -k2,2n', inbed, '|', 'bedtools', 'merge', '-c', '4,5', '-o', 'distinct,median', '-d', str(d)]
    cmd.extend(options)
    cmd.append(f' > {tem_bed}')
    try:
        # logging.info(f'Codes: {ut.covert_cmd(cmd)}')
        res = ut.run_biotools(ut.covert_cmd(cmd))

        # ut.print_cc(res)
        # logging.info('Build fasta index done')

        exitcode = res.returncode
        if exitcode == 0:
            pro_out = pd.read_csv(tem_bed, names=['chr', 'start', 'end', 'core sequence', 'score'], sep='\t')
            pro_out.to_csv(out, sep='\t', index=False)
            if ut.check_files_exist(out):
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


def intersect_gene(intar, tempath, tss_bed,
                   target_up=50000, target_down=50000,
                   min_score=0.5, tss_up=500, tss_down=500):
    """
    Extract genes that may be regulated by lncRNA near the DNA target regions
    :param tss_down: the distance of gene TSS downstream
    :param tss_up: the distance of gene TSS upstream
    :param tss_bed: gene TSS bed file
    :param tempath: temporary path
    :param min_score: DNA target regions score
    :param intar: DNA target regions file
    :param target_up: the distance of DNA target regions upstream
    :param target_down: the distance of DNA target regions downstream
    :return: run state
    """
    # check input
    if not ut.check_files_exist(intar) or not ut.check_files_exist(tss_bed):
        logging.error('Please check input for intersect_gene function')
        raise FileNotFoundError
    if not ut.check_paths_exist(tempath):
        if not ut.mkdir(tempath):
            raise OSError(f'Failed to create {tempath} directory')

    awk_region = "awk" + r" '" + """BEGIN{OFS="\t"}""" + f" $NF>{min_score}&&NR>1"
    awk_region = awk_region + """{if(($2-up)>0) print "chr"$1,$2-up,$3+down,$(NF-2),$(NF-1),"chr"$1,$2,$3;else print "chr"$1,0,$3+down,$(NF-2),$(NF-1),"chr"$1,$2,$3}'"""
    regions = f'{tempath}/tem_regions_for_gene_predict.bed'
    awk_region = f'{awk_region} up={target_up} down={target_down} {intar} > {regions}'
    awk_region = str(awk_region)
    logging.info('Prepare temporary bed file for target genes')
    logging.info(f'Tem bed file path: {regions}')
    # logging.info(f'Codes: {repr(awk_region)}')
    awk_region_run = ut.run_biotools(awk_region)

    tem_pro = f'{tempath}/tem_process.txt'
    if not ut.check_files_exist(regions):
        logging.error('Please check intersect_gene awk server')
        return False
    if ut.check_files_exist(regions):
        cmd = f'bedtools intersect -a {regions} -b {tss_bed} -wo | '
        awk_e = "awk" + r" '" + """BEGIN{OFS="\t"}"""
        awk_e = awk_e + "{if($7>$11) print $0,$6,$11,$7;if($8<$10) print $0,$6,$8,$10;if($7<$11&&$11<$8);print $0,$1,$2,$11}' "
        cmd = f"{cmd}{awk_e}> {tem_pro}"
        cmd = str(cmd)
        logging.info(f'Tem process file path: {tem_pro}')
        # logging.info(f'Codes: {repr(cmd)}')
        res = ut.run_biotools(cmd)

    getfasta_bed = f"{tempath}/for_getfasta_linker_sequence.bed"
    if not ut.check_files_exist(tem_pro):
        logging.error('Please check intersect_gene bedtools intersect server')
        return False
    if ut.check_files_exist(tem_pro):
        awk_linker = "awk" + r" '" + """BEGIN{OFS="\t"}"""
        awk_linker = awk_linker + "($(NF-1)<$NF){print $(NF-2),$(NF-1),$NF}' " + f"{tem_pro} | sed 's/chr//g' | sort | uniq "
        awk_linker = f"{awk_linker}> {getfasta_bed}"
        awk_linker = str(awk_linker)
        logging.info(f'Tem linker bed file path: {getfasta_bed}')
        # logging.info(f'Codes: {repr(awk_linker)}')
        awk_linker_run = ut.run_biotools(awk_linker)

    tss_getfa_bed = f"{tempath}/for_getfasta_tss_region.bed"
    if not ut.check_files_exist(getfasta_bed):
        logging.error('Please check intersect_gene awk server')
        return False
    if ut.check_files_exist(getfasta_bed):
        awk_tss = "awk" + r" '" + """BEGIN{OFS="\t"}"""
        awk_tss = awk_tss + "{print $9,$10-up,$11+down,$12}' " + f" up={tss_up} down={tss_down} {tem_pro} | sed 's/chr//g' | sort | uniq "
        awk_tss = f"{awk_tss}> {tss_getfa_bed}"
        awk_tss = str(awk_tss)
        logging.info(f'Tem tss bed file path: {tss_getfa_bed}')
        # logging.info(f'Codes: {repr(awk_tss)}')
        awk_tss_run = ut.run_biotools(awk_tss)

    if not ut.check_files_exist(tss_getfa_bed):
        logging.error('Please check intersect_gene awk server')
        return False

    return True


def getfastaTgene(fi, tempath):
    """
    Extract gene TSS region and linker sequence
    :param fi: reference genome file
    :param tempath: temporary path
    :return: run state
    """
    # check input
    if not ut.check_files_exist(fi) or not ut.check_files_exist(fi):
        logging.error('Please check input for intersect_gene function')
        raise FileNotFoundError
    if not ut.check_paths_exist(tempath):
        if not ut.mkdir(tempath):
            raise OSError(f'Failed to create {tempath} directory')

    bed_linker = f"{tempath}/for_getfasta_linker_sequence.bed"
    out_linker = f"{tempath}/linker_sequence_tab.txt"
    cmd_linker = f'bedtools getfasta -fi {fi} -bed {bed_linker} -bedOut | '
    cmd_linker = cmd_linker + "awk" + r" '" + """{print "chr"$0"\t"length($4)}' """
    cmd_linker = f"{cmd_linker}> {out_linker}"
    cmd_linker = str(cmd_linker)
    logging.info(f'Tem linker sequence file path: {out_linker}')
    # logging.info(f'Codes: {repr(cmd_linker)}')
    linker_run = ut.run_biotools(cmd_linker)

    bed_tss = f"{tempath}/for_getfasta_tss_region.bed"
    out_tss = f"{tempath}/tss_sequence_tab.txt"
    if not ut.check_files_exist(out_linker):
        logging.info('Please check getfastaTgene bedtools server')
        return False
    if ut.check_files_exist(out_linker):
        cmd_tss = f'bedtools getfasta -fi {fi} -bed {bed_tss} -bedOut | '
        cmd_tss = cmd_tss + "awk" + r" '" + """{print "chr"$0}' """ + f'> {out_tss}'
        cmd_tss = str(cmd_tss)
        logging.info(f'Tem TSS sequence file path: {out_tss}')
        # logging.info(f'Codes: {repr(cmd_tss)}')
        linker_run = ut.run_biotools(cmd_tss)
    if not ut.check_files_exist(out_tss):
        logging.info('Please check getfastaTgene bedtools server')
        return False

    return True
