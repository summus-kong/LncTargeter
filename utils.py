#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.9
@Author: Weibo Hou
@Description: utils
"""

import datetime
import pathlib
import subprocess
import logging

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def check_files_exist(*files):
    """
    Check if files exist
    :param files: a list files to check
    :return: return true only if all files exist
    """

    for f in files:
        if not f:
            return False
        elif not pathlib.Path(f).is_file():
            return False

    return True


def check_paths_exist(*paths):
    """
    Check if paths exist
    :param paths: a list paths to check
    :return: return true only if all files exist
    """

    for p in paths:
        if not p:
            return False
        elif not pathlib.Path(p).is_dir():
            return False

    return True


def check_bowtieindex(index):
    """
    Check if Bowtie index exist
    :param index: path to the index
    :return: return TRUE if Bowtie index exist
    """

    if not index:
        return False
    bidx = ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']
    index = [index + x for x in bidx]
    if check_files_exist(*index):
        return True


def get_timestamp(shorten=False):
    """
    Return current timestamp.
    :param shorten: return short version without space, dash and colons
    :return: timestamp as string
    """

    timestamp = str(datetime.datetime.now()).split(".")[0].replace(" ", "-")
    if shorten:
        timestamp = timestamp.replace("-", "").replace(" ", "").replace(":", "")
    return timestamp


def get_file_dir(file):
    """
    Get the folder of input file
    :param file: input file path
    :return: folder
    """

    return pathlib.Path(file).parent


def get_filename(file):
    """
    Return filename
    :param file: input file path
    :return: filename
    """

    return pathlib.Path(file).name


def mkdir(dir_path):
    """
    Create a directory
    :param dir_path: input folder path
    :return: true is folder created
    """
    try:
        pathlib.Path(dir_path).mkdir(parents=True)
    except OSError:
        return False
    return True


def covert_cmd(input_cmd):
    """
    Covert command list to string
    :param input_cmd:
    :return:
    """
    if input_cmd:
        if isinstance(input_cmd, list):
            return ' '.join(input_cmd)

    return input_cmd


def parse_biotools_options(*args, **kw):
    """
    Parse Bioinformatics toolkits options arguments
    example: parse_biotools_options(*('-1', '-2'), **{'--threads': 10})
    :param args: tuple positional arguments
    :param kw: dict keywords arguments
    :return: list
    """
    # parse
    pa = [str(x) for x in args]
    ka = []
    for i in kw.keys():
        ka.append(str(i))
        ka.append(str(kw[i]))

    # combine
    internal = []
    if pa:
        internal.extend(pa)
    if ka:
        internal.extend(ka)

    return internal


def run_biotools(command):
    """
    Run shell command by subprocess
    :param command: string shell command
    :return: CalledCompleted object
    """

    if not command:
        raise ValueError('Please check input arguments')

    res = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if res.stdout:
        res.stdout = res.stdout.decode('utf-8')
    else:
        res.stdout = ''
    if res.stderr:
        res.stderr = res.stderr.decode('utf-8')
    else:
        res.stdout = ''

    return res


def print_cc(res):
    """
    Print CalledCompleted object stderr and stdout
    :param res: CalledCompleted object
    :return: none
    """
    if not isinstance(res, subprocess.CompletedProcess):
        logging.error('input is not subprocess.CompletedProcess object')
        raise ValueError('Please check input')

    if isinstance(res, subprocess.CompletedProcess):
        if res.stdout:
            print("STDOUT:\n" + res.stdout)
        if res.stderr:
            print("STDERR:\n" + res.stderr)


def rmdir(dir_path):
    """
    Remove a directory
    :param dir_path: input directory path
    :return: true is folder created
    """
    cmd = ['rm', '-r', dir_path]
    try:
        rmres = run_biotools(covert_cmd(cmd))
        exitcode = rmres.returncode
        if exitcode == 0:
            return True
    except:
        return False


def get_version():
    """
    Get LncTargeter software version
    :return: version
    """
    version = 'LncTargeter version 0.1.0 \n '

    return version
