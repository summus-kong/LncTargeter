#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date: 2021.10.9
@Author: Weibo Hou
@Description: utils
"""

import pathlib
import logging
import xgboost as xgb
import numpy as np
from src import utils as ut

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s: %(asctime)s %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S"
                    )


def run_xgboost(inmatrix, xgb_model):
    """
    Run XGBoost model to predicted LncRNA DNA binding regions
    :param xgb_model: XGBoost model path
    :param inmatrix: DMatrix object
    :return: probability np.array
    """
    # check input
    if not inmatrix:
        logging.error('Please provide valid parameter')
        raise ValueError
    if not ut.check_files_exist(xgb_model):
        logging.error('XGBoost model not exist')
        raise FileNotFoundError('Please provide valid XGBoost model path')

    logging.info('Run XGBoost for predicting LncRNA DNA binding regions')
    lncrna_xgb = xgb.Booster()
    file = pathlib.Path(xgb_model)
    file = file.resolve()
    lncrna_xgb.load_model(file)
    xgb_prob = np.around(lncrna_xgb.predict(inmatrix), 5)

    if xgb_prob.size:
        logging.info('XGBoost done')
        return xgb_prob
    if not xgb_prob.size:
        logging.error('Please check XGBoost predicted model')
        return None
