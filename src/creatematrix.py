#!/usr/bin/env python
"""
Creating the explanatory matrix for GLM.
"""

import sys
import numpy as np

def create_matrix(data, model='H0'):

    width = len(data.exper[:, 0])
    if model == 'H0':
        explanatory = np.zeros((width, 3))
        explanatory[data.idxCtl, 0] = 1
        explanatory[data.idxTrt, 1] = 1
        explanatory[data.idxRna, 2] = 1
    elif model == 'H1':
        explanatory = np.zeros((width, 4))
        explanatory[data.idxCtl, 0] = 1
        explanatory[data.idxTrt, 1] = 1
        explanatory[data.idxRna, 2] = 1
        explanatory[np.intersect1d(data.idxRibo, data.idxTrt), 3] = 1
    else:
        sys.stderr.write('ERROR: The parameter \'model\' in create_sparse() can only accept \'H0\' or \'H1\' as its input.\n')

    return explanatory
