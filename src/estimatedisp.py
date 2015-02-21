#!/usr/bin/env python 
"""
Three steps to estimate dispersion.
""" 

import creatematrix as cm
import rawdisp as rd
import fitdisp as fd
import adjdisp as ad
import cPickle as pickle

def estimate_disp(data, opts):
    """ Create explanatory matrix and estimate dispersion. 
        Temporarily save data in ./TmpData.pkl file.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input arguments to the main TE function 
    @type opts: Instance
    """

    explanatory = cm.create_matrix(data, model='H1')
    data.matrix = explanatory

    outpath = opts.resPath
    pklFile = outpath + 'TmpData.pkl'

    if opts.dispDiff:
        data = rd.disper_raw(data, opts)
    else:
        data = rd.disper_raw_scalar(data, opts)
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

    print '*'*25

    data = fd.disper_fit(data, opts)
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

    print '*'*25

    if opts.dispDiff:
        data = ad.disper_adj(data, opts)
    else:
        data = ad.disper_adj_scalar(data, opts)
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

    return data
