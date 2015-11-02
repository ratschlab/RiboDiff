#!/usr/bin/env python 
"""
Statistical test.
"""

import sys
import numpy as np
import creatematrix as cm
import statsmodels.api as sm
from scipy.stats import chi2
import statsmodels.sandbox as sms

def test_count(data, opts):
    """
    Make a test for all genes iteratively.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input argument to the main TE function 
    @type opts: Instance
    """

    print 'Start the statistical test.'

    num = len(data.geneIDs)
    pval = np.empty((num, 1))
    pval.fill(np.nan)

    explanatory0 = cm.create_matrix(data, model='H0')
    explanatory1 = cm.create_matrix(data, model='H1')
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    lenSampleRibo = data.idxRibo.size
    lenSampleRna  = data.idxRna.size

    errorCnt = 0

    for i in range(num):
        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if opts.dispDiff and np.isnan(data.dispAdjRibo[i]):
            continue
        if not opts.dispDiff and np.isnan(data.dispAdj[i]):
            continue

        response = np.hstack([data.countRibo[i, :], data.countRna[i, :]])

        if opts.dispDiff:
            disp = np.hstack([np.repeat(data.dispAdjRibo[i], lenSampleRibo), np.repeat(data.dispAdjRna[i], lenSampleRna)])
        else:
            disp = data.dispAdj[i]

        try:
            modNB0 = sm.GLM(response, explanatory0, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
            modNB1 = sm.GLM(response, explanatory1, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
            result0 = modNB0.fit()
            result1 = modNB1.fit()
        except sm.tools.sm_exceptions.PerfectSeparationError:
            errorCnt += 1
        else:
            if not opts.dispDiff:
                pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, explanatory1.shape[1] - explanatory0.shape[1])
            elif opts.dispDiff:
                pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, (explanatory1.shape[1] - explanatory0.shape[1]) / 2.5)
            else:
                pass

    data.pval = pval

    sys.stdout.write('Warning: Failed to do test: %i genes. P value set to \'nan\'.\n' % errorCnt)

    return data

def adj_pval(data, opts):
    """
    Perform multiple test correction.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input arguments to the main TE function 
    @type opts: Instance
    """

    pval = data.pval.copy()
    idx = ~np.isnan(pval)

    if opts.multiTest == 'BH':
        method = 'fdr_bh'
    elif opts.multiTest == 'Bonferroni':
        method = 'bonferroni'
    elif opts.multiTest == 'Holm':
        method = 'holm'
    elif opts.multiTest == 'Hochberg':
        method = 'simes-hochberg'
    elif opts.multiTest == 'Hommel':
        method = 'hommel'
    elif opts.multiTest == 'BY':
        method = 'fdr_by'
    elif opts.multiTest == 'TSBH':
        method = 'tsbh'
    else:
        sys.stderr.write('ERROR: The methods for multiple test correction can only accept \'Bonferroni\', \'Holm\', \'Hochberg\', \'Hommel\', \'BH\', \'BY\' or \'TSBH\' as its input.\n')
        sys.exit()

    mtc = sms.stats.multicomp.multipletests(pval[idx], alpha=0.1, method=method, returnsorted=False)

    padj = pval.copy()
    padj[idx] = mtc[1]
    data.padj = padj

    return data

def cal_TEchange(data):
    """
    Calculate Translational Efficiency changes.

    @args data: Store all input data and results
    @type data: Class object
    """

    const = 1.0

    count = np.hstack([data.countRibo, data.countRna])
    libSizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    idxRiboCtl = np.intersect1d(data.idxRibo, data.idxCtl)
    idxRnaCtl  = np.intersect1d(data.idxRna,  data.idxCtl)
    meanRiboCtl = np.mean((count/libSizes)[:, idxRiboCtl], axis=1)
    meanRnaCtl  = np.mean((count/libSizes)[:, idxRnaCtl], axis=1)
    idxC = np.logical_or(meanRiboCtl == 0.0, meanRnaCtl == 0.0)
    meanRiboCtl[idxC] = meanRiboCtl[idxC] + const
    meanRnaCtl[idxC]  = meanRnaCtl[idxC]  + const
    TEctl = meanRiboCtl / meanRnaCtl
    TEctl = np.reshape(TEctl, (TEctl.size, 1))
    data.TEctl = TEctl

    idxRiboTrt = np.intersect1d(data.idxRibo, data.idxTrt)
    idxRnaTrt  = np.intersect1d(data.idxRna,  data.idxTrt)
    meanRiboTrt = np.mean((count/libSizes)[:, idxRiboTrt], axis=1)
    meanRnaTrt  = np.mean((count/libSizes)[:, idxRnaTrt], axis=1)
    idxT = np.logical_or(meanRiboTrt == 0.0, meanRnaTrt == 0.0)
    meanRiboTrt[idxT] = meanRiboTrt[idxT] + const
    meanRnaTrt[idxT]  = meanRnaTrt[idxT]  + const
    TEtrt = meanRiboTrt / meanRnaTrt
    TEtrt = np.reshape(TEtrt, (TEtrt.size, 1))
    data.TEtrt = TEtrt

    logFoldChangeTE = np.log2(TEtrt) - np.log2(TEctl)
    logFoldChangeTE = np.reshape(logFoldChangeTE, (logFoldChangeTE.size, 1))

    data.logFoldChangeTE = logFoldChangeTE

    return data

