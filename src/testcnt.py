import sys
import numpy as np
import creatematrix as cm
import statsmodels.api as sm
from scipy.stats import chi2
import statsmodels.sandbox as sms
import cPickle as pickle

def test_count(data, opts):

    print 'Start the statistical test.'

    num = len(data.geneIDs)
    pval = np.empty((num, 1))
    pval.fill(np.nan)

    explanatory0 = cm.create_matrix(data, model='H0')
    explanatory1 = cm.create_matrix(data, model='H1')
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    lenSampleRibo = len(data.idxRibo)
    lenSampleRNA  = len(data.idxRNA)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num
            
        if opts.dispDiff and np.isnan(data.disperAdjRibo[i]):
            continue
        if not opts.dispDiff and np.isnan(data.disperAdj[i]):
            continue

        response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])

        if opts.dispDiff:
            disper = np.hstack([np.repeat(data.disperAdjRibo[i], lenSampleRibo), np.repeat(data.disperAdjRNA[i], lenSampleRNA)])
        else:
            disper = data.disperAdj[i]

        modNB0 = sm.GLM(response, explanatory0, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
        modNB1 = sm.GLM(response, explanatory1, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
        result0 = modNB0.fit()
        result1 = modNB1.fit()
        pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, explanatory1.shape[1] - explanatory0.shape[1])

    data.pval = pval

    return data

def adj_pval(data, opts):

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

    mtc = sms.stats.multicomp.multipletests(pval[idx], alpha=0.1, method=method, returnsorted=False)

    padj = pval.copy()
    padj[idx] = mtc[1]
    data.padj = padj

    return data

def cal_TEchange(data):

    const = 1e-5

    count = np.hstack([data.countRibo, data.countRNA])
    libSizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    idxRiboCtl = np.intersect1d(data.idxRibo, data.idxCtl)
    idxRNACtl  = np.intersect1d(data.idxRNA,  data.idxCtl)
    TEctl = (const + np.mean((count/libSizes)[:, idxRiboCtl], axis=1)) / (const + np.mean((count/libSizes)[:, idxRNACtl], axis=1))
    TEctl = np.reshape(TEctl, (len(TEctl), 1))
    data.TEctl = TEctl

    idxRiboTrt = np.intersect1d(data.idxRibo, data.idxTrt)
    idxRNATrt  = np.intersect1d(data.idxRNA,  data.idxTrt)
    TEtrt = (const + np.mean((count/libSizes)[:, idxRiboTrt], axis=1)) / (const + np.mean((count/libSizes)[:, idxRNATrt], axis=1))
    TEtrt = np.reshape(TEtrt, (len(TEtrt), 1))
    data.TEtrt = TEtrt

    logFoldChangeTE = np.log2(TEtrt) - np.log2(TEctl)
    logFoldChangeTE = np.reshape(logFoldChangeTE, (len(logFoldChangeTE), 1))

    data.logFoldChangeTE = logFoldChangeTE

    return data
