#!/usr/bin/env python
"""
Fitting the raw dispersion to a gamma regression.
"""

import numpy as np
import statsmodels.api as sm

def do_fitting(data, obj):

    if obj == 'Ribo':
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        dispRaw = data.dispRawRibo
    elif obj == 'mRNA':
        countMean = np.mean(data.countRna  / data.libSizesRna,  axis=1)
        dispRaw = data.dispRawRna
    elif obj == 'RR':
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        dispRaw = data.dispRaw
    else:
        pass

    countMean = np.reshape(countMean, (countMean.size, 1))
    dispRawConv = data.dispRawConv
    dispRawMthd = data.dispRawMthd

    if obj == 'RR':
        index = np.nonzero(dispRawConv == True)[0]
    else:
        index = np.logical_and(dispRawConv == True, dispRawMthd != 'Bounded').nonzero()[0]

    lowerBound = np.percentile(np.unique(dispRaw[index]),  1)
    upperBound = np.percentile(np.unique(dispRaw[index]), 99)

    idx = np.logical_and(dispRaw > lowerBound, dispRaw < upperBound).nonzero()[0]

    matrix = np.empty((idx.size, 2))
    matrix.fill(np.nan)
    matrix[:, 0] = 1 / countMean[idx].flatten()
    matrix[:, 1] = 1

    modGamma = sm.GLM(dispRaw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
    result = modGamma.fit()
    Lambda = result.params

    dispFitted = dispRaw.copy()
    IDX = ~np.isnan(dispRaw)
    dispFitted[IDX] = Lambda[0] / countMean[IDX] + Lambda[1]
    if np.nonzero(dispFitted < 0)[0].size > 0:
        print 'Negative fitted dispersion exist!' 

    if obj == 'Ribo':
        data.dispFittedRibo = dispFitted
        data.LambdaRibo = Lambda
        data.dispFittedRiboIdx = idx
    elif obj == 'mRNA':
        data.dispFittedRna = dispFitted
        data.LambdaRna = Lambda
        data.dispFittedRnaIdx = idx
    elif obj == 'RR':
        data.dispFitted = dispFitted
        data.Lambda = Lambda
        data.dispFittedIdx = idx
    else:
        pass

    if obj != 'RR':
        print 'Fit dispersion for %s: Done.' % obj
    else:
        print 'Fit dispersion: Done.'

    return data

def disper_fit(data, opts):

    if opts.dispDiff:
        data = do_fitting(data, obj='Ribo')
        print '*'*25
        data = do_fitting(data, obj='mRNA')
    else:
        data = do_fitting(data, obj='RR')

    return data
