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
    beta = np.array([0.1, 0.01])

    iter = 5

    if obj == 'RR':
        index = np.nonzero(dispRawConv == True)[0]
    else:
        index = np.logical_and(dispRawConv == True, dispRawMthd != 'Bounded').nonzero()[0]

    lowerBound = np.percentile(np.unique(dispRaw[index]),  1)
    upperBound = np.percentile(np.unique(dispRaw[index]), 99)

    idx = np.logical_and(dispRaw > lowerBound, dispRaw < upperBound).nonzero()[0]

    for i in range(iter):

        matrix = np.empty((idx.size, 2))
        matrix.fill(np.nan)
        matrix[:, 0] = 1 / countMean[idx].flatten()
        matrix[:, 1] = 1

        modGamma = sm.GLM(dispRaw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
        result = modGamma.fit(start_params=beta)

        betaBef = beta
        beta = result.params

        if sum(np.log(beta / betaBef)**2) < 1e-4:
            dispFittedConv = True
            break
        if i == iter - 1:
            dispFittedConv = False
            if obj != 'RR':
                print 'Fitting dispersion for %s does not converge.' % obj
            else:
                print 'Fitting dispersion does not converge.'

    dispFitted = dispRaw.copy()
    IDX = ~np.isnan(dispRaw)
    dispFitted[IDX] = beta[0] / countMean[IDX] + beta[1]
    if np.nonzero(dispFitted < 0)[0].size > 0:
        print 'Negative fitted dispersion exist!' 

    if obj == 'Ribo':
        data.dispFittedRibo = dispFitted
        data.betaRibo = beta
        data.dispFittedRiboConv = dispFittedConv
        data.dispFittedRiboIdx = idx
    elif obj == 'mRNA':
        data.dispFittedRna = dispFitted
        data.betaRna = beta
        data.dispFittedRnaConv = dispFittedConv
        data.dispFittedRnaIdx = idx
    elif obj == 'RR':
        data.dispFitted = dispFitted
        data.beta = beta
        data.dispFittedConv = dispFittedConv
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
