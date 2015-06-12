#!/usr/bin/env python
"""
Estimating the raw dispersion.
"""

import sys
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import adjlik as al

def disper_raw(data, opts):

    cntCutoff = opts.sumCntCutoff

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    muRaw = np.empty((num, data.idxRibo.size + data.idxRna.size))
    muRawRibo = np.empty((num, data.idxRibo.size))
    muRawRna  = np.empty((num, data.idxRna.size))
    dispRawRibo = np.empty((num, 1))
    dispRawRibo.fill(np.nan)
    dispRawRna  = dispRawRibo.copy()
    dispRawConv = dispRawRibo.copy()
    dispRawMthd = np.empty((num, 1), dtype='S10')
    dispRawMthd.fill('nan')

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    lenSampleRibo = data.idxRibo.size
    lenSampleRna  = data.idxRna.size

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if sum(data.countRibo[i, :] / data.libSizesRibo) >= cntCutoff and sum(data.countRna[i, :] / data.libSizesRna) >= cntCutoff:

            response = np.hstack([data.countRibo[i, :], data.countRna[i, :]])

            dispInitialRibo = opts.dispInitial
            dispInitialRna  = opts.dispInitial
            disp = np.hstack([np.repeat(dispInitialRibo, lenSampleRibo), np.repeat(dispInitialRna, lenSampleRna)])

            optimize_scalar = False

            mthd = 'SLSQP'
            j = 0
            while j < 10:
                try:
                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    x0 = [dispBef[0], dispBef[-1]]
                    yhat = result.mu
                    sign = -1.0

                    if mthd == 'SLSQP':
                        res = minimize(al.adj_loglikelihood, x0, args=(lenSampleRibo, lenSampleRna, explanatory, response, yhat, sign), method='SLSQP', bounds=((0, 10.0), (0, 10.0)), tol=1e-5)

                    if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                        if mthd == 'SLSQP':
                            j = 0
                            x0 = [dispInitialRibo, dispInitialRna]
                            dispBef = np.hstack([np.repeat(dispInitialRibo, lenSampleRibo), np.repeat(dispInitialRna, lenSampleRna)])
                        res = minimize(al.adj_loglikelihood, x0, args=(lenSampleRibo, lenSampleRna, explanatory, response, yhat, sign), method='Nelder-Mead', tol=1e-5)
                        mthd = 'Nelder'

                    if res.success:
                        disp = np.hstack([np.repeat(res.x[0], lenSampleRibo), np.repeat(res.x[1], lenSampleRna)])
                        if abs(np.log(disp[0]) - np.log(dispBef[0])) < 0.01 and abs(np.log(disp[-1]) - np.log(dispBef[-1])) < 0.01:
                            dispRawRibo[i] = disp[0]
                            dispRawRna[i]  = disp[-1]
                            dispRawConv[i] = True
                            dispRawMthd[i] = mthd
                            muRaw[i, :] = yhat
                            muRawRibo[i, :] = yhat[data.idxRibo]
                            muRawRna[i, :] = yhat[data.idxRna]
                            break
                        elif j == 9:
                            dispRawRibo[i] = disp[0]
                            dispRawRna[i] = disp[-1]
                            dispRawConv[i] = False
                            dispRawMthd[i] = mthd
                            muRaw[i, :] = yhat
                            muRawRibo[i, :] = yhat[data.idxRibo]
                            muRawRna[i, :] = yhat[data.idxRna]
                        else:
                            pass
                    else:
                        optimize_scalar = True
                        break

                except sm.tools.sm_exceptions.PerfectSeparationError:
                    dispRawRibo[i] = disp[0]
                    dispRawRna[i] = disp[-1]
                    dispRawConv[i] = False
                    dispRawMthd[i] = mthd
                    muRaw[i, :] = yhat
                    muRawRibo[i, :] = yhat[data.idxRibo]
                    muRawRna[i, :] = yhat[data.idxRna]

                j += 1

            if optimize_scalar:

                disp = opts.dispInitial
                for k in range(10):
                    try:
                        modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                        result = modNB.fit()

                        dispBef = disp
                        yhat = result.mu
                        sign = -1.0
                        res  = minimize_scalar(al.adj_loglikelihood_scalar, args=(explanatory, response, yhat, sign), method='Bounded', bounds=(0.0, 10.0), tol=1e-5)
                        disp = res.x
                        if abs(np.log(disp) - np.log(dispBef)) < 1e-4:
                            dispRawRibo[i] = disp
                            dispRawRna[i]  = disp
                            dispRawConv[i] = True
                            dispRawMthd[i] = 'Bounded'
                            muRaw[i, :] = yhat
                            muRawRibo[i, :] = yhat[data.idxRibo]
                            muRawRna[i, :] = yhat[data.idxRna]
                            break
                        elif k == 9:
                            dispRawRibo[i] = disp
                            dispRawRna[i]  = disp
                            dispRawConv[i] = False
                            dispRawMthd[i] = 'Bounded'
                            muRaw[i, :] = yhat
                            muRawRibo[i, :] = yhat[data.idxRibo]
                            muRawRna[i, :] = yhat[data.idxRna]
                        else:
                            pass
                    except sm.tools.sm_exceptions.PerfectSeparationError:
                        dispRawRibo[i] = disp
                        dispRawRna[i]  = disp
                        dispRawConv[i] = False
                        dispRawMthd[i] = 'Bounded'
                        muRaw[i, :] = yhat
                        muRawRibo[i, :] = yhat[data.idxRibo]
                        muRawRna[i, :] = yhat[data.idxRna]

    data.muRaw     = muRaw
    data.muRawRibo = muRawRibo
    data.muRawRna  = muRawRna
    data.dispRawRibo = dispRawRibo
    data.dispRawRna  = dispRawRna
    data.dispRawConv = dispRawConv
    data.dispRawMthd = dispRawMthd

    return data

def disper_raw_scalar(data, opts):

    cntCutoff = opts.sumCntCutoff

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    muRaw = np.empty((num, data.idxRibo.size + data.idxRna.size))
    muRawRibo = np.empty((num, data.idxRibo.size))
    muRawRna = np.empty((num, data.idxRna.size))
    dispRaw = np.empty((num, 1))
    dispRaw.fill(np.nan)
    dispRawConv = dispRaw.copy()

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if sum(data.countRibo[i, :] / data.libSizesRibo) >= cntCutoff and sum(data.countRna[i, :] / data.libSizesRna) >= cntCutoff:

            response = np.hstack([data.countRibo[i, :], data.countRna[i, :]])

            disp = opts.dispInitial

            for k in range(10):
                try:
                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    yhat = result.mu
                    sign = -1.0
                    res = minimize_scalar(al.adj_loglikelihood_scalar, args=(explanatory, response, yhat, sign), method='Bounded', bounds=(0, 10.0), tol=1e-5)
                    disp = res.x

                    if abs(np.log(disp) - np.log(dispBef)) < 1e-4:
                        dispRaw[i] = disp
                        dispRawConv[i] = True
                        muRaw[i, :] = yhat
                        muRawRibo[i, :] = yhat[data.idxRibo]
                        muRawRna[i, :] = yhat[data.idxRna]
                        break
                    elif k == 9:
                        dispRaw[i] = disp
                        dispRawConv[i] = False
                        muRaw[i, :] = yhat
                        muRawRibo[i, :] = yhat[data.idxRibo]
                        muRawRna[i, :] = yhat[data.idxRna]
                    else:
                        pass
                except sm.tools.sm_exceptions.PerfectSeparationError:
                        dispRaw[i] = disp
                        dispRawConv[i] = False
                        muRaw[i, :] = yhat
                        muRawRibo[i, :] = yhat[data.idxRibo]
                        muRawRna[i, :] = yhat[data.idxRna]

    data.muRaw     = muRaw
    data.muRawRibo = muRawRibo
    data.muRawRna  = muRawRna
    data.dispRaw = dispRaw
    data.dispRawConv = dispRawConv

    return data

