#!/usr/bin/env python
"""
Estimating raw dispersion.
"""

import sys
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.special import polygamma
import adjlik as al

def calculate_logprior(disp, dispFitted, varPrior):
    """
    """

    logprior = (np.log(disp) - np.log(dispFitted)) ** 2 / (2 * varPrior ** 2)

    return logprior

def adj_loglikelihood_shrink(x, lenSampleRibo, lenSampleRna, explanatory, response, yhat, dispFittedRibo, dispFittedRna, varPriorRibo, varPriorRna, sign):
    """
    """

    loglik_adj = al.adj_loglikelihood(x, lenSampleRibo, lenSampleRna, explanatory, response, yhat, 1.0)
    logpriorRibo = calculate_logprior(x[0], dispFittedRibo, varPriorRibo)
    logpriorRna  = calculate_logprior(x[1], dispFittedRna, varPriorRna)
    loglik_adj_shrk = loglik_adj - logpriorRibo - logpriorRna

    return loglik_adj_shrk * sign

def adj_loglikelihood_shrink_scalar(disp, explanatory, response, yhat, dispFittedRibo, dispFittedRna, varPriorRibo, varPriorRna, sign):
    """
    """

    loglik_adj = al.adj_loglikelihood_scalar(disp, explanatory, response, yhat, 1.0)
    logpriorRibo = calculate_logprior(disp, dispFittedRibo, varPriorRibo)
    logpriorRna  = calculate_logprior(disp, dispFittedRna, varPriorRna)
    loglik_adj_shrk = loglik_adj - logpriorRibo - logpriorRna

    return loglik_adj_shrk * sign

def adj_loglikelihood_shrink_scalar_onedisper(disp, explanatory, response, yhat, dispFitted, varPrior, sign):
    """
    """

    loglik_adj = al.adj_loglikelihood_scalar(disp, explanatory, response, yhat, 1.0)
    logprior = calculate_logprior(disp, dispFitted, varPrior)
    loglik_adj_shrk = loglik_adj - logprior

    return loglik_adj_shrk * sign

def calculate_varPrior(dispRaw, dispFitted, dispFittedIdx, varLogDispSamp):
    """
    """

    logResidule = np.log(dispRaw[dispFittedIdx]) - np.log(dispFitted[dispFittedIdx])
    stdLogResidule = np.median(np.abs(logResidule - np.median(logResidule))) * 1.4826

    varLogResidule = stdLogResidule ** 2
    varPrior = varLogResidule - varLogDispSamp

    varPrior = max(varPrior, 0.1)

    return varPrior

def disper_adj(data, opts):
    """
    """

    print 'Start to estimate adjusted dispersions.'

    num = len(data.geneIDs)
    muAdj = np.empty((num, data.idxRibo.size + data.idxRna.size))
    muAdjRibo = np.empty((num, data.idxRibo.size))
    muAdjRna  = np.empty((num, data.idxRna.size))
    dispAdjRibo = np.empty((num, 1))
    dispAdjRibo.fill(np.nan)
    dispAdjRna  = dispAdjRibo.copy()
    dispAdjConv = dispAdjRibo.copy()
    dispAdjMthd = np.empty((num, 1), dtype='S10')
    dispAdjMthd.fill('nan')

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    lenSampleRibo = data.idxRibo.size
    lenSampleRna  = data.idxRna.size

    matrix = data.matrix
    numSample = matrix.shape[0]
    numCoef = matrix.shape[1]
    varLogDispSamp = polygamma(1, (numSample - numCoef)/2)

    varPriorRibo = calculate_varPrior(data.dispRawRibo, data.dispFittedRibo, data.dispFittedRiboIdx, varLogDispSamp)
    varPriorRna  = calculate_varPrior(data.dispRawRna,  data.dispFittedRna,  data.dispFittedRnaIdx,  varLogDispSamp)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if not np.isnan(data.dispRawRibo[i]):

            response = np.hstack([data.countRibo[i, :], data.countRna[i, :]])

            dispInitialRibo = opts.dispInitial
            dispInitialRna  = opts.dispInitial
            disp = np.hstack([np.repeat(dispInitialRibo, lenSampleRibo), np.repeat(dispInitialRna, lenSampleRna)])

            dispFittedRibo = data.dispFittedRibo[i]
            dispFittedRna  = data.dispFittedRna[i]
            dispRawRibo = data.dispRawRibo[i]
            dispRawRna  = data.dispRawRna[i]

            optimize_scalar = False
            mthd = 'SLSQP'
            j = 0
            while j < 10:
                try:
                    modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    x0 = [dispBef[0], dispBef[-1]]
                    yhat = result.mu
                    sign = -1.0

                    if mthd == 'SLSQP':
                        res = minimize(adj_loglikelihood_shrink, x0, args=(lenSampleRibo, lenSampleRna, explanatory, response, yhat, dispFittedRibo, dispFittedRna, varPriorRibo, varPriorRna, sign), method='SLSQP', bounds=((0, 10.0), (0, 10.0)), tol=1e-5)

                    if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                        if mthd == 'SLSQP':
                            j = 0
                            x0 = [dispInitialRibo, dispInitialRna]
                            dispBef = np.hstack([np.repeat(dispInitialRibo, lenSampleRibo), np.repeat(dispInitialRna, lenSampleRna)])
                        res = minimize(adj_loglikelihood_shrink, x0, args=(lenSampleRibo, lenSampleRna, explanatory, response, yhat, dispFittedRibo, dispFittedRna, varPriorRibo, varPriorRna, sign), method='Nelder-Mead', tol=1e-5)
                        mthd = 'Nelder'

                    if res.success:
                        disp = np.hstack([np.repeat(res.x[0], lenSampleRibo), np.repeat(res.x[1], lenSampleRna)])
                        if abs(np.log(disp[0]) - np.log(dispBef[0])) < 0.01 and abs(np.log(disp[-1]) - np.log(dispBef[-1])) < 0.01:
                            dispAdjRibo[i] = disp[0]
                            dispAdjRna[i]  = disp[-1]
                            dispAdjConv[i] = True
                            dispAdjMthd[i] = mthd
                            muAdj[i, :] = yhat
                            muAdjRibo[i, :] = yhat[data.idxRibo]
                            muAdjRna[i, :] = yhat[data.idxRna]
                            break
                        elif j == 9:
                            dispAdjRibo[i] = disp[0]
                            dispAdjRna[i]  = disp[-1]
                            dispAdjConv[i] = False
                            dispAdjMthd[i] = mthd
                            muAdj[i, :] = yhat
                            muAdjRibo[i, :] = yhat[data.idxRibo]
                            muAdjRna[i, :] = yhat[data.idxRna]
                        else:
                            pass
                    else:
                        optimize_scalar = True
                        break

                except sm.tools.sm_exceptions.PerfectSeparationError:
                    dispAdjRibo[i] = disp[0]
                    dispAdjRna[i] = disp[-1]
                    dispAdjConv[i] = False
                    dispAdjMthd[i] = mthd
                    muAdj[i, :] = yhat
                    muAdjRibo[i, :] = yhat[data.idxRibo]
                    muAdjRna[i, :] = yhat[data.idxRna]

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
                        res  = minimize_scalar(adj_loglikelihood_shrink_scalar, args=(explanatory, response, yhat, dispFittedRibo, dispFittedRna, varPriorRibo, varPriorRna, sign), method='Bounded', bounds=(0.0, 10.0), tol=1e-5)
                        disp = res.x

                        if abs(np.log(disp) - np.log(dispBef)) < 1e-4:
                            dispAdjRibo[i] = disp
                            dispAdjRna[i]  = disp
                            dispAdjConv[i] = True
                            dispAdjMthd[i] = 'Bounded'
                            muAdj[i, :] = yhat
                            muAdjRibo[i, :] = yhat[data.idxRibo]
                            muAdjRna[i, :] = yhat[data.idxRna]
                            break
                        elif k == 9:
                            dispAdjRibo[i] = disp
                            dispAdjRna[i]  = disp
                            dispAdjConv[i] = False
                            dispAdjMthd[i] = 'Bounded'
                            muAdj[i, :] = yhat
                            muAdjRibo[i, :] = yhat[data.idxRibo]
                            muAdjRna[i, :] = yhat[data.idxRna]
                        else:
                            pass
                    except sm.tools.sm_exceptions.PerfectSeparationError:
                        dispAdjRibo[i] = disp
                        dispAdjRna[i]  = disp
                        dispAdjConv[i] = False
                        dispAdjMthd[i] = 'Bounded'
                        muAdj[i, :] = yhat
                        muAdjRibo[i, :] = yhat[data.idxRibo]
                        muAdjRna[i, :] = yhat[data.idxRna]

    data.muAdj     = muAdj
    data.muAdjRibo = muAdjRibo
    data.muAdjRna  = muAdjRna
    data.dispAdjRibo = dispAdjRibo
    data.dispAdjRna  = dispAdjRna
    data.dispAdjConv = dispAdjConv
    data.dispAdjMthd = dispAdjMthd

    return data

def disper_adj_scalar(data, opts):

    print 'Start to estimate adjusted dispersions.'

    num = len(data.geneIDs)
    muAdj = np.empty((num, data.idxRibo.size + data.idxRna.size))
    muAdjRibo = np.empty((num, data.idxRibo.size))
    muAdjRna = np.empty((num, data.idxRna.size))
    dispAdj = np.empty((num, 1))
    dispAdj.fill(np.nan)
    dispAdjConv = dispAdj.copy()

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRna])

    matrix = data.matrix
    numSample = matrix.shape[0]
    numCoef = matrix.shape[1]
    varLogDispSamp = polygamma(1, (numSample - numCoef)/2)

    varPrior = calculate_varPrior(data.dispRaw, data.dispFitted, data.dispFittedIdx, varLogDispSamp)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if not np.isnan(data.dispRaw[i]):

            response = np.hstack([data.countRibo[i, :], data.countRna[i, :]])

            disp = opts.dispInitial

            for j in range(10):
                try:
                    modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disp), offset=np.log(librarySizes))
                    result = modNB.fit()

                    dispBef = disp
                    yhat = result.mu
                    dispFitted = data.dispFitted[i]
                    dispRaw = data.dispRaw[i]
                    sign = -1.0
                    res = minimize_scalar(adj_loglikelihood_shrink_scalar_onedisper, args=(explanatory, response, yhat, dispFitted, varPrior, sign), method='Bounded', bounds=(0, 10.0), tol=1e-5)
                    disp = res.x

                    if abs(np.log(disp) - np.log(dispBef)) < 1e-4:
                        dispAdj[i] = disp
                        dispAdjConv[i] = True
                        muAdj[i, :] = yhat
                        muAdjRibo[i, :] = yhat[data.idxRibo]
                        muAdjRna[i, :] = yhat[data.idxRna]
                        break
                    elif j == 9:
                        dispAdj[i] = disp
                        dispAdjConv[i] = False
                        muAdj[i, :] = yhat
                        muAdjRibo[i, :] = yhat[data.idxRibo]
                        muAdjRna[i, :] = yhat[data.idxRna]
                    else:
                        pass
                except sm.tools.sm_exceptions.PerfectSeparationError:
                        dispAdj[i] = disp
                        dispAdjConv[i] = False
                        muAdj[i, :] = yhat
                        muAdjRibo[i, :] = yhat[data.idxRibo]
                        muAdjRna[i, :] = yhat[data.idxRna]

    data.muAdj = muAdj
    data.muAdjRibo = muAdjRibo
    data.muAdjRna  = muAdjRna
    data.dispAdj = dispAdj
    data.dispAdjConv = dispAdjConv

    return data

