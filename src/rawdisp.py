import sys
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import adjlik as al
import pdb

def disper_raw(data, opts):

    cntCutoff = opts.sumCntCutoff

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    disperRawRibo = np.empty((num, 1))
    disperRawRibo.fill(np.nan)
    disperRawRNA  = disperRawRibo.copy()
    disperRawConv = disperRawRibo.copy()
    disperRawMthd = np.empty((num, 1), dtype='S10')
    disperRawMthd.fill('nan')

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    lenSampleRibo = len(data.idxRibo)
    lenSampleRNA  = len(data.idxRNA)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if sum(data.countRibo[i, :] / data.libSizesRibo) >= cntCutoff and sum(data.countRNA[i, :] / data.libSizesRNA) >= cntCutoff:

            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])

            disperInitialRibo = opts.dispInitial
            disperInitialRNA  = opts.dispInitial
            disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])

            if opts.dispDiff == 2:

                optimize_scalar = False

                mthd = 'SLSQP'
                j = 0
                while j < 10:

                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                    result = modNB.fit()

                    disperBef = disper
                    x0 = [disperBef[0], disperBef[-1]]
                    yhat = result.mu
                    sign = -1.0

                    if mthd == 'SLSQP':
                        res = minimize(al.adj_loglikelihood, x0, args=(lenSampleRibo, lenSampleRNA, explanatory, response, yhat, sign), method='SLSQP', bounds=((0, None), (0, None)), tol=1e-5)

                    if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                        if mthd == 'SLSQP':
                            j = 0
                            x0 = [disperInitialRibo, disperInitialRNA]
                            disperBef = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                        res = minimize(al.adj_loglikelihood, x0, args=(lenSampleRibo, lenSampleRNA, explanatory, response, yhat, sign), method='Nelder-Mead', tol=1e-5)
                        mthd = 'Nelder'

                    if res.success:
                        disper = np.hstack([np.repeat(res.x[0], lenSampleRibo), np.repeat(res.x[1], lenSampleRNA)])
                        if abs(np.log(disper[0]) - np.log(disperBef[0])) < 0.01 and abs(np.log(disper[-1]) - np.log(disperBef[-1])) < 0.01:

                            if mthd == 'SLSQP' and (disper[0] > 5.0 or disper[-1] > 5.0):
                                j = 0
                                mthd = 'Nelder'
                                disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                                continue

                            if np.log10(disper[0]) < -5.0 or np.log10(disper[-1]) < -5.0 or abs(abs(np.log10(disper[0])) - abs(np.log10(disper[-1]))) > 5.0:
                                optimize_scalar = True
                                break

                            disperRawRibo[i] = disper[0]
                            disperRawRNA[i]  = disper[-1]
                            disperRawConv[i] = True
                            disperRawMthd[i] = mthd
                            break
                        elif j == 9:
                            if mthd == 'SLSQP' and (disper[0] > 5.0 or disper[-1] > 5.0):
                                j = 0
                                mthd = 'Nelder'
                                disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                                continue
                            disperRawRibo[i] = disper[0]
                            disperRawRNA[i] = disper[-1]
                            disperRawConv[i] = False
                            disperRawMthd[i] = mthd
                        else:
                            pass
                    else:
                        optimize_scalar = True
                        break

                    j += 1

            if opts.dispDiff == 1 or optimize_scalar:

                disper = opts.dispInitial
                for k in range(10):
                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                    result = modNB.fit()

                    disperBef = disper
                    yhat = result.mu
                    sign = -1.0
                    res  = minimize_scalar(al.adj_loglikelihood_scalar, args=(explanatory, response, yhat, sign), method='Bounded', bounds=(0, 20), tol=1e-5)
                    disper = res.x
                    if abs(np.log(disper) - np.log(disperBef)) < 0.01:
                        disperRawRibo[i] = disper
                        disperRawRNA[i]  = disper
                        disperRawConv[i] = True
                        disperRawMthd[i] = 'Bounded'
                        break
                    elif k == 9:
                        disperRawRibo[i] = disper
                        disperRawRNA[i]  = disper
                        disperRawConv[i] = False
                        disperRawMthd[i] = 'Bounded'
                    else:
                        pass

    data.disperRawRibo = disperRawRibo
    data.disperRawRNA  = disperRawRNA
    data.disperRawConv = disperRawConv
    data.disperRawMthd = disperRawMthd

    return data

def disper_raw_scalar(data, opts):

    cntCutoff = opts.sumCntCutoff

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    disperRaw = np.empty((num, 1))
    disperRaw.fill(np.nan)
    disperRawConv = disperRaw.copy()

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if sum(data.countRibo[i, :] / data.libSizesRibo) >= cntCutoff and sum(data.countRNA[i, :] / data.libSizesRNA) >= cntCutoff:

            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])

            disper = opts.dispInitial

            for k in range(10):
                modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                result = modNB.fit()

                disperBef = disper
                yhat = result.mu
                sign = -1.0
                res = minimize_scalar(al.adj_loglikelihood_scalar, args=(explanatory, response, yhat, sign), method='Bounded', bounds=(0, 20), tol=1e-5)
                disper = res.x

                if abs(np.log(disper) - np.log(disperBef)) < 0.01:
                    disperRaw[i] = disper
                    disperRawConv[i] = True
                    break
                elif k == 9:
                    disperRaw[i] = disper
                    disperRawConv[i] = False
                else:
                    pass

    data.disperRaw = disperRaw
    data.disperRawConv = disperRawConv

    return data
