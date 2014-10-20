import sys
import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.special import polygamma
import adjlik as al

def calculate_logprior(disper, disperFitted, varPrior):

    logprior = (np.log(disper) - np.log(disperFitted)) ** 2 / (2 * varPrior ** 2)

    return logprior

def adj_loglikelihood_shrink(x, lenSampleRibo, lenSampleRNA, explanatory, response, yhat, disperFittedRibo, disperFittedRNA, varPriorRibo, varPriorRNA, sign):

    loglik_adj = al.adj_loglikelihood(x, lenSampleRibo, lenSampleRNA, explanatory, response, yhat, 1.0)
    logpriorRibo = calculate_logprior(x[0], disperFittedRibo, varPriorRibo)
    logpriorRNA  = calculate_logprior(x[1], disperFittedRNA, varPriorRNA)
    loglik_adj_shrk = loglik_adj - logpriorRibo - logpriorRNA

    return loglik_adj_shrk * sign

def adj_loglikelihood_shrink_scalar(disper, explanatory, response, yhat, disperFittedRibo, disperFittedRNA, varPriorRibo, varPriorRNA, sign):

    loglik_adj = al.adj_loglikelihood_scalar(disper, explanatory, response, yhat, 1.0)
    logpriorRibo = calculate_logprior(disper, disperFittedRibo, varPriorRibo)
    logpriorRNA  = calculate_logprior(disper, disperFittedRNA, varPriorRNA)
    loglik_adj_shrk = loglik_adj - logpriorRibo - logpriorRNA

    return loglik_adj_shrk * sign

def adj_loglikelihood_shrink_scalar_onedisper(disper, explanatory, response, yhat, disperFitted, varPrior, sign):

    loglik_adj = al.adj_loglikelihood_scalar(disper, explanatory, response, yhat, 1.0)
    logprior = calculate_logprior(disper, disperFitted, varPrior)
    loglik_adj_shrk = loglik_adj - logprior

    return loglik_adj_shrk * sign

def calculate_varPrior(disperRaw, disperFitted, disperFittedIdx, varLogDispSamp):

    LogResidule = np.log(disperRaw[disperFittedIdx]) - np.log(disperFitted[disperFittedIdx])
    stdLogResidule = np.median(np.abs(LogResidule - np.median(LogResidule))) * 1.4826

    varLogResidule = stdLogResidule ** 2
    varPrior = varLogResidule - varLogDispSamp

    varPrior = max(varPrior, 0.1)

    return varPrior

def disper_adj(data, opts):

    print 'Start to estimate adjusted dispersions.'

    num = len(data.geneIDs)
    disperAdjRibo = np.empty((num, 1))
    disperAdjRibo.fill(np.nan)
    disperAdjRNA  = disperAdjRibo.copy()
    disperAdjConv = disperAdjRibo.copy()
    disperAdjMthd = np.empty((num, 1), dtype='S10')
    disperAdjMthd.fill('nan')

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    lenSampleRibo = len(data.idxRibo)
    lenSampleRNA  = len(data.idxRNA)

    matrix = data.matrix
    numSample = matrix.shape[0]
    numCoef = matrix.shape[1]
    varLogDispSamp = polygamma(1, (numSample - numCoef)/2)

    varPriorRibo = calculate_varPrior(data.disperRawRibo, data.disperFittedRibo, data.disperFittedRiboIdx, varLogDispSamp)
    varPriorRNA  = calculate_varPrior(data.disperRawRNA,  data.disperFittedRNA,  data.disperFittedRNAIdx,  varLogDispSamp)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if not np.isnan(data.disperRawRibo[i]):

            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])

            disperInitialRibo = opts.dispInitial
            disperInitialRNA  = opts.dispInitial
            disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])

            optimize_scalar = False
            mthd = 'SLSQP'
            j = 0
            while j < 10:

                modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                result = modNB.fit()

                disperBef = disper
                x0 = [disperBef[0], disperBef[-1]]
                yhat = result.mu
                disperFittedRibo = data.disperFittedRibo[i]
                disperFittedRNA  = data.disperFittedRNA[i]
                sign = -1.0

                if mthd == 'SLSQP':
                    res = minimize(adj_loglikelihood_shrink, x0, args=(lenSampleRibo, lenSampleRNA, explanatory, response, yhat, disperFittedRibo, disperFittedRNA, varPriorRibo, varPriorRNA, sign), method='SLSQP', bounds=((0, None), (0, None)), tol=1e-5)

                if mthd == 'Nelder' or (mthd == 'SLSQP' and (not res.success or np.isnan(res.fun))):
                    if mthd == 'SLSQP':
                        j = 0
                        x0 = [disperInitialRibo, disperInitialRNA]
                        disperBef = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                    res = minimize(adj_loglikelihood_shrink, x0, args=(lenSampleRibo, lenSampleRNA, explanatory, response, yhat, disperFittedRibo, disperFittedRNA, varPriorRibo, varPriorRNA, sign), method='Nelder-Mead', tol=1e-5)
                    mthd = 'Nelder'

                if res.success:
                    disper = np.hstack([np.repeat(res.x[0], lenSampleRibo), np.repeat(res.x[1], lenSampleRNA)])
                    if abs(np.log(disper[0]) - np.log(disperBef[0])) < 0.01 and abs(np.log(disper[-1]) - np.log(disperBef[-1])) < 0.01:
                        if mthd == 'SLSQP' and (disper[0] > 5.0 or disper[-1] > 5.0):
                            j = 0
                            mthd = 'Nelder'
                            disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                            continue
                        disperAdjRibo[i] = disper[0]
                        disperAdjRNA[i]  = disper[-1]
                        disperAdjConv[i] = True
                        disperAdjMthd[i] = mthd
                        break
                    elif j == 9:
                        if mthd == 'SLSQP' and (disper[0] > 5.0 or disper[-1] > 5.0):
                            j = 0
                            mthd = 'Nelder'
                            disper = np.hstack([np.repeat(disperInitialRibo, lenSampleRibo), np.repeat(disperInitialRNA, lenSampleRNA)])
                            continue
                        disperAdjRibo[i] = disper[0]
                        disperAdjRNA[i]  = disper[-1]
                        disperAdjConv[i] = False
                        disperAdjMthd[i] = mthd
                    else:
                        pass
                else:
                    optimize_scalar = True
                    break

                j += 1

            if optimize_scalar:

                disper = opts.dispInitial
                for k in range(10):
                    modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                    result = modNB.fit()

                    disperBef = disper
                    yhat = result.mu
                    sign = -1.0
                    res  = minimize_scalar(adj_loglikelihood_shrink_scalar, args=(explanatory, response, yhat, disperFittedRibo, disperFittedRNA, varPriorRibo, varPriorRNA, sign), method='Bounded', bounds=(0, 20), tol=1e-5)
                    disper = res.x
                    if abs(np.log(disper) - np.log(disperBef)) < 0.01:
                        disperAdjRibo[i] = disper
                        disperAdjRNA[i]  = disper
                        disperAdjConv[i] = True
                        disperAdjMthd[i] = 'Bounded'
                        break
                    elif k == 9:
                        disperAdjRibo[i] = disper
                        disperAdjRNA[i]  = disper
                        disperAdjConv[i] = False
                        disperAdjMthd[i] = 'Bounded'
                    else:
                        pass

    data.disperAdjRibo = disperAdjRibo
    data.disperAdjRNA  = disperAdjRNA
    data.disperAdjConv = disperAdjConv
    data.disperAdjMthd = disperAdjMthd

    return data

def disper_adj_scalar(data, opts):

    print 'Start to estimate adjusted dispersions.'

    num = len(data.geneIDs)
    disperAdj = np.empty((num, 1))
    disperAdj.fill(np.nan)
    disperAdjConv = disperAdj.copy()

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    matrix = data.matrix
    numSample = matrix.shape[0]
    numCoef = matrix.shape[1]
    varLogDispSamp = polygamma(1, (numSample - numCoef)/2)

    varPrior = calculate_varPrior(data.disperRaw, data.disperFitted, data.disperFittedIdx, varLogDispSamp)

    for i in range(num):

        sys.stdout.flush()

        if i % 50 == 0:
            print '\r%i genes finished...' % i ,
        if i+1 == num:
            print '\r%i genes finished.' % num

        if not np.isnan(data.disperRaw[i]):
            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])

            disper = opts.dispInitial

            for j in range(10):
                modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                result = modNB.fit()

                disperBef = disper
                yhat = result.mu
                disperFitted = data.disperFitted[i]
                sign = -1.0

                res  = minimize_scalar(adj_loglikelihood_shrink_scalar_onedisper, args=(explanatory, response, yhat, disperFitted, varPrior, sign), method='Bounded', bounds=(0, 20), tol=1e-5)
                disper = res.x
                if abs(np.log(disper) - np.log(disperBef)) < 0.01:
                    disperAdj[i] = disper
                    disperAdjConv[i] = True
                    break
                elif j == 9:
                    disperAdj[i] = disper
                    disperAdjConv[i] = False
                else:
                    pass

    data.disperAdj = disperAdj
    data.disperAdjConv = disperAdjConv

    return data
