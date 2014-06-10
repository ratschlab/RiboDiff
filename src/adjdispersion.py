import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from scipy.special import polygamma
import adjlikelihood as al

def calculate_logprior(disper, disperFitted, varPrior):

    logprior = (np.log(disper) - np.log(disperFitted)) ** 2 / (2 * varPrior ** 2)

    return logprior

def adj_loglikelihood_shrink(disper, explanatory, response, yhat, disperFitted, varPrior, sign):

    loglik_adj = al.adj_loglikelihood(disper, explanatory, response, yhat, 1.0)
    logprior = calculate_logprior(disper, disperFitted, varPrior)
    loglik_adj_shrk = loglik_adj - logprior

    return loglik_adj_shrk * sign

def disper_adj(data):

    print 'Start to estimate adjusted dispersions.'

    num = len(data.geneIDs)
    disperAdj = np.empty((num, 1))
    disperAdj.fill(np.nan)

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    matrix = data.matrix
    numSample = matrix.shape[0]
    numCoef = matrix.shape[1]
    varLogDispSamp = polygamma(1, (numSample - numCoef)/2)

    index = ~np.isnan(data.disperRaw)
    LogResidule = np.log(data.disperRaw[index]) - np.log(data.disperFitted[index])
    varLogResidule = (np.median(np.abs(LogResidule) - np.median(LogResidule)) * 1.4826 ) ** 2

    varPrior = min((varLogResidule - varLogDispSamp), 0.25)

    for i in range(num):

        if i % 50 == 0:
            print '%i genes finished...' % i
        if i+1 == num:
            print '%i genes finished...' % num

        if not np.isnan(data.disperRaw[i]):
            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
            disper = 0.1
            for j in range(10):
                modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                result = modNB.fit()
                disperBef = disper
                yhat = result.mu
                disperFitted = data.disperFitted[i]
                sign = -1.0
                res = minimize_scalar(adj_loglikelihood_shrink, bounds=(0, 1), args=(explanatory, response, yhat, disperFitted, varPrior, sign), tol=1e-5, method='Bounded')
                disper = res.x
                if abs(np.log(disper) - np.log(disperBef)) < 0.03:
                    break
            disperAdj[i] = disper

    data.disperAdj = disperAdj

    return data
