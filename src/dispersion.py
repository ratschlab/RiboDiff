import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from scipy.stats import nbinom

def adj_loglikelihood_simp(disper, explanatory, response, yhat, sign):

    n = 1 / disper
    p = n / (n + yhat)
    loglik = sum(nbinom.logpmf(response, n, p))

    mu = yhat
    var = mu + disper * mu ** 2
    w = 1 / ((1 / mu) ** 2 * var)
    w = w.reshape(len(w), 1)
    q, r = np.linalg.qr(explanatory * w ** 0.5)
    coxreid = sum(np.log(abs(np.diag(r)[range(np.linalg.matrix_rank(explanatory))])))

    return (loglik - coxreid) * sign

def disper_simp(data):

    cntCutoff = 10

    num = len(data.geneIDs)
    disperRaw = np.empty((num, 1))
    disperRaw.fill(np.nan)

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    for i in range(num):
        if sum(data.countRibo[i, :]) >= cntCutoff and sum(data.countRNA[i, :]) >= cntCutoff:
            response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
            disper = 0.1
            for j in range(10):
                modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
                #modNB = sm.NegativeBinomial(response, explanatory, loglike_method='nb2')
                result = modNB.fit()
                #print result.summary()
                #print result.mu
                disperBef = disper
                yhat = result.mu
                sign = -1.0
                res  = minimize_scalar(adj_loglikelihood_simp, bounds=(0, 1), args=(explanatory, response, yhat, sign), tol=1e-5, method='Bounded')
                disper = res.x
                if abs(np.log(disper) - np.log(disperBef)) < 0.03:
                    break
            disperRaw[i] = disper

    data.disperRaw = disperRaw

    return data

def fit_disper_gamma(data):

    countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countRiboMean=np.reshape(countRiboMean, (len(countRiboMean), 1))
    disperRaw = data.disperRaw
    #disperRaw = disperRaw.flatten()
    beta = np.array([0.1, 1])

    iter = 10
    for i in range(iter):
        ratio = disperRaw / (beta[0] + beta[1]/countRiboMean)
        idx = np.logical_and(ratio>1e-3, ratio<10)

        matrix = np.empty((len(np.nonzero(idx)[0]), 2))
        matrix.fill(np.nan)
        matrix[:, 0] = 1
        matrix[:, 1] = 1 / countRiboMean[idx]

        modGamma = sm.GLM(disperRaw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
        result = modGamma.fit(start_params=beta)

        #print result.summary()
        #print result.params

        betaBef = beta
        beta = result.params

        if sum(np.log(beta / betaBef)**2) < 0.001:
            break

        if i == iter - 1:
            print 'Fitting dispersion does not converge.'

    # Add codes for dispersion shrinkage.
    disperFitted = disperRaw.copy()
    index = ~np.isnan(disperRaw)
    disperFitted[index] = beta[0] / countRiboMean[index] + beta[1]
    data.disperFitted = disperFitted

    return data
