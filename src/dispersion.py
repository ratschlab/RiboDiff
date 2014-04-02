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

    num = len(data.geneIDs)
    disperRaw = np.empty((num, 1))
    disperRaw.fill(np.nan)

    for i in range(num):
    #for i in range(1):
        response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
        explanatory = data.matrix
        librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])
        disper = 0.1
        for j in range(10):
            modNB  = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            #modNB = sm.NegativeBinomial(response, explanatory, loglike_method='nb2')
            result = modNB.fit()
            #print result.summary()
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
