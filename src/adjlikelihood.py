import numpy as np
from scipy.stats import nbinom

def adj_loglikelihood(disper, explanatory, response, yhat, sign):

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

