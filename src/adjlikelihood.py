import numpy as np
from scipy.stats import nbinom

def adj_loglikelihood(disper, explanatory, response, yhat, sign):

    # The NB distribution has two shape parameters: n and p. n is the size parameter, 
    # namely the number of successful trails. It also can be calculated as 1 / dispersion. 
    # Must be strictly positive, need not be integer. p stands for the probability of 
    # success in each trial. p = size / (size + mu). mu is the estimated count, yhat.
    # The variance of NB is mu + disper * mu^2, or mu + mu^2 / size

    n = 1 / disper
    p = n / (n + yhat)
    loglik = sum(nbinom.logpmf(response, n, p))

    mu = yhat

    #var = mu + disper * mu ** 2
    #w = 1 / ((1 / mu) ** 2 * var)
    #w = w.reshape(len(w), 1)
    #q, r = np.linalg.qr(explanatory * w ** 0.5)
    #coxreid = sum(np.log(abs(np.diag(r)[range(np.linalg.matrix_rank(explanatory))])))

    diagVec = mu / (1 + mu * disper)
    diagWM = np.diag(diagVec)
    xtwx = np.dot(np.dot(np.transpose(explanatory), diagWM), explanatory)
    coxreid = 0.5 * np.log(np.linalg.det(xtwx))

    return (loglik - coxreid) * sign

