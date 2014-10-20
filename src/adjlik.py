import numpy as np
import scipy as sp
from scipy.stats import nbinom
from scipy.special import digamma

def adj_loglikelihood(xVec, lenSampleRibo, lenSampleRNA, X, y, mu, sign):

    # The NB distribution has two shape parameters: n and p. n is the size parameter, 
    # namely the number of successful trails. It also can be calculated as 1 / dispersion. 
    # Must be strictly positive, need not be integer. p stands for the probability of 
    # success in each trial. p = size / (size + mu). mu is the estimated count.
    # The variance of NB is mu + disper * mu^2, or mu + mu^2 / size

    disper = np.hstack([np.repeat(xVec[0], lenSampleRibo), np.repeat(xVec[1], lenSampleRNA)])
    n = 1 / disper
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))
    #loglik = np.log(sp.special.gamma(y + n)/(sp.special.gamma(y+1)*sp.special.gamma(n)) * p**n * (1-p)**y)

    #var = mu + disper * mu ** 2
    #w = 1 / ((1 / mu) ** 2 * var)
    #w = w.reshape(len(w), 1)
    #q, r = np.linalg.qr(explanatory * w ** 0.5)
    #coxreid = sum(np.log(abs(np.diag(r)[range(np.linalg.matrix_rank(explanatory))])))
    
    diagVec = mu / (1 + np.dot(mu.transpose(), disper))
    diagWM = np.diagflat(diagVec)
    xtwx = np.dot(np.dot(np.transpose(X), diagWM), X)
    coxreid = 0.5 * np.log(np.linalg.det(xtwx))

    return (loglik - coxreid) * sign

def adj_loglikelihood_gradient(xVec, lenSampleRibo, lenSampleRNA, X, y, mu, sign):

    #disper = np.hstack([np.repeat(xVec[0], lenSampleRibo), np.repeat(xVec[1], lenSampleRNA)])
    #Gradient = np.zeros_like(disper)

    ##Iterate over the dimensions of the dispersion and compute the components of the gradient
    #for i in range(len(disper)):
    #    f1 = (digamma((1 / disper[i])) - digamma( (1 / disper[i]) + y[i])) / (disper[i] ** 2)
    #    f2 = -((disper[i] * mu[i] + (1 + disper[i] * mu[i]) * np.log(1 / (1 + disper[i] * mu[i]))) / ((disper[i] ** 2) * (1 + disper[i] * mu[i])))
    #    f3 = y[i] / (disper[i] + (disper[i] ** 2) * mu[i])
    #    f4 = 0.5 * X.shape[1] * (mu[i] / (1 + np.dot(mu.transpose(), disper)))
    #    Gradient[i] = f1 + f2 + f3 + f4
    #    #Gradient[i] = f1 + f2 + f3

    disper = np.hstack([np.repeat(xVec[0], lenSampleRibo), np.repeat(xVec[1], lenSampleRNA)])
    Gradient = np.zeros_like(xVec)

    #Iterate over the dimensions of the dispersion and compute the components of the gradient
    for i in range(len(xVec)):
        f1 = (digamma((1 / xVec[i])) - digamma( (1 / xVec[i]) + y[i])) / (xVec[i] ** 2)
        f2 = -((xVec[i] * mu[i] + (1 + xVec[i] * mu[i]) * np.log(1 / (1 + xVec[i] * mu[i]))) / ((xVec[i] ** 2) * (1 + xVec[i] * mu[i])))
        f3 = y[i] / (xVec[i] + (xVec[i] ** 2) * mu[i])
        f4 = 0.5 * X.shape[1] * (mu[i] / (1 + np.dot(mu.transpose(), disper)))
        Gradient[i] = f1 + f2 + f3 + f4
        #Gradient[i] = f1 + f2 + f3

    return Gradient

def adj_loglikelihood_scalar(disper, X, y, mu, sign):

    n = 1 / disper
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + mu * disper)
    diagWM = np.diag(diagVec)
    xtwx = np.dot(np.dot(np.transpose(X), diagWM), X)
    coxreid = 0.5 * np.log(np.linalg.det(xtwx))

    return (loglik - coxreid) * sign
