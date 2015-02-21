#!/usr/bin/env python
"""
Negatibe binomial likelihood function.
"""

import numpy as np
import scipy as sp
from scipy.stats import nbinom
from scipy.special import digamma

def adj_loglikelihood(xVec, lenSampleRibo, lenSampleRna, X, y, mu, sign):

    disp = np.hstack([np.repeat(xVec[0], lenSampleRibo), np.repeat(xVec[1], lenSampleRna)])
    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + np.dot(mu.transpose(), disp))
    diagWM = np.diagflat(diagVec)
    xtwx = np.dot(np.dot(np.transpose(X), diagWM), X)
    coxreid = 0.5 * np.log(np.linalg.det(xtwx))

    return (loglik - coxreid) * sign

def adj_loglikelihood_gradient(xVec, lenSampleRibo, lenSampleRna, X, y, mu, sign):

    disp = np.hstack([np.repeat(xVec[0], lenSampleRibo), np.repeat(xVec[1], lenSampleRna)])
    Gradient = np.zeros_like(xVec)

    for i in range(len(xVec)):
        f1 = (digamma((1 / xVec[i])) - digamma( (1 / xVec[i]) + y[i])) / (xVec[i] ** 2)
        f2 = -((xVec[i] * mu[i] + (1 + xVec[i] * mu[i]) * np.log(1 / (1 + xVec[i] * mu[i]))) / ((xVec[i] ** 2) * (1 + xVec[i] * mu[i])))
        f3 = y[i] / (xVec[i] + (xVec[i] ** 2) * mu[i])
        f4 = 0.5 * X.shape[1] * (mu[i] / (1 + np.dot(mu.transpose(), disp)))
        Gradient[i] = f1 + f2 + f3 + f4

    return Gradient

def adj_loglikelihood_scalar(disp, X, y, mu, sign):

    n = 1 / disp
    p = n / (n + mu)
    loglik = sum(nbinom.logpmf(y, n, p))

    diagVec = mu / (1 + mu * disp)
    diagWM = np.diag(diagVec)
    xtwx = np.dot(np.dot(np.transpose(X), diagWM), X)
    coxreid = 0.5 * np.log(np.linalg.det(xtwx))

    return (loglik - coxreid) * sign
