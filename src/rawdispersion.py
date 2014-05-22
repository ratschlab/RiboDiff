import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
import adjlikelihood as al

def disper_raw(data):

    cntCutoff = 10

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    disperRaw = np.empty((num, 1))
    disperRaw.fill(np.nan)

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    for i in range(num):

        if i % 50 == 0:
            print '%i genes finished...' % i
        if i+1 == num:
            print '%i genes finished...' % num

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
                res  = minimize_scalar(al.adj_loglikelihood, bounds=(0, 1), args=(explanatory, response, yhat, sign), tol=1e-5, method='Bounded')
                disper = res.x
                if abs(np.log(disper) - np.log(disperBef)) < 0.03:
                    break
            disperRaw[i] = disper

    data.disperRaw = disperRaw

    return data

