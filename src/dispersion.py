import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from scipy.stats import nbinom

def adj_loglikelihood_simp(disper, explanatory, response, yhat, sign):

    

def disper_simp(data, explanatory):

    num = len(data.geneIDs)
    #for i in range(num):
    for i in range(1):
        response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
        librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])
        disper = 0.1
        #for j in range(10):
        for j in range(1):
            modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            #modNB = sm.NegativeBinomial(response, explanatory, loglike_method='nb2')
            result = modNB.fit()
            print result.summary()
            disperBef = disper
            yhat = result.mu
            sign = -1.0
            minimize_scalar(adj_loglikelihood_simp, bounds=(0, 1), args=(param, explanatory, response, yhat, sign), tol=0.01, method='bounded')
            
