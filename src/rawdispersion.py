import numpy as np
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
import adjlikelihood as al
import matplotlib.pyplot as plt

def disper_raw(data):

    cntCutoff = 10

    print 'Start to estimate raw dispersions.'

    num = len(data.geneIDs)
    disperRaw = np.empty((num, 1))
    disperRaw.fill(np.nan)

    explanatory = data.matrix
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    ###create varibles for plotting
    emptyArray = np.empty((5, 10))
    emptyArray.fill(np.nan)
    dispersion = emptyArray.copy()
    deviance = emptyArray.copy()
    likelihood = emptyArray.copy()
    cnt = 0
    ###############################

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
                result = modNB.fit()
                #print result.summary()
                #print result.mu

                ###assign values to varibles for plotting.
                if cnt < 5:
                    dispersion[cnt, j] = disper
                    deviance[cnt, j] = result.deviance
                    likelihood[cnt, j] = result.llf
                ####################################################################

                disperBef = disper
                yhat = result.mu
                sign = -1.0
                res  = minimize_scalar(al.adj_loglikelihood, bounds=(0, 1), args=(explanatory, response, yhat, sign), tol=1e-5, method='Bounded')
                disper = res.x
                if abs(np.log(disper) - np.log(disperBef)) < 0.03:
                    break
            disperRaw[i] = disper

            ######plot the trends of changes of dispersion, deviance and likelihood.
            cnt += 1
            if cnt == 5:
                fig, ax = plt.subplots()
                for k in range(cnt):
                    ax.plot(np.arange(10), dispersion[k,], 'o-')
                    plt.xlabel('iteration #')
                    plt.ylabel('dispersion')
                plt.savefig('/Users/zhongyi/Yi/Work/RFSeq/src/plot_dispersion.pdf')
                
                fig, ax = plt.subplots()
                for k in range(cnt):
                    ax.plot(np.arange(10), deviance[k,], 'o-')
                    plt.xlabel('iteration #')
                    plt.ylabel('deviance')
                plt.savefig('/Users/zhongyi/Yi/Work/RFSeq/src/plot_deviance.pdf')

                fig, ax = plt.subplots()
                for k in range(cnt):
                    ax.plot(np.arange(10), likelihood[k,], 'o-')
                    plt.xlabel('iteration #')
                    plt.ylabel('loglikelihood')
                plt.savefig('/Users/zhongyi/Yi/Work/RFSeq/src/plot_likelihood.pdf')
                print 'Plotting: Done!'
                cnt == -1
            #######################################################################

    data.disperRaw = disperRaw

    return data

