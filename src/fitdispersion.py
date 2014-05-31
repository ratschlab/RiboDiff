import numpy as np
import statsmodels.api as sm

def disper_fit(data):

    countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countRiboMean=np.reshape(countRiboMean, (len(countRiboMean), 1))
    disperRaw = data.disperRaw
    #disperRaw = disperRaw.flatten()
    beta = np.array([0.1, 1])

    iter = 10
    for i in range(iter):
        ratio = disperRaw / (beta[0] + beta[1]/countRiboMean)
        idx = np.logical_and(ratio>1e-4, ratio<15)

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

        if sum(np.log(beta / betaBef)**2) < 1e-6:
            break

        if i == iter - 1:
            print 'Fitting dispersion does not converge.'

    disperFitted = disperRaw.copy()
    index = ~np.isnan(disperRaw)
    disperFitted[index] = beta[0] / countRiboMean[index] + beta[1]
    data.disperFitted = disperFitted

    print '*'*25
    print 'Fit dispersion: Done.'
    print '*'*25

    return data

