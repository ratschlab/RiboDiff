import numpy as np
import statsmodels.api as sm
import pdb

def do_fitting(data, obj):

    if obj == 'Ribo':
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        disperRaw = data.disperRawRibo
    elif obj == 'mRNA':
        countMean = np.mean(data.countRNA  / data.libSizesRNA,  axis=1)
        disperRaw = data.disperRawRibo
    elif obj == 'RR':
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        disperRaw = data.disperRaw
    else:
        pass

    disperRawMthd = data.disperRawMthd

    countMean = np.reshape(countMean, (len(countMean), 1))
    disperRawConv = data.disperRawConv
    beta = np.array([1, 0.1])

    iter = 20

    index = np.nonzero(disperRawConv == True)[0]

    for i in range(iter):

        ratioForBound = disperRaw[index] / (beta[0] / countMean[index] + beta[1])
        
        lowerBound = np.percentile(ratioForBound[~np.isnan(ratioForBound)], 5)
        upperBound = np.percentile(ratioForBound[~np.isnan(ratioForBound)], 95)

        ratio = disperRaw / (beta[0] / countMean + beta[1])
        idx = np.logical_and(ratio > lowerBound, ratio < upperBound)

        matrix = np.empty((len(np.nonzero(idx)[0]), 2))
        matrix.fill(np.nan)
        matrix[:, 0] = 1 / countMean[idx]
        matrix[:, 1] = 1

        modGamma = sm.GLM(disperRaw[idx], matrix, family=sm.families.Gamma(sm.families.links.identity))
        result = modGamma.fit(start_params=beta)

        betaBef = beta
        beta = result.params

        if sum(np.log(beta / betaBef)**2) < 1e-4:
            disperFittedConv = True
            break
        if i == iter - 1:
            disperFittedConv = False
            if obj != 'RR':
                print 'Fitting dispersion for %s does not converge.' % obj
            else:
                print 'Fitting dispersion does not converge.'

    disperFitted = disperRaw.copy()
    IDX = ~np.isnan(disperRaw)
    disperFitted[IDX] = beta[0] / countMean[IDX] + beta[1]
    disperFitted[disperFitted < 0] = disperRaw[disperFitted < 0]

    if obj == 'Ribo':
        data.disperFittedRibo = disperFitted
        data.betaRibo = beta
        data.disperFittedRiboConv = disperFittedConv
        data.disperFittedRiboIdx = idx
    elif obj == 'mRNA':
        data.disperFittedRNA = disperFitted
        data.betaRNA = beta
        data.disperFittedRNAConv = disperFittedConv
        data.disperFittedRNAIdx = idx
    elif obj == 'RR':
        data.disperFitted = disperFitted
        data.beta = beta
        data.disperFittedConv = disperFittedConv
        data.disperFittedIdx = idx
    else:
        pass

    if obj != 'RR':
        print 'Fit dispersion for %s: Done.' % obj
    else:
        print 'Fit dispersion: Done.'

    return data

def disper_fit(data, opts):

    if opts.dispDiff:
        data = do_fitting(data, obj='Ribo')
        print '*'*25
        data = do_fitting(data, obj='mRNA')
    else:
        data = do_fitting(data, obj='RR')

    return data
