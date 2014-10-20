#!/usr/bin/env python

import sys
import numpy as np
from scipy.stats import norm
from scipy.stats import nbinom
import random
import pdb

def parse_options(argv):

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()

    required = OptionGroup(parser, 'REQUIRED')

    required.add_option('--numEntry', action='store', type='int', dest='numEntry', help='How many entries (genes / proteins) to be generated.')
    required.add_option('--numSampleConA', action='store', type='int', dest='numSampleConA', help='How many samples / replicates to be generated for condition A.')
    required.add_option('--numSampleConB', action='store', type='int', dest='numSampleConB', help='How many samples / replicates to be generated for condition B.')
    required.add_option('--output', action='store', type='string', dest='output', help='Output file name.')

    optional = OptionGroup(parser, 'OPTIONAL')

    optional.add_option('--nParamNB', action='store', type='float', dest='nParamNB', default=1, help='Assume the mean read count follows negative binomial distribution. This argument is for the n parameter. [default: 1]')
    optional.add_option('--pParamNB', action='store', type='float', dest='pParamNB', default=0.01, help='Assume the mean read count follows negative binomial distribution. This argument is for the p parameter. This parameter controls the read count level. Usually use 0.01 for Ribo-seq, 0.001 for RNA-seq. [default: 0.01]')
    optional.add_option('--beta1', action='store', type='float', dest='beta1', default=0.200, help='Assume disper = beta1 / MeanCount + beta2 (beta1 must > 0.0). [default: 0.2]')
    optional.add_option('--beta2', action='store', type='float', dest='beta2', default=0.001, help='Assume disper = beta1 / MeanCount + beta2 (beta2 must > 0.0). [default: 0.001]')
    optional.add_option('--addDisperError', action='store', type='float', dest='addDisperError', help='Add standard Gaussian distributed noise to log(dispersion) and use this dispersion to generate read count. The value of this argument is the standard deviation of the Gaussian distribution. [defualt: no noise added to dispersion]')
    optional.add_option('--numDiff', action='store', type='int', dest='numDiff', help='How many entries show different mean read count. This script generates half as increased mean count and half as decreased mean count. For example, 1001 out of 2000 genes showing different fold change, the program will generate 501 showing increase and 499 showing decrease. [default: no gene showing different]')
    optional.add_option('--meanFoldDiff', action='store', type='float', dest='meanFoldDiff', help='The MEAN fold change of mean read count in condition B compare to condition A. You only need to specify either increased or decreased mean fold change, the program will calculate the other. For example, by specifying 1.5 fold change, 501 genes have mean fold change equal to 3/2, and 499 genes have mean fold change equal to 2/3.  Fold change values for genes are Gaussian distributed. [defult: no gene showing different]')

    parser.add_option_group(required)
    parser.add_option_group(optional)

    (options, args) = parser.parse_args(argv)

    if len(argv) < 2:
        parser.print_help()
        sys.exit()

    mandatories = ['numEntry', 'numSampleConA', 'numSampleConB', 'output']
    for eachM in mandatories:
        if not options.__dict__[eachM]:
            sys.stderr.write('\nError: --%s is a required option.\n\n' % eachM)
            sys.exit()

    if (options.numDiff and not options.meanFoldDiff) or (options.meanFoldDiff and not options.numDiff):
        sys.stderr.write('\n--numDiff and --foldDiff should be used together.\n\n')
        sys.exit()

    return options

def generate_count(options):

    numGene = options.numEntry
    numSampleConA = options.numSampleConA
    numSampleConB = options.numSampleConB
    nParamNB = options.nParamNB
    pParamNB = options.pParamNB
    beta1 = options.beta1
    beta2 = options.beta2
    output = options.output

    mu = nbinom.rvs(nParamNB, pParamNB, 0.0, numGene)
    idx = np.nonzero(mu == 0.0)[0]
    mu[idx] = 1
    disper = beta1 / mu + beta2

    if options.addDisperError:
        std = options.addDisperError
        errorNorm = norm.rvs(loc=0.0, scale=std, size=numGene)
        disper = np.exp(np.log(disper) + errorNorm)

    muA = mu.copy()
    muB = mu.copy()
    if options.numDiff and options.meanFoldDiff:
        if options.numDiff > options.numEntry:
            print 'numDiff should be smaller than numGene!'
            sys.exit()

        numDiff = options.numDiff
        numDiffUp = round(numDiff / 2.0)
        numDiffDn = numDiff - numDiffUp
        meanFoldDiffUp = options.meanFoldDiff
        meanFoldDiffDn = 1 / meanFoldDiffUp

        randomNormUp = norm.rvs(loc=0.0, scale=3**0.5, size=numDiffUp)
        probMax = norm.pdf(0.0, loc=0.0, scale=3**0.5)
        foldDiffUp = meanFoldDiffUp + randomNormUp / abs(randomNormUp) * (meanFoldDiffUp - 1) * (1 - norm.pdf(randomNormUp, loc=0.0, scale=3**0.5) / probMax)

        randomNormDn = norm.rvs(loc=0.0, scale=3**0.5, size=numDiffDn)
        probMax = norm.pdf(0.0, loc=0.0, scale=3**0.5)
        foldDiffDn = meanFoldDiffDn + randomNormDn / abs(randomNormDn) * (1 - meanFoldDiffDn) * (1 - norm.pdf(randomNormDn, loc=0.0, scale=3**0.5) / probMax)

        # fold change genes are randomly selected.
        idx = random.sample(range(numGene), numDiff)

        # fold change genes are drawn from the NB distribution as well.
        #idx = nbinom.rvs(nParamNB, pParamNB, 0.0, numDiff)

        foldDiff = np.hstack([foldDiffUp, foldDiffDn])
        muA[idx] = 2 * mu[idx] / (foldDiff + 1)
        muB[idx] = muA[idx] * foldDiff

    numDigits = len(str(numGene))
    with open(output, 'w') as FileOut:
        FileOut.write('Entry\t' + 'ConditionA\t'*numSampleConA + 'ConditionB\t'*numSampleConB + 'Dispersion\t' + 'MeanCondA\t' + 'MeanCondB\t' + 'MeanFoldChange\n')
        for i in range(numGene):
            z = numDigits - len(str(i+1))
            name = 'G' + '0'*z + str(i+1)
            n = 1.0 / disper[i]
            pA = n / (n + muA[i])
            pB = n / (n + muB[i])
            countListA = nbinom.rvs(n, pA, size=numSampleConA).tolist()
            countListB = nbinom.rvs(n, pB, size=numSampleConB).tolist()
            countList  = countListA + countListB
            countString = '\t'.join(str(element) for element in countList)
            FileOut.write(name + '\t' + countString + '\t' + str(disper[i]) + '\t' + str(np.mean(countListA)) + '\t' + str(np.mean(countListB)) + '\t' + str(np.mean(countListB)/np.mean(countListA)) + '\n') 

if __name__ == '__main__':

    options = parse_options(sys.argv)
    generate_count(options)
    print 'Done!'
