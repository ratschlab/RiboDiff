#!/usr/bin/env python

""" This script generates two groups of negative binomial distributed discrete random variables. The two groups (condition A and condition B) contain equal number of  
    entries (GENE0001, GENE0002, ... GENE1000, ...), and each entry can have multiple replicates in each condition. The mean of each entry across replicates is negative 
    binomial distributed as well. By using arguments '--beta1' and '--beta2', users can change the dispersion-mean relationship in order to generate different data set, 
    for instance, Ribo-Seq and RNA-Seq data. """

import sys
import numpy as np
from scipy.stats import gamma
from scipy.stats import nbinom
import random

def parse_options(argv):

    from optparse import OptionParser, OptionGroup

    parser = OptionParser()

    required = OptionGroup(parser, 'REQUIRED')

    required.add_option('--numEntry', action='store', type='int', dest='numEntry', help='How many entries (genes / proteins) to be generated.')
    required.add_option('--numSampleConA', action='store', type='int', dest='numSampleConA', help='How many samples / replicates to be generated for condition A.')
    required.add_option('--numSampleConB', action='store', type='int', dest='numSampleConB', help='How many samples / replicates to be generated for condition B.')
    required.add_option('--output', action='store', type='string', dest='output', help='Output file name.')

    optional = OptionGroup(parser, 'OPTIONAL')

    optional.add_option('--nParamNB', action='store', type='float', dest='nParamNB', default=1, help='Assume the mean read count of each entry follows negative binomial distribution. ' \
                        'This argument is for the n parameter. [default: 1]')
    optional.add_option('--pParamNB', action='store', type='float', dest='pParamNB', default=0.01, help='Assume the mean read count of each entry follows negative binomial distribution. ' \
                        'This argument is for the p parameter. It controls the read count level. For instance, use 0.001 for RNA-Seq, and 0.001~0.01 for Ribo-Seq. [default: 0.01]')
    optional.add_option('--beta1', action='store', type='float', dest='beta1', default=0.100, help='Assume disper = beta1 / MeanCount + beta2 (beta1 must > 0.0). [default: 0.1]')
    optional.add_option('--beta2', action='store', type='float', dest='beta2', default=0.001, help='Assume disper = beta1 / MeanCount + beta2 (beta2 must > 0.0). [default: 0.001]')
    optional.add_option('--dispFile', action='store', type='string', dest='dispFile', help='A text file contains dispersions from which the read count will be generated for every gene. ' \
                        'The number of dispersions should equal to the number of entries. One dispersion per line in this text file. If this argument is given, the \'beta1\' and \'beta2\' ' \
                        'arguments will not be used.')
    optional.add_option('--addDisperError', action='store', type='float', dest='addDisperError', help='Add standard Gaussian distributed noise to log(dispersion) and use this dispersion ' \
                        'to generate read count. The value of this argument is the standard deviation of the Gaussian distribution. [defualt: no noise]')
    optional.add_option('--numDiff', action='store', type='int', dest='numDiff', help='How many entries show different mean read count. This script generates half entries as increased ' \
                        'mean count and another half as decreased mean count. For example, 1001 out of 20000 genes showing different fold change, the program will generate 500 showing ' \
                        'increase and 501 showing decrease. [default: no gene showing different]')
    optional.add_option('--diffFile', action='store', type='string', dest='diffFile', help='A text file contains integers 1, 2 or 0 to indicate which entries should be generated with mean ' \
                        'count UP, DOWN or NO CHANGE. The number of integers should equal to the number of entries. One integer per line in this text file. The program will randomly ' \
                        'choose which entries to alter their mean counts, if this argument is not given.')
    optional.add_option('--shapeGamma', action='store', type='float', dest='shapeGamma', default=1.5, help='Assume fold change of mean read count between two conditions follows gamma ' \
                        'distribution. Use this argument to set the shape paramter. [default: 1.5]')
    optional.add_option('--scaleGamma', action='store', type='float', dest='scaleGamma', default=0.5, help='Assume fold change of mean read count between two conditions follows gamma ' \
                        'distribution. Use this argument to set the scale paramter. [default: 0.5]')

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

    # First generate the mean read count for each gene. Assume this mean value follows NB distribution (Observed from real data).
    mu = nbinom.rvs(nParamNB, pParamNB, loc=0.0, size=numGene)

    # If the mean of certain genes are 0, change them as 1.
    idx = np.nonzero(mu == 0.0)[0]
    mu[idx] = 1.0

    # Generate dispersions for all genes.
    if not options.dispFile:
        # Generate dispersions as a function of mean count for all genes.
        disper = beta1 / mu + beta2
    else:
        # Load the dispersions to generate the count.
        disper = np.loadtxt(options.dispFile, dtype=float, skiprows=0, usecols=(0,))
        if disper.size != numGene:
            sys.stderr.write('\nError: The number of specified dispersions is not the same with number of genes!\n\n')
            sys.exit()

    # Add Gaussian distributed noise to log(dispersion).
    if options.addDisperError:
        std = options.addDisperError
        errorNorm = norm.rvs(loc=0.0, scale=std, size=numGene)
        disper = np.exp(np.log(disper) + errorNorm)

    muA = mu.copy()
    muB = mu.copy()

    # For some genes, generate read count with different mean value in different conditions.
    if options.numDiff or options.diffFile:

        # Fold change genes are randomly selected, or are chosen as indicated by file.
        if not options.diffFile:
            numDiff = options.numDiff

            # The number of genes showing increased and decreased mean count is equal or 1 less. 
            numDiffUp = numDiff / 2
            numDiffDn = numDiff - numDiffUp

            idx = random.sample(range(numGene), numDiff)
            idxUp = random.sample(idx, numDiffUp)
            idxDn = np.setdiff1d(idx, idxUp)
        else:
            diffInfo = np.loadtxt(options.diffFile, dtype=int, skiprows=0, usecols=(0,))
 
            idxUp = (diffInfo==2).nonzero()[0]
            idxDn = (diffInfo==1).nonzero()[0]

            numDiffUp = idxUp.size
            numDiffDn = idxDn.size

            numDiff = numDiffUp + numDiffDn

        if numDiff > options.numEntry:
            print 'numDiff should be smaller than numGene!'
            sys.exit()

        shapeParam = options.shapeGamma
        scaleParam = options.scaleGamma

        # Assume fold changes of mean count of different genes follow gamma distribution. If fold change value of increased gene set is x, the decreased set is 1/x.
        foldDiffUp = gamma.rvs(a=shapeParam, scale=scaleParam, loc=1.0, size=numDiffUp)
        foldDiffDn = 1.0 / gamma.rvs(a=shapeParam, scale=scaleParam, loc=1.0, size=numDiffDn)

        # Change the mean count of condition A and condition B without changing the overall mean count across the two conditions.
        # (MeanCountA + MeanCountB) / 2 = MeanCountOrigin & MeanCountA * FoldChange = MeanCountB
        # Assume there is a negative correlation between mean count and fold change.
        idxUpMem = np.searchsorted(np.sort(mu[idxUp]), mu[idxUp])
        idxDnMem = np.searchsorted(np.sort(mu[idxDn]), mu[idxDn])

        muAnewUp = 2 * np.sort(mu[idxUp]) / (np.sort(foldDiffUp)[::-1] + 1)
        muBnewUp = muAnewUp * np.sort(foldDiffUp)[::-1]
        muA[idxUp] = muAnewUp[idxUpMem]
        muB[idxUp] = muBnewUp[idxUpMem]

        muAnewDn = 2 * np.sort(mu[idxDn]) / (np.sort(foldDiffDn) + 1)
        muBnewDn = muAnewDn * np.sort(foldDiffDn)
        muA[idxDn] = muAnewDn[idxDnMem]
        muB[idxDn] = muBnewDn[idxDnMem]

    n = 1.0 / disper
    pA = n / (n + muA)
    pB = n / (n + muB)

    numDigits = len(str(numGene))
    with open(output, 'w') as FileOut:
        FileOut.write('Entry\t' + 'conditionA\t'*numSampleConA + 'conditionB\t'*numSampleConB + 'Dispersion\t' + 'MeanCondA\t' + 'MeanCondB\t' + 'MeanFoldChange\t' + 'SetAsDiff\n')
        for i in range(numGene):
            z = numDigits - len(str(i+1))
            name = 'G' + '0'*z + str(i+1)

            # The dispersion parameter (1/n) is the same for both conditions. The probability parameters are different if there is fold change in mean count for different conditions.
            countListA = nbinom.rvs(n[i], pA[i], size=numSampleConA).tolist()
            countListB = nbinom.rvs(n[i], pB[i], size=numSampleConB).tolist()
            countList  = countListA + countListB

            countString = '\t'.join(str(element) for element in countList)

            if not options.numDiff:
                setAsDiff = '-1'
            elif i in idxUp:
                setAsDiff = '1'
            elif i in idxDn:
                setAsDiff = '2'
            else:
                setAsDiff = '0'

            FileOut.write(name + '\t' + countString + '\t' + str(disper[i]) + '\t' + str(np.mean(countListA)) + '\t' + str(np.mean(countListB)) + '\t' + str((np.mean(countListB)+1e-5)/(np.mean(countListA)+1e-5)) + '\t' + setAsDiff + '\n')

if __name__ == '__main__':

    options = parse_options(sys.argv)
    generate_count(options)
    print '%s: Done!' % options.output

