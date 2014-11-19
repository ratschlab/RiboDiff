#!/usr/bin/env python

import sys
import os
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser, OptionGroup
import pdb

def empDisp_scatter(data, fileOutName):

    cntRiboNorm = data.countRibo / data.libSizesRibo
    cntRnaNorm  = data.countRna  / data.libSizesRna 

    idx = np.logical_and(np.sum(cntRiboNorm, axis=1)/data.libSizesRibo.size > 1, np.sum(cntRnaNorm, axis=1)/data.libSizesRna.size > 1).nonzero()[0]

    cntRiboMean = np.mean(cntRiboNorm[idx], axis=1)
    cntRnaMean  = np.mean(cntRnaNorm[idx],  axis=1)

    varRibo = np.var(cntRiboNorm[idx], axis=1, ddof=0)
    varRna  = np.var(cntRnaNorm[idx],  axis=1, ddof=0)

    dispRibo = (varRibo - cntRiboMean) / cntRiboMean ** 2
    dispRna  = (varRna  - cntRnaMean ) / cntRnaMean  ** 2

    stdDispRibo = np.std(np.log10(dispRibo[dispRibo > 0]), ddof=1)
    stdDispRna  = np.std(np.log10(dispRna[dispRna > 0]), ddof=1)

    if np.percentile(np.log2(cntRnaMean), 99.0) >= np.percentile(np.log2(cntRiboMean), 99.0):
        maxCnt = np.percentile(np.log2(cntRnaMean), 99.0)
    else:
        maxCnt = np.percentile(np.log2(cntRiboMean), 99.0)
    if np.percentile(np.log2(cntRnaMean), 1.0) <= np.percentile(np.log2(cntRiboMean), 1.0):
        minCnt = np.percentile(np.log2(cntRnaMean), 1.0)
    else:
        minCnt = np.percentile(np.log2(cntRiboMean), 1.0)

    winSize = (maxCnt - minCnt) / 15.0
    dispWinMedRibo = []
    dispWinMedRna  = []
    cntWinRibo  = []
    cntWinRna   = []

    for i in np.arange(minCnt, maxCnt, winSize):
        IDX1 = np.logical_and(np.logical_and(np.log2(cntRiboMean) > i, np.log2(cntRiboMean) < i + winSize), dispRibo > 0).nonzero()[0]
        IDX2 = np.logical_and(np.logical_and(np.log2(cntRnaMean)  > i, np.log2(cntRnaMean)  < i + winSize), dispRna  > 0).nonzero()[0]
        if i + winSize / 2.0 >= np.percentile(np.log2(cntRiboMean), 1.0) and i + winSize / 2.0 <= np.percentile(np.log2(cntRiboMean), 99.0):
            dispWinMedRibo.extend([np.median(dispRibo[IDX1])])
            cntWinRibo.extend([i + winSize / 2.0])
        if i + winSize / 2.0 >= np.percentile(np.log2(cntRnaMean), 1.0) and i + winSize / 2.0 <= np.percentile(np.log2(cntRnaMean), 99.0):
            dispWinMedRna.extend([np.median(dispRna[IDX2])])
            cntWinRna.extend([i + winSize / 2.0])

    fig, ax = plt.subplots()

    ax.scatter(np.log2(cntRnaMean[dispRna   > 0]), np.log10(dispRna[dispRna   > 0]), marker='o', color='lightsalmon',  s=0.5, lw=0, label='RNA-Seq' )
    ax.scatter(np.log2(cntRiboMean[dispRibo > 0]), np.log10(dispRibo[dispRibo > 0]), marker='o', color='lightskyblue', s=0.5, lw=0, label='Ribo-Seq')

    ax.plot(cntWinRna,  np.log10(dispWinMedRna),  color='tomato', linestyle='-', marker='^', markersize=4, markeredgewidth=0, label='RNA-Seq, trend' )
    ax.plot(cntWinRibo, np.log10(dispWinMedRibo), color='blue',   linestyle='-', marker='*', markersize=5, markeredgewidth=0, label='Ribo-Seq, trend')

    smallestDisp = min(np.hstack([np.log10(dispRna[dispRna > 0]), np.log10(dispRibo[dispRibo > 0])]))
    largestDisp  = max(np.hstack([np.log10(dispRna[dispRna > 0]), np.log10(dispRibo[dispRibo > 0])]))

    lowerBound = np.floor(smallestDisp) - 3
    upperBound = np.ceil(largestDisp) + 4

    ax.scatter(np.log2(cntRnaMean[dispRna   <= 0]), np.repeat(lowerBound + 1.0, cntRnaMean[dispRna   <= 0].size), marker='o', color='lightsalmon',  s=0.5, lw=0)
    ax.scatter(np.log2(cntRiboMean[dispRibo <= 0]), np.repeat(lowerBound + 0.8, cntRiboMean[dispRibo <= 0].size), marker='o', color='lightskyblue', s=0.5, lw=0)

    if np.mod(np.floor(smallestDisp), 2) == 1:
        lowerEndTick = np.floor(smallestDisp) - 1
    else:
        lowerEndTick = np.floor(smallestDisp)

    if np.mod(upperBound, 2) == 1:
        upperEndTick = upperBound - 1
    else:
        upperEndTick = upperBound

    ax.set_ylim(lowerBound, upperBound)
    ax.set_xlim(0, None)

    ax.spines['left'].set_visible(False)
    breakPoint = lowerEndTick + 1
    ax.plot((0, 0), (breakPoint+0.15, upperBound), color='black', lw=1.5)
    ax.plot((0, 0), (breakPoint-0.1,  lowerBound), color='black', lw=1.5)
    ax.plot((-0.1, 0.1), (breakPoint, breakPoint+0.2), color='black', lw=1, clip_on=False)
    ax.plot((-0.1, 0.1), (breakPoint-0.2, breakPoint), color='black', lw=1, clip_on=False)

    lowerEndTick = lowerEndTick.astype(int)
    upperEndTick = upperEndTick.astype(int)
    plt.yticks(np.arange(lowerEndTick, upperEndTick+1, 2))
    tklabels = ax.axes.get_yticks().tolist()
    tklabels[0] = '-Inf'
    ax.axes.set_yticklabels(tklabels)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    ax.legend(loc='upper right', prop={'size':9})
    ax.set_xlabel(r'$log_{2}(mean\/counts)$', fontsize=15)
    ax.set_ylabel(r'$log_{10}(dispersion)$', fontsize=15)
    ax.set_title('Empirical dispersion')

    ax.text(0.03, 0.96, r'$\sigma_{Ribo\/}=%1.2f$' % stdDispRibo, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=11)
    ax.text(0.03, 0.92, r'$\sigma_{RNA}=%1.2f$' % stdDispRna, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=11)

    plt.savefig(fileOutName, format='pdf')

def empDisp_hist(data, fileOutName):

    cntRiboNorm = data.countRibo / data.libSizesRibo
    cntRnaNorm  = data.countRna  / data.libSizesRna 

    idx = np.logical_and(np.sum(cntRiboNorm, axis=1)/data.libSizesRibo.size > 1, np.sum(cntRnaNorm, axis=1)/data.libSizesRna.size > 1).nonzero()[0]

    cntRiboMean = np.mean(cntRiboNorm[idx], axis=1)
    cntRnaMean  = np.mean(cntRnaNorm[idx],  axis=1)

    varRibo = np.var(cntRiboNorm[idx], axis=1, ddof=0)
    varRna  = np.var(cntRnaNorm[idx],  axis=1, ddof=0)

    dispRibo = (varRibo - cntRiboMean) / cntRiboMean ** 2
    dispRna  = (varRna  - cntRnaMean ) / cntRnaMean  ** 2

    index = np.logical_and(dispRibo > 0, dispRna > 0).nonzero()[0]
    deltaLogDisps = np.log10(dispRibo[index]) - np.log10(dispRna[index])

    num = index.size
    meanDelta = np.mean(deltaLogDisps)
    stdDelta  = np.std(deltaLogDisps, ddof=1)

    fig, ax = plt.subplots()

    minDelta = min(deltaLogDisps)
    maxDelta = max(deltaLogDisps)
    lowerBound = np.floor(minDelta)
    upperBound = np.ceil(maxDelta)
    winSize = (upperBound - lowerBound) / 40.0
    ax.hist(deltaLogDisps, np.arange(lowerBound, upperBound, winSize), histtype='bar', align='mid', color='limegreen', lw=0.5, edgecolor='white')

    if lowerBound < 0.0 and upperBound > 0.0:
        ax.set_xlim(lowerBound, upperBound)
    elif lowerBound >= 0.0:
        ax.set_xlim(-0.1 * (upperBound - lowerBound), upperBound)
    elif maxDelta <= 0.0:
        ax.set_xlim(lowerBound, 0.1 * (upperBound - lowerBound))
    ymin, ymax = plt.ylim()
    ax.plot((0, 0), (ymax * 0.02, ymax * 0.98), linestyle='--', color='orange', lw=1)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    ax.set_xlabel(r'$log_{10}\alpha_{Ribo}-log_{10}\alpha_{RNA}$', fontsize=15)
    ax.set_title(r'Distribution of dispersion difference')

    ax.text(0.88, 0.96, r'$n=%i$' % num, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=11)
    ax.text(0.88, 0.92, r'$\mu=%1.2f$' % meanDelta, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=11)
    ax.text(0.88, 0.88, r'$\sigma=%1.2f$' % stdDelta, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=11)

    plt.savefig(fileOutName, format='pdf')

def dotplot_disp(countMean, dispRaw, dispFitted, dispAdj, beta, fileOutName, whatSeq):

    fig, ax = plt.subplots()
    idx = np.logical_and(countMean>1, ~np.isnan(dispAdj))
    ax.scatter(np.log2(countMean[idx]), np.log10(dispRaw[idx]), marker='o', color='silver', s=4, lw=0, label='Raw')

    ax.scatter(np.log2(countMean[idx]), np.log10(dispAdj[idx]), marker='o', color='mediumturquoise', s=4, lw=0, label='Adjusted')

    countMeanSorted = np.sort(countMean[idx])
    dispFittedCal = beta[0]/countMeanSorted + beta[1]
    ax.plot(np.log2(countMeanSorted[dispFittedCal > 0]), np.log10(dispFittedCal[dispFittedCal > 0]), color='tomato', label='Fitted')

    ax.legend(loc='upper right', prop={'size':10})

    #xLowBound = 0.0
    #xUpBound  = np.percentile(countMean[idx], 95.0)
    #yLowBound = np.percentile(dispRaw[idx], 2.5)
    #yUpBound  = np.percentile(dispRaw[idx], 97.5)
    #ax.set_xlim(xLowBound, xUpBound)
    #ax.set_ylim(yLowBound, yUpBound)

    if whatSeq:
        ax.set_title('Dispersion estimation for %s' % whatSeq)
    else:
        ax.set_title('Dispersion estimation')
    plt.xlabel('Mean count')
    plt.ylabel('Dispersion')

    plt.savefig(fileOutName, format='pdf')

def dotplot_TEchange(countMean, logFoldChangeTE, dispAdj, padj, threshold, fileOutName):

    index = ~np.isnan(dispAdj)
    countMeanEff = countMean[index]
    logFoldChangeTeEff = logFoldChangeTE[index]
    padjEff = padj[index]

    idx = padjEff < threshold
    countMeanSig = countMeanEff[idx]
    logFoldChangeTeSig = logFoldChangeTeEff[idx]

    fig, ax = plt.subplots()
    ax.scatter(countMeanEff, logFoldChangeTeEff, marker='o', color='silver', s=3, lw=0, label='Tested genes')
    ax.scatter(countMeanSig, logFoldChangeTeSig, marker='o', color='orange', s=3, lw=0, label='Significant genes')

    ax.legend(loc='upper right', prop={'size':10})

    ax.set_title('Overview of TE fold change')
    plt.xlabel('Mean count')
    plt.ylabel('log2 TE fold change')

    xLowBound = (np.percentile(countMeanSig, 97.5) - np.min(countMeanSig)) * -0.02
    xUpBound  = np.percentile(countMeanSig, 97.5)

    ax.set_xlim(xLowBound, xUpBound)
    ax.set_ylim(-5, 5)

    plt.savefig(fileOutName, format='pdf')

def make_plots(data, opts):

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        outputNamePrefix = opts.outFile
    else:
        pos = opts.outFile.rfind('.')
        outputNamePrefix = opts.outFile[:pos]

    if opts.dispDiff:
        fileOutNameBoxplot = outputNamePrefix + '.disp.boxplot.pdf'
        boxplot_disp(data.dispAdjRibo, data.dispAdjRna, fileOutNameBoxplot)

    if opts.dispDiff:
        countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        countRiboMean = np.reshape(countRiboMean, (countRiboMean.size, 1))
        fileOutNameDotplotRibo = outputNamePrefix + '.dispRibo.dotplot.pdf'
        dotplot_disp(countRiboMean, data.dispRawRibo, data.dispAdjRibo, data.betaRibo, fileOutNameDotplotRibo, 'Ribo-Seq')

        countRnaMean = np.mean(data.countRna / data.libSizesRna, axis=1)
        countRnaMean = np.reshape(countRnaMean, (countRnaMean.size, 1))
        fileOutNameDotplotRna = outputNamePrefix + '.dispRNA.dotplot.pdf'
        dotplot_disp(countRnaMean, data.dispRawRna, data.dispAdjRna, data.betaRna, fileOutNameDotplotRna, 'RNA-Seq')
    else:
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        countMean = np.reshape(countMean, (countMean.size, 1))
        fileOutNameDotplot = outputNamePrefix + '.disp.dotplot.pdf'
        dotplot_disp(countMean, data.dispRaw, data.dispAdj, data.beta, fileOutNameDotplot, '')

    countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countMean = np.reshape(countMean, (countMean.size, 1))
    logFoldChangeTE = data.logFoldChangeTE
    if opts.dispDiff:
        dispAdj = data.dispAdjRibo
    else:
        dispAdj = data.dispAdj
    padj = data.padj
    threshold = 0.1
    fileOutNameDotplotTE = outputNamePrefix + '.TEchange.dotplot.pdf'
    dotplot_TEchange(countMean, logFoldChangeTE, dispAdj, padj, threshold, fileOutNameDotplotTE)

def parse_options(argv):

    parser = OptionParser()

    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-k', action='store', type='string', dest='dataPkl', help='Binary file containing all source data and results of TE change analysis.')
    required.add_option('-o', action='store', type='string', dest='outputPrefix', help='Specify the prefix of the output file name.')

    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-p', action='store', type='string', dest='plotWhich', default='All', help='Which figure to be plotted. Options: EmpDisp, EstDisp, TE or All. [default: All]')

    parser.add_option_group(required)
    parser.add_option_group(optional)

    (opts, args) = parser.parse_args(argv)

    if len(argv) < 2:
        parser.print_help()
        sys.exit()

    mandatories = ['dataPkl', 'outputPrefix']
    for eachOpt in mandatories:
        if not opts.__dict__[eachOpt]:
            parser.error('-%s is a required option.\n' % eachOpt[0])
        if eachOpt in mandatories[:1] and not os.path.exists(opts.__dict__[eachOpt]):
            sys.stderr.write('\nError: File \'%s\' does not exist.\n\n' % opts.__dict__[eachOpt])
            sys.exit()

        if eachOpt == mandatories[1]:
            if not os.path.dirname(opts.__dict__[eachOpt]):
                opts.__dict__[eachOpt] = os.getcwd() + os.getcwd()[0] + opts.__dict__[eachOpt]
            if not os.path.exists(os.path.dirname(opts.__dict__[eachOpt])):
                sys.stderr.write('\nError: Path \'%s\' of the output file does not exist.\n\n' % os.path.dirname(opts.__dict__[eachOpt]))
                sys.exit()

    if opts.plotWhich not in ['EmpDisp', 'EstDisp', 'TE', 'All']:
        parser.error('-p option can only take \'EmpDisp\', \'EstDisp\', \'TE\' or \'All\' as argument.\n')

    return opts

if __name__ == '__main__':

    opts = parse_options(sys.argv)

    with open(opts.dataPkl, 'rb') as FileIn:
        data = pickle.load(FileIn)

    if opts.plotWhich in ['EmpDisp', 'All']:
        fileOutName = opts.outputPrefix + '.EmpDisp.scatter.pdf'
        empDisp_scatter(data, fileOutName)

        fileOutName = opts.outputPrefix + '.EmpDisp.hist.pdf'
        empDisp_hist(data, fileOutName)

#        if opts.plotWhich in ['DotDisp', 'All']:
#            countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
#            countRiboMean = np.reshape(countRiboMean, (countRiboMean.size, 1))
#            fileOutNameDotplotRibo = opts.outputPrefix + '.dispRibo.dotplot.pdf'
#            dotplot_disp(countRiboMean, data.dispRawRibo, data.dispFittedRibo, data.dispAdjRibo, data.betaRibo, fileOutNameDotplotRibo, 'Ribo-Seq')
#
#            countRnaMean = np.mean(data.countRna / data.libSizesRna, axis=1)
#            countRnaMean = np.reshape(countRnaMean, (countRnaMean.size, 1))
#            fileOutNameDotplotRna = opts.outputPrefix + '.dispRNA.dotplot.pdf'
#            dotplot_disp(countRnaMean, data.dispRawRna, data.dispFittedRna, data.dispAdjRna, data.betaRna, fileOutNameDotplotRna, 'RNA-Seq')
#    else:
#        if opts.plotWhich in ['DotDisp', 'All']:
#            countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
#            countMean = np.reshape(countMean, (countMean.size, 1))
#            fileOutNameDotplot = opts.outputPrefix + '.disp.dotplot.pdf'
#            dotplot_disp(countMean, data.dispRaw, data.dispFitted, data.dispAdj, data.beta, fileOutNameDotplot, '')
#
#    if opts.plotWhich in ['DotTE', 'All']:
#        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
#        countMean = np.reshape(countMean, (countMean.size, 1))
#        logFoldChangeTE = data.logFoldChangeTE
#        if data.dispDiff:
#            dispAdj = data.dispAdjRibo
#        else:
#            dispAdj = data.dispAdj
#        padj = data.padj
#        threshold = 0.1
#        fileOutNameDotplotTE = opts.outputPrefix + '.TEchange.dotplot.pdf'
#        dotplot_TEchange(countMean, logFoldChangeTE, dispAdj, padj, threshold, fileOutNameDotplotTE)

