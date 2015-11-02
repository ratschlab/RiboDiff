#!/usr/bin/env python
"""
Plotting the data and results.
"""

import os
import sys
import cPickle as pickle
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser, OptionGroup

def parse_options(argv):

    parser = OptionParser()

    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-k', action='store', type='string', dest='dataPkl', help='Binary file containing all source data and results of TE change analysis.')
    required.add_option('-o', action='store', type='string', dest='outputPrefix', help='Specify the prefix of the output file name.')

    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-p', action='store', type='string', dest='plotWhich', default='All', help='Which figure to be plotted. Options: EmpDisp, TE or All. [default: All]')
    optional.add_option('-q', action='store', type='float', dest='cutoffFDR', default=0.1, help='Set the FDR cutoff for significant case to plot. [default: 0.1]')

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
                try:
                    os.makedirs(os.path.dirname(opts.__dict__[eachOpt]))
                except OSError:
                    sys.stderr.write('\nError: Failed to create directory: \'%s\' \n\n' % os.path.dirname(opts.__dict__[eachOpt]))
                    sys.exit()

    if opts.plotWhich not in ['EmpDisp', 'TE', 'All']:
        parser.error('-p option can only take \'EmpDisp\', \'TE\' or \'All\' as argument.\n')

    return opts

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
        if i + winSize / 2.0 >= np.percentile(np.log2(cntRiboMean), 2.5) and i + winSize / 2.0 <= np.percentile(np.log2(cntRiboMean), 97.5):
            dispWinMedRibo.extend([np.median(dispRibo[IDX1])])
            cntWinRibo.extend([i + winSize / 2.0])
        if i + winSize / 2.0 >= np.percentile(np.log2(cntRnaMean), 2.5) and i + winSize / 2.0 <= np.percentile(np.log2(cntRnaMean), 97.5):
            dispWinMedRna.extend([np.median(dispRna[IDX2])])
            cntWinRna.extend([i + winSize / 2.0])

    fig, ax = plt.subplots()

    ax.scatter(np.log2(cntRnaMean[dispRna   > 0]), np.log10(dispRna[dispRna   > 0]), marker='o', color='lightsalmon',  s=0.5, lw=0, label='dispersion, RNA-Seq' )
    ax.scatter(np.log2(cntRiboMean[dispRibo > 0]), np.log10(dispRibo[dispRibo > 0]), marker='o', color='lightskyblue', s=0.5, lw=0, label='dispersion, Ribo-Seq')

    ax.plot(cntWinRna,  np.log10(dispWinMedRna),  color='crimson', linestyle='-', marker='^', markersize=4, markeredgewidth=0, label='window mean, RNA-Seq' )
    ax.plot(cntWinRibo, np.log10(dispWinMedRibo), color='dodgerblue',   linestyle='-', marker='*', markersize=5, markeredgewidth=0, label='window mean, Ribo-Seq')

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

    ax.legend(loc='upper right', prop={'size':10})
    ax.set_xlabel(r'$log_{2}(mean\/counts)$', fontsize=15)
    ax.set_ylabel(r'$log_{10}(dispersion)$', fontsize=15)
    ax.set_title('Empirical Dispersion')

    ax.text(0.03, 0.96, r'$\sigma_{Ribo\/}=\,%1.2f$' % stdDispRibo, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=12)
    ax.text(0.03, 0.92, r'$\sigma_{RNA}=\,%1.2f$' % stdDispRna, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=12)

    plt.savefig(fileOutName, format='pdf', bbox_inches='tight')

def cnt_deltaTE_scatter(data, fdr, fileOutName):

    cntRiboNorm = data.countRibo / data.libSizesRibo
    cntRnaNorm  = data.countRna  / data.libSizesRna

    padj = data.padj.flatten()

    with np.errstate(invalid='ignore'):
        idx = np.logical_and(~np.isnan(padj), np.logical_and(np.sum(cntRiboNorm, axis=1)/data.libSizesRibo.size > 2, np.sum(cntRnaNorm, axis=1)/data.libSizesRna.size > 2)).nonzero()[0]
        cntRiboMean = np.mean(cntRiboNorm[idx], axis=1)
        logFoldChangeTE = data.logFoldChangeTE[idx]

        index = np.logical_and(padj < fdr, np.logical_and(np.sum(cntRiboNorm, axis=1)/data.libSizesRibo.size > 2, np.sum(cntRnaNorm, axis=1)/data.libSizesRna.size > 2)).nonzero()[0]
        cntRiboMeanSig = np.mean(cntRiboNorm[index], axis=1)
        logFoldChangeTEsig = data.logFoldChangeTE[index]

    fig, ax = plt.subplots()

    ax.scatter(cntRiboMean, logFoldChangeTE, marker='o', color='silver', s=1, lw=0, label='Tested genes')
    ax.scatter(cntRiboMeanSig, logFoldChangeTEsig, marker='o', color='darkorange', s=1, lw=0, label='Significant genes')

    ax.legend(loc='upper right', prop={'size':10})

    xLowerBound = (np.percentile(cntRiboMean, 99.0) - min(cntRiboMean)) * -0.02
    xUpperBound = np.percentile(cntRiboMean, 99.0)

    ax.set_xlim(xLowerBound, xUpperBound)
    ax.set_ylim(np.percentile(logFoldChangeTE, 0.5), np.percentile(logFoldChangeTE, 99.5))

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    ax.set_title(r'Translation Efficiency Change')
    ax.set_xlabel(r'$Mean\/count\/of\/Ribo$-$Seq$', fontsize=15)
    ax.set_ylabel(r'$log_{2}(TE_{%s}/TE_{%s})$' % (data.nameCondB, data.nameCondA), fontsize=15)

    plt.savefig(fileOutName, format='pdf', bbox_inches='tight')

def deltaTE_hist(data, fdr, fileOutName):

    cntRiboNorm = data.countRibo / data.libSizesRibo
    cntRnaNorm  = data.countRna  / data.libSizesRna

    padj = data.padj.flatten()
    
    deltaTE = data.logFoldChangeTE.flatten()

    with np.errstate(invalid='ignore'):
        idxNaN   = np.nonzero(~np.isnan(padj))[0]
        idxSigDn = np.nonzero(np.logical_and(padj <= fdr, deltaTE < 0))[0]
        idxSigUp = np.nonzero(np.logical_and(padj <= fdr, deltaTE > 0))[0]

    num = idxNaN.size
    muDeltaTE  = np.mean(deltaTE[idxNaN])
    stdDeltaTE = np.std(deltaTE[idxNaN], ddof=0)

    fig, ax = plt.subplots()

    maxExtreme = max(deltaTE)
    minExtreme = min(deltaTE)
    stepSize   = (np.percentile(deltaTE, 97.5) - np.percentile(deltaTE, 2.5)) / 25.0

    ax.hist(deltaTE[idxNaN], np.arange(minExtreme, maxExtreme, stepSize), histtype='bar', color='darkgrey', rwidth=1.0, linewidth=0.5, edgecolor='white', align='mid', label='All')
    ax.hist(deltaTE[idxSigDn], np.arange(minExtreme, maxExtreme, stepSize), histtype='bar', color='crimson', rwidth=1.0, linewidth=0.5, edgecolor='white', align='mid', label='TE down')
    ax.hist(deltaTE[idxSigUp], np.arange(minExtreme, maxExtreme, stepSize), histtype='bar', color='dodgerblue', rwidth=1.0, linewidth=0.5, edgecolor='white', align='mid', label='TE up')

    ax.legend(loc='upper right', prop={'size':10})

    ax.set_title(r'Histogram of Translation Efficiency Change')
    ax.set_xlabel(r'$log_{2}(TE_{%s}/TE_{%s})$' % (data.nameCondB, data.nameCondA), fontsize=15)
    ax.set_ylabel(r'$Frequency$', fontsize=15)

    ax.text(0.02, 0.97, r'$n\,=\,%i$' % num, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=12)
    ax.text(0.02, 0.93, r'$\mu\,=\,%1.2f$' % muDeltaTE, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=12)
    ax.text(0.02, 0.89, r'$\sigma\,=\,%1.2f$' % stdDeltaTE, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontsize=12)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    plt.savefig(fileOutName, format='pdf', bbox_inches='tight')

def make_plots(data, opts):

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        outputNamePrefix = opts.outFile
    else:
        pos = opts.outFile.rfind('.')
        if opts.outFile[pos:pos+2] == './':
            outputNamePrefix = opts.outFile
        else:
            outputNamePrefix = opts.outFile[:pos]

    fileOutName = outputNamePrefix + '.EmpDisp.scatter.pdf'
    empDisp_scatter(data, fileOutName)

    if opts.__dict__['cutoffFDR']:
        fdr = opts.cutoffFDR
    else:
        fdr = 0.1

    fileOutName = outputNamePrefix + '.TEchange.scatter.pdf'
    cnt_deltaTE_scatter(data, fdr, fileOutName)

    fileOutName = outputNamePrefix + '.TEchange.hist.pdf'
    deltaTE_hist(data, fdr, fileOutName)

if __name__ == '__main__':

    opts = parse_options(sys.argv)

    with open(opts.dataPkl, 'rb') as FileIn:
        data = pickle.load(FileIn)

    if opts.plotWhich in ['EmpDisp', 'All']:
        fileOutName = opts.outputPrefix + '.EmpDisp.scatter.pdf'
        empDisp_scatter(data, fileOutName)

    if opts.plotWhich in ['TE', 'All']:
        fdr = opts.cutoffFDR
        fileOutName = opts.outputPrefix + '.TEchange.scatter.pdf'
        cnt_deltaTE_scatter(data, fdr, fileOutName)

        fileOutName = opts.outputPrefix + '.TEchange.hist.pdf'
        deltaTE_hist(data, fdr, fileOutName)

