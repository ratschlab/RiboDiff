#!/usr/bin/env python

import sys
import os
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser, OptionGroup
import pdb

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python plot.py data.pkl FigDispersionRibo.pdf FigDispersionRNA.pdf FigTEchange.pdf' + '\n')

def boxplot_disp(disperRibo, disperRNA, fileOutName):

    idx1 = ~np.isnan(disperRibo)
    idx2 = ~np.isnan(disperRNA)

    disp = [np.log10(disperRibo[idx1]), np.log10(disperRNA[idx2])]

    fig, ax = plt.subplots()
    bp = plt.boxplot(disp, notch=False, sym='b+', vert=True, whis=1.5, widths=0.25)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=1.0)
    ax.set_axisbelow(True)
    ax.set_title('Dispersion comparison', fontsize=15)
    ax.set_xlabel('')
    ax.set_ylabel(r'$log_{10}$(Dispersion)', fontsize=12)
    ax.set_xticklabels(['Ribo-Seq', 'RNA-Seq'])
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=10)

    plt.savefig(fileOutName, format='pdf')

def dotplot_disp(countMean, disperRaw, disperAdj, beta, fileOutName, whatSeq):

    fig, ax = plt.subplots()
    idx = np.logical_and(countMean>1, ~np.isnan(disperAdj))
    ax.scatter(np.log2(countMean[idx]), np.log10(disperRaw[idx]), marker='o', color='silver', s=4, lw=0, label='Raw')
    #IDX = np.logical_and(idx, mthd=='Nelder')
    #ax.scatter(np.log2(countMean[IDX]), np.log10(disperRaw[IDX]), marker='o', color='blue', s=4, lw=0, label='Raw')
    ax.scatter(np.log2(countMean[idx]), np.log10(disperAdj[idx]), marker='o', color='mediumturquoise', s=4, lw=0, label='Adjusted')
    countMeanSorted = np.sort(countMean[idx])
    disperFittedCal = beta[0]/countMeanSorted + beta[1]

    ax.plot(np.log2(countMeanSorted[disperFittedCal > 0]), np.log10(disperFittedCal[disperFittedCal > 0]), color='tomato', label='Fitted')

    ax.legend(loc='upper right', prop={'size':10})

    #xLowBound = 0.0
    #xUpBound  = np.percentile(countMean[idx], 95.0)
    #yLowBound = np.percentile(disperRaw[idx], 2.5)
    #yUpBound  = np.percentile(disperRaw[idx], 97.5)
    #ax.set_xlim(xLowBound, xUpBound)
    #ax.set_ylim(yLowBound, yUpBound)

    if whatSeq:
        ax.set_title('Dispersion estimation for %s' % whatSeq)
    else:
        ax.set_title('Dispersion estimation')
    plt.xlabel('Mean count')
    plt.ylabel('Dispersion')

    plt.savefig(fileOutName, format='pdf')

def dotplot_TEchange(countMean, logFoldChangeTE, disperAdj, padj, threshold, fileOutName):

    index = ~np.isnan(disperAdj)
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
    yLowBound = np.percentile(logFoldChangeTeSig, 2.5)
    yUpBound  = np.percentile(logFoldChangeTeSig, 97.5)

    ax.set_xlim(xLowBound, xUpBound)
    ax.set_ylim(yLowBound, yUpBound)

    plt.savefig(fileOutName, format='pdf')

def make_plots(data, opts):

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        outputNamePrefix = opts.outFile
    else:
        pos = opts.outFile.rfind('.')
        outputNamePrefix = opts.outFile[:pos]

    if opts.dispDiff:
        fileOutNameBoxplot = outputNamePrefix + '.disp.boxplot.pdf'
        boxplot_disp(data.disperAdjRibo, data.disperAdjRNA, fileOutNameBoxplot)

    if opts.dispDiff:
        countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        countRiboMean = np.reshape(countRiboMean, (countRiboMean.size, 1))
        fileOutNameDotplotRibo = outputNamePrefix + '.dispRibo.dotplot.pdf'
        dotplot_disp(countRiboMean, data.disperRawRibo, data.disperAdjRibo, data.betaRibo, fileOutNameDotplotRibo, 'Ribo-Seq')

        countRNAMean = np.mean(data.countRNA / data.libSizesRNA, axis=1)
        countRNAMean = np.reshape(countRNAMean, (countRNAMean.size, 1))
        fileOutNameDotplotRNA = outputNamePrefix + '.dispRNA.dotplot.pdf'
        dotplot_disp(countRNAMean, data.disperRawRNA, data.disperAdjRNA, data.betaRNA, fileOutNameDotplotRNA, 'RNA-Seq')
    else:
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        countMean = np.reshape(countMean, (countMean.size, 1))
        fileOutNameDotplot = outputNamePrefix + '.disp.dotplot.pdf'
        dotplot_disp(countMean, data.disperRaw, data.disperAdj, data.beta, fileOutNameDotplot, '')

    countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
    countMean = np.reshape(countMean, (countMean.size, 1))
    logFoldChangeTE = data.logFoldChangeTE
    if opts.dispDiff:
        disperAdj = data.disperAdjRibo
    else:
        disperAdj = data.disperAdj
    padj = data.padj
    threshold = 0.1
    fileOutNameDotplotTE = outputNamePrefix + '.TEchange.dotplot.pdf'
    dotplot_TEchange(countMean, logFoldChangeTE, disperAdj, padj, threshold, fileOutNameDotplotTE)

def parse_options(argv):

    parser = OptionParser()

    required = OptionGroup(parser, 'REQUIRED')
    required.add_option('-k', action='store', type='string', dest='dataPkl', help='Binary file containing all source data and results of TE change analysis.')
    required.add_option('-o', action='store', type='string', dest='outputPrefix', help='Specify the prefix of the output file name.')

    optional = OptionGroup(parser, 'OPTIONAL')
    optional.add_option('-p', action='store', type='string', dest='plotWhich', default='All', help='Which figure to be plotted. Options: BoxDisp, DotDisp, DotTE or All. [default: All]')

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

    if opts.plotWhich not in ['BoxDisp', 'DotDisp', 'DotTE', 'All']:
        parser.error('-p option can only take \'BoxDisp\', \'DotDisp\', \'DotTE\' or \'All\' as argument.\n')

    return opts

if __name__ == '__main__':

    opts = parse_options(sys.argv)

    with open(opts.dataPkl, 'rb') as FileIn:
        data = pickle.load(FileIn)

    if data.dispDiff:
        if opts.plotWhich in ['BoxDisp', 'All']:
            fileOutNameBoxplot = opts.outputPrefix + '.disp.boxplot.pdf'
            boxplot_disp(data.disperAdjRibo, data.disperAdjRNA, fileOutNameBoxplot)

        if opts.plotWhich in ['DotDisp', 'All']:
            countRiboMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
            countRiboMean = np.reshape(countRiboMean, (countRiboMean.size, 1))
            fileOutNameDotplotRibo = opts.outputPrefix + '.dispRibo.dotplot.pdf'
            dotplot_disp(countRiboMean, data.disperRawRibo, data.disperAdjRibo, data.betaRibo, fileOutNameDotplotRibo, 'Ribo-Seq')

            countRNAMean = np.mean(data.countRNA / data.libSizesRNA, axis=1)
            countRNAMean = np.reshape(countRNAMean, (countRNAMean.size, 1))
            fileOutNameDotplotRNA = opts.outputPrefix + '.dispRNA.dotplot.pdf'
            dotplot_disp(countRNAMean, data.disperRawRNA, data.disperAdjRNA, data.betaRNA, fileOutNameDotplotRNA, 'RNA-Seq')
    else:
        if opts.plotWhich in ['DotDisp', 'All']:
            countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
            countMean = np.reshape(countMean, (countMean.size, 1))
            fileOutNameDotplot = opts.outputPrefix + '.disp.dotplot.pdf'
            dotplot_disp(countMean, data.disperRaw, data.disperAdj, data.beta, fileOutNameDotplot, '')

    if opts.plotWhich in ['DotTE', 'All']:
        countMean = np.mean(data.countRibo / data.libSizesRibo, axis=1)
        countMean = np.reshape(countMean, (countMean.size, 1))
        logFoldChangeTE = data.logFoldChangeTE
        if data.dispDiff:
            disperAdj = data.disperAdjRibo
        else:
            disperAdj = data.disperAdj
        padj = data.padj
        threshold = 0.1
        fileOutNameDotplotTE = opts.outputPrefix + '.TEchange.dotplot.pdf'
        dotplot_TEchange(countMean, logFoldChangeTE, disperAdj, padj, threshold, fileOutNameDotplotTE)

