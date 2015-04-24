#!/usr/bin/env python
"""
Command line argument parser.
"""

import os
import sys
from optparse import OptionParser, OptionGroup

def parse_options(argv):

    parser = OptionParser(usage='usage: %prog [options] arguments')

    required = OptionGroup(parser, 'REQUIRED')

    required.add_option('-e', action='store', type='string', dest='exptOutline', help='Text file describing experiment Outline. Must follow required format, please see the manual.')
    required.add_option('-c', action='store', type='string', dest='cntFile', help='Text file containing the count data. Header line must be consistent with information in experiment Outline.')
    required.add_option('-o', action='store', type='string', dest='outFile', help='Tab delimited text file containing the results.')

    optional = OptionGroup(parser, 'OPTIONAL')

    optional.add_option('-d', action='store', type='int', dest='dispDiff', default=0, help='Allow different dispersions for Ribo-seq and RNA-seq count data. Off: 0; On: 1. [default: 0]')
    optional.add_option('-s', action='store', type='int', dest='sumCntCutoff', default=10, help='Set the sum of normalized read count as the threshold to do the test. This option applies for both Ribo-seq and RNA-seq data. [default: 10]')
    optional.add_option('-i', action='store', type='float', dest='dispInitial', default=0.01, help='Set the initial dispersion to start the estimation. [default: 0.01]')
    optional.add_option('-m', action='store', type='string', dest='multiTest', default='BH', help='Method for multiple test correction. Options: BH (Benjamini-Hochberg); Bonferroni. [default: BH]')
    optional.add_option('-r', action='store', type='int', dest='rankResult', default=0, help='Rank the result table in ascending order by a specific column. Adjusted p value: 1; TE change: 2; Gene id: 3; Keep the order as in count file: 0. [default: 0]')
    optional.add_option('-p', action='store', type='int', dest='plots', default=0, help='Make plots to show the data and results. Plots are in pdf format. On: 1; Off: 0. [default: 0]')
    optional.add_option('-q', action='store', type='float', dest='cutoffFDR', default=0.1, help='Set the FDR cutoff for significant case to plot. [default: 0.1]')

    parser.add_option_group(required)
    parser.add_option_group(optional)

    (opts, args) = parser.parse_args(argv)

    if len(argv) < 2:
        parser.print_help()
        sys.exit()

    mandatories = ['exptOutline', 'cntFile', 'outFile']
    for eachOpt in mandatories:
        if not opts.__dict__[eachOpt]:
            parser.error('-%s is a required option.\n' % eachOpt[0])

        if eachOpt in mandatories[:2] and not os.path.exists(opts.__dict__[eachOpt]):
            sys.stderr.write('\nError: File \'%s\' does not exist.\n\n' % opts.__dict__[eachOpt])
            sys.exit()

        if eachOpt == mandatories[2]:
            if not os.path.dirname(opts.__dict__[eachOpt]):
                opts.__dict__[eachOpt] = os.getcwd() + os.getcwd()[0] + opts.__dict__[eachOpt]
            if not os.path.exists(os.path.dirname(opts.__dict__[eachOpt])):
                try:
                    os.makedirs(os.path.dirname(opts.__dict__[eachOpt]))
                except OSError:
                    sys.stderr.write('\nError: Failed to create directory: \'%s\' \n\n' % os.path.dirname(opts.__dict__[eachOpt]))
                    sys.exit()

            opts.__dict__['resPath'] = os.path.dirname(opts.__dict__[eachOpt]) + os.getcwd()[0]

    if opts.dispDiff not in [0, 1]:
        parser.error('-d option can only take either 0 or 1 as argument.\n')

    if opts.sumCntCutoff < 1.0:
        parser.error('-s option cannot take values smaller than 1 as argument.\n')

    if opts.dispInitial <= 0.0:
        parser.error('-i option can only take positive value as argument.\n')

    if opts.multiTest not in ['BH', 'Bonferroni', 'Holm', 'Hochberg', 'Hommel', 'BY', 'TSBH']:
        parser.error('-m option can only take \'BH\' or \'Bonferroni\' as argument.\n')

    if opts.rankResult not in [0, 1, 2, 3]:
        parser.error('-r option can only take 0, 1, 2 or 3 as argument.\n')

    if opts.plots not in [0, 1]:
        parser.error('-p option can only take either 0 or 1 as argument.\n')

    if opts.cutoffFDR < 0.0 or opts.cutoffFDR > 1.0:
        parser.error('-q option can oly take values between zero and one.\n')

    return opts

