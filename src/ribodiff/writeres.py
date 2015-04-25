#!/usr/bin/env python 
"""
Output result files.
"""

import numpy as np
import cPickle as pickle

def write_result(data, opts):
    """
    Output the result to plain text format.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input argument to the main TE function 
    @type opts: Instance
    """

    geneIDs = data.geneIDs
    pval = data.pval.astype(str)
    padj = data.padj.astype(str)
    TEctl = data.TEctl.astype(str)
    TEtrt = data.TEtrt.astype(str)
    nameCondA = data.nameCondA
    nameCondB = data.nameCondB
    logFoldChangeTE = data.logFoldChangeTE.astype(str)

    if opts.dispDiff:
        dispAdjRibo = data.dispAdjRibo.astype(str)
        dispAdjRna  = data.dispAdjRna.astype(str)
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdjRibo, dispAdjRna, pval, padj, TEctl, TEtrt, logFoldChangeTE])
        header = 'geneIDs\tdisperRibo\tdisperRNA\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\t' % (nameCondA, nameCondB, nameCondB, nameCondA)
    else:
        dispAdj = data.dispAdj.astype(str)
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdj, pval, padj, TEctl, TEtrt, logFoldChangeTE])
        header = 'geneIDs\tdisper\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\t' % (nameCondA, nameCondB, nameCondB, nameCondA)

    if opts.rankResult == 0:
        outNdarray = outNdarrayUnsorted.copy()
    else:
        if opts.rankResult == 1:
            idx = np.argsort(data.padj, axis=None)
        elif opts.rankResult == 2:
            idx = np.argsort(data.logFoldChangeTE, axis=None)
        elif opts.rankResult == 3:
            idx = np.argsort(data.geneIDs, axis=None)
        else:
            pass

        outNdarray = outNdarrayUnsorted[idx]

    np.savetxt(opts.outFile, outNdarray, fmt='%s', delimiter='\t', header=header, comments='')

def save_data(data, opts):
    """
    Saving in python object form.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input argument to the main TE function 
    @type opts: Instance
    """

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        pklFile = opts.outFile + '.pkl'
    else:
        pos = opts.outFile.rfind('.')
        if opts.outFile[pos:pos+2] == './':
            outputNamePrefix = opts.outFile
        else:
            outputNamePrefix = opts.outFile[:pos]
        pklFile = outputNamePrefix + '.pkl'
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

