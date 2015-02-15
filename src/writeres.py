import numpy as np
import cPickle as pickle

def write_result(data, opts):

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
        dispAdjConv = data.dispAdjConv.astype(str)
        dispAdjConv[dispAdjConv=='1.0'] = 'Y'
        dispAdjConv[dispAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdjRibo, dispAdjRna, dispAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE])
        header = 'geneIDs\tdisperRibo\tdisperRNA\tdisperConv\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\t' % (nameCondA, nameCondB, nameCondB, nameCondA)
    else:
        dispAdj = data.dispAdj.astype(str)
        dispAdjConv = data.dispAdjConv.astype(str)
        dispAdjConv[dispAdjConv=='1.0'] = 'Y'
        dispAdjConv[dispAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdj, dispAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE])
        header = 'geneIDs\tdisper\tdisperConv\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\t' % (nameCondA, nameCondB, nameCondB, nameCondA)

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

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        pklFile = opts.outFile + '.pkl'
    else:
        pos = opts.outFile.rfind('.')
        outputNamePrefix = opts.outFile[:pos]
        pklFile = outputNamePrefix + '.pkl'
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

