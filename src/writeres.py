import numpy as np
import cPickle as pickle

def write_result(data, opts):

    geneIDs = data.geneIDs
    pval = data.pval.astype(str)
    padj = data.padj.astype(str)
    TEctl = data.TEctl.astype(str)
    TEtrt = data.TEtrt.astype(str)
    logFoldChangeTE = data.logFoldChangeTE.astype(str)
    countRibo = np.around(data.countRibo / data.libSizesRibo).astype(str)
    countRna  = np.around(data.countRna  / data.libSizesRna ).astype(str)
    sampleNameRibo = '\t'.join(data.headerRibo.tolist())
    sampleNameRna  = '\t'.join(data.headerRna.tolist())

    if opts.dispDiff:
        dispRawRibo = data.dispRawRibo.astype(str)
        dispRawRna  = data.dispRawRna.astype(str)
        dispRawMthd = data.dispRawMthd.astype(str)
        dispRawConv = data.dispRawConv.astype(str)
        dispRawConv[dispRawConv=='1.0'] = 'Y'
        dispRawConv[dispRawConv=='0.0'] = 'N'
        dispFittedRibo = data.dispFittedRibo.astype(str)
        dispFittedRna  = data.dispFittedRna.astype(str)
        dispAdjRibo = data.dispAdjRibo.astype(str)
        dispAdjRna  = data.dispAdjRna.astype(str)
        dispAdjMthd = data.dispAdjMthd.astype(str)
        dispAdjConv = data.dispAdjConv.astype(str)
        dispAdjConv[dispAdjConv=='1.0'] = 'Y'
        dispAdjConv[dispAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, dispRawRibo, dispRawRna, dispRawMthd, dispRawConv, dispFittedRibo, dispFittedRna, dispAdjRibo, dispAdjRna, dispAdjMthd, dispAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE, countRibo, countRna])
        header = 'geneIDs\tdisperRawRibo\tdisperRawRNA\tdisperRawMthd\tdisperRawConverge\tdisperFittedRibo\tdisperFittedRNA\tdisperAdjRibo\tdisperAdjRNA\tdisperAdjMthd\tdisperAdjConverge\tpval\tpadj\tTEctl\tTEtrt\tlog2FoldChangeTE\t' + sampleNameRibo + '\t' + sampleNameRna
    else:
        dispRaw = data.dispRaw.astype(str)
        dispRawConv = data.dispRawConv.astype(str)
        dispRawConv[dispRawConv=='1.0'] = 'Y'
        dispRawConv[dispRawConv=='0.0'] = 'N'
        dispFitted = data.dispFitted.astype(str)
        dispAdj = data.dispAdj.astype(str)
        dispAdjConv = data.dispAdjConv.astype(str)
        dispAdjConv[dispAdjConv=='1.0'] = 'Y'
        dispAdjConv[dispAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, dispRaw, dispRawConv, dispFitted, dispAdj, dispAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE, countRibo, countRna])
        header = 'geneIDs\tdisperRaw\tdisperRawConverge\tdisperFitted\tdisperAdj\tdisperAdjConverge\tpval\tpadj\tTEctl\tTEtrt\tlog2FoldChangeTE\t' + sampleNameRibo + '\t' + sampleNameRna

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

