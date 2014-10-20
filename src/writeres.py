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
    countRNA  = np.around(data.countRNA  / data.libSizesRNA ).astype(str)
    sampleNameRibo = '\t'.join(data.headerRibo.tolist())
    sampleNameRNA  = '\t'.join(data.headerRNA.tolist())

    if opts.dispDiff:
        disperRawRibo = data.disperRawRibo.astype(str)
        disperRawRNA  = data.disperRawRNA.astype(str)
        disperRawMthd = data.disperRawMthd.astype(str)
        disperRawConv = data.disperRawConv.astype(str)
        disperRawConv[disperRawConv=='1.0'] = 'Y'
        disperRawConv[disperRawConv=='0.0'] = 'N'
        disperFittedRibo = data.disperFittedRibo.astype(str)
        disperFittedRNA  = data.disperFittedRNA.astype(str)
        disperAdjRibo = data.disperAdjRibo.astype(str)
        disperAdjRNA  = data.disperAdjRNA.astype(str)
        disperAdjMthd = data.disperAdjMthd.astype(str)
        disperAdjConv = data.disperAdjConv.astype(str)
        disperAdjConv[disperAdjConv=='1.0'] = 'Y'
        disperAdjConv[disperAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, disperRawRibo, disperRawRNA, disperRawMthd, disperRawConv, disperFittedRibo, disperFittedRNA, disperAdjRibo, disperAdjRNA, disperAdjMthd, disperAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE, countRibo, countRNA])
        header = 'geneIDs\tdisperRawRibo\tdisperRawRNA\tdisperRawMthd\tdisperRawConverge\tdisperFittedRibo\tdisperFittedRNA\tdisperAdjRibo\tdisperAdjRNA\tdisperAdjMthd\tdisperAdjConverge\tpval\tpadj\tTEctl\tTEtrt\tlog2FoldChangeTE\t' + sampleNameRibo + '\t' + sampleNameRNA
    else:
        disperRaw = data.disperRaw.astype(str)
        disperRawConv = data.disperRawConv.astype(str)
        disperRawConv[disperRawConv=='1.0'] = 'Y'
        disperRawConv[disperRawConv=='0.0'] = 'N'
        disperFitted = data.disperFitted.astype(str)
        disperAdj = data.disperAdj.astype(str)
        disperAdjConv = data.disperAdjConv.astype(str)
        disperAdjConv[disperAdjConv=='1.0'] = 'Y'
        disperAdjConv[disperAdjConv=='0.0'] = 'N'
        outNdarrayUnsorted = np.hstack([geneIDs, disperRaw, disperRawConv, disperFitted, disperAdj, disperAdjConv, pval, padj, TEctl, TEtrt, logFoldChangeTE, countRibo, countRNA])
        header = 'geneIDs\tdisperRaw\tdisperRawConverge\tdisperFitted\tdisperAdj\tdisperAdjConverge\tpval\tpadj\tTEctl\tTEtrt\tlog2FoldChangeTE\t' + sampleNameRibo + '\t' + sampleNameRNA

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

