#!/usr/bin/env python

import sys
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import pdb
 
def cal_tpr_fpr(simFileName, data):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    for fdr in np.arange(0.05, 0.151, 0.025):
        with np.errstate(invalid='ignore'):
            detectPosID = data.geneIDs[data.padj <= fdr]
        tpr = np.in1d(detectPosID, truePosID).nonzero()[0].size / float(truePosID.size)
        fpr = np.in1d(detectPosID, truePosID, invert=True).nonzero()[0].size / float(trueNegID.size)
        print 'FDR < %1.3f, True Positive Rate = %1.3f; False Positive Rate = %1.3f' % (fdr, tpr, fpr)

    print '*'*20

def cal_tpr_fpr_babel(simFileName, babelFileName):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    babelFDRstr = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(6,))
    babelFDRstr[babelFDRstr=='NA'] = 'nan'
    babelFDR = babelFDRstr.astype(float)
    babelID = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(1,))

    for fdr in np.arange(0.05, 0.151, 0.025):
        with np.errstate(invalid='ignore'):
            detectPosID = babelID[babelFDR <= fdr]
        tpr = np.in1d(detectPosID, truePosID).nonzero()[0].size / float(truePosID.size)
        fpr = np.in1d(detectPosID, truePosID, invert=True).nonzero()[0].size / float(trueNegID.size)
        print 'FDR < %1.3f, True Positive Rate = %1.3f; False Positive Rate = %1.3f' % (fdr, tpr, fpr)

    print '*'*20

def cal_TPR_FPR(simFileName, data):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    idx = np.argsort(data.pval, axis=None)
    geneIDs = data.geneIDs[idx]
    TPR = []
    FPR = []
    for i in range(geneIDs.size):
        tpr = np.in1d(geneIDs[0:i+1], truePosID).nonzero()[0].size / float(truePosID.size)
        fpr = np.in1d(geneIDs[0:i+1], truePosID, invert=True).nonzero()[0].size / float(trueNegID.size)
        TPR.extend([tpr])
        FPR.extend([fpr])

    print 'calculate TPR, FPR: Done!'

    return (TPR, FPR)

def cal_TPR_FPR_babel(simFileName, babelFileName):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    babelPvalstr = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(5,))
    babelPvalstr[babelPvalstr=='NA'] = 'nan'
    babelPval = babelPvalstr.astype(float)
    babelID = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(1,))

    idx = np.argsort(babelPval, axis=None)
    geneIDs = babelID[idx]
    TPR = []
    FPR = []
    for i in range(geneIDs.size):
        tpr = np.in1d(geneIDs[0:i+1], truePosID).nonzero()[0].size / float(truePosID.size)
        fpr = np.in1d(geneIDs[0:i+1], truePosID, invert=True).nonzero()[0].size / float(trueNegID.size)
        TPR.extend([tpr])
        FPR.extend([fpr])

    print 'calculate TPR, FPR: Done!'

    return (TPR, FPR)

def plot_roc(TPR1, FPR1, TPR2, FPR2, TPR3, FPR3, rocFileName):

    fig, ax = plt.subplots()

    ax.plot(FPR1, TPR1, linestyle='-', color='tomato', label='Single dispersion')
    ax.plot(FPR2, TPR2, linestyle='-', color='blue', label='Two dispersions')
    ax.plot(FPR3, TPR3, linestyle='-', color='violet', label='babel')

    ax.legend(loc='lower right', prop={'size':11})

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    ax.set_title(r'ROC curve')
    ax.set_xlabel(r'$False\/positive\/rate$', fontsize=15)
    ax.set_ylabel(r'$True\/positive\/rate$', fontsize=15)

    plt.savefig(rocFileName, format='pdf')

    print 'plot: Done!'

if __name__ == '__main__':

    simFileName = sys.argv[1]

    with open(sys.argv[2], 'rb') as FileIn:
        data1 = pickle.load(FileIn)

    with open(sys.argv[3], 'rb') as FileIn:
        data2 = pickle.load(FileIn)

    babelFileName = sys.argv[4]

    rocFileName = sys.argv[-1]

    cal_tpr_fpr(simFileName, data1)
    cal_tpr_fpr(simFileName, data2)
    cal_tpr_fpr_babel(simFileName, babelFileName)
    TPR1, FPR1 = cal_TPR_FPR(simFileName, data1)
    TPR2, FPR2 = cal_TPR_FPR(simFileName, data2)
    TPR3, FPR3 = cal_TPR_FPR_babel(simFileName, babelFileName)
    plot_roc(TPR1, FPR1, TPR2, FPR2, TPR3, FPR3, rocFileName)
