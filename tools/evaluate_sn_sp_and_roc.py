#!/usr/bin/env python

import sys
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import pdb

def cal_SN_SP(simFileName, data):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    fdr = 0.15
    with np.errstate(invalid='ignore'):
        observedPosID = data.geneIDs[data.padj <= fdr]
        observedNegID = data.geneIDs[data.padj >  fdr]
    sn = np.intersect1d(observedPosID, truePosID).size / float(truePosID.size)
    sp = np.intersect1d(observedNegID, trueNegID).size / float(trueNegID.size)
    print 'FDR < %1.2f, Sensitivity = %1.3f, Specificity = %1.3f' % (fdr, sn, sp)

    print '*'*20

    idx = np.argsort(data.pval, axis=None)
    geneIDs = data.geneIDs[idx]
    SN = np.zeros(geneIDs.size)
    SP = np.zeros(geneIDs.size)
    P  = np.zeros(geneIDs.size)
    for i in range(geneIDs.size):
        sn = np.intersect1d(geneIDs[:i+1], truePosID).size / float(truePosID.size)
        sp = np.intersect1d(geneIDs[i+1:], trueNegID).size / float(trueNegID.size)
        SN[i] = sn
        SP[i] = sp
        P[i]  = data.pval[idx[i]][0]

    return (SN, SP, P)

def plot_sn_sp(SN1, SP1, P1, SN2, SP2, P2, SN3, SP3, P3, SN4, SP4, P4, SN5, SP5, P5, snspFigName):

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))

    ax1.plot(P1[~np.isnan(P1)], SN1[~np.isnan(P1)], linestyle='-', color='orange', label=r'Ribo, Gamma $\alpha$=1.5')
    ax1.plot(P2[~np.isnan(P2)], SN2[~np.isnan(P2)], linestyle='-', color='skyblue', label=r'RNA, Gamma $\alpha$=1.5')
    ax1.plot(P3[~np.isnan(P3)], SN3[~np.isnan(P3)], linestyle='--', color='orange', label=r'Ribo, Gamma $\alpha$=0.8')
    ax1.plot(P4[~np.isnan(P4)], SN4[~np.isnan(P4)], linestyle='--', color='skyblue', label=r'RNA, Gamma $\alpha$=0.6')
    ax1.plot(P5[~np.isnan(P5)], SN5[~np.isnan(P5)], linestyle='-', color='tomato', label=r'Ribo & RNA, combined')

    ax1.legend(loc='lower right', handlelength=3, prop={'size':7})

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=8)

    ax1.set_title(r'Sensitivity', fontsize=10)
    ax1.set_xlabel(r'$P\/value$', fontsize=10)
    ax1.set_ylabel(r'$Sensitivity$', fontsize=10)

    ax2.plot(P1[~np.isnan(P1)], SP1[~np.isnan(P1)], linestyle='-', color='orange', label=r'Ribo, Gamma $\alpha$=1.5')
    ax2.plot(P2[~np.isnan(P2)], SP2[~np.isnan(P2)], linestyle='-', color='skyblue', label=r'RNA, Gamma $\alpha$=1.5')
    ax2.plot(P3[~np.isnan(P3)], SP3[~np.isnan(P3)], linestyle='--', color='orange', label=r'Ribo, Gamma $\alpha$=0.8')
    ax2.plot(P4[~np.isnan(P4)], SP4[~np.isnan(P4)], linestyle='--', color='skyblue', label=r'RNA, Gamma $\alpha$=0.6')
    ax2.plot(P5[~np.isnan(P5)], SP5[~np.isnan(P5)], linestyle='-', color='tomato', label=r'Ribo & RNA, combined')

    #ax2.legend(loc='lower left', handlelength=3, prop={'size':7})

    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    ax2.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='y', labelsize=8)

    ax2.set_title(r'Specificity', fontsize=10)
    ax2.set_xlabel(r'$P\/value$', fontsize=10)
    #ax2.set_xlabel(r'$\log_{10}(P\/value)$', fontsize=7)
    ax2.set_ylabel(r'$Specificity$', fontsize=10)

    plt.savefig(snspFigName, format='pdf', bbox_inches='tight') 

def cal_TPR_FPR(simFileName, data):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    fdr = 0.15
    with np.errstate(invalid='ignore'):
        observedPosID = data.geneIDs[data.padj <= fdr]
    tpr = np.intersect1d(observedPosID, truePosID).size / float(truePosID.size)
    fpr = np.setdiff1d(observedPosID, truePosID).size / float(trueNegID.size)
    print 'FDR < %1.2f, True positive rate = %1.3f, False positive rate = %1.3f' % (fdr, tpr, fpr)

    print '*'*20

    idx = np.argsort(data.pval, axis=None)
    geneIDs = data.geneIDs[idx]
    TPR = np.zeros(geneIDs.size)
    FPR = np.zeros(geneIDs.size)
    for i in range(geneIDs.size):
        tpr = np.intersect1d(geneIDs[:i+1], truePosID).size / float(truePosID.size)
        fpr = np.setdiff1d(geneIDs[:i+1], truePosID).size / float(trueNegID.size)
        TPR[i] = tpr
        FPR[i] = fpr

    return (TPR, FPR)

def cal_TPR_FPR_babel(simFileName, babelFileName):

    diffMarker = np.loadtxt(simFileName, dtype=int, delimiter='\t', skiprows=1, usecols=(-1,))
    entryID = np.loadtxt(simFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    truePosID = entryID[diffMarker > 0]
    trueNegID = entryID[diffMarker == 0]

    babelPvalstr = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(5,))
    babelPvalstr[babelPvalstr=='NA'] = 'nan'
    babelPval = babelPvalstr.astype(float)
    babelPadjstr = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(6,))
    babelPadjstr[babelPadjstr=='NA'] = 'nan'
    babelPadj = babelPadjstr.astype(float)
    babelID = np.loadtxt(babelFileName, dtype=str, delimiter='\t', skiprows=1, usecols=(1,))

    fdr = 0.15
    with np.errstate(invalid='ignore'):
        observedPosID = babelID[babelPadj <= fdr]
    tpr = np.intersect1d(observedPosID, truePosID).size / float(truePosID.size)
    fpr = np.setdiff1d(observedPosID, truePosID).size / float(trueNegID.size)
    print 'FDR < %1.2f, True positive rate = %1.3f, False positive rate = %1.3f' % (fdr, tpr, fpr)

    print '*'*20

    idx = np.argsort(babelPval, axis=None)
    geneIDs = babelID[idx]
    TPR = np.zeros(geneIDs.size)
    FPR = np.zeros(geneIDs.size)
    for i in range(geneIDs.size):
        tpr = np.intersect1d(geneIDs[:i+1], truePosID).size / float(truePosID.size)
        fpr = np.setdiff1d(geneIDs[:i+1], truePosID).size / float(trueNegID.size)
        TPR[i] = tpr
        FPR[i] = fpr

    return (TPR, FPR)

def plot_roc(TPR1, FPR1, TPR2, FPR2, TPR3, FPR3, rocFigName):

    fig, ax = plt.subplots()

    ax.plot(FPR1, TPR1, linestyle='-', color='tomato', label='Single dispersion')
    ax.plot(FPR2, TPR2, linestyle='-', color='skyblue', label='Two dispersions')
    ax.plot(FPR3, TPR3, linestyle='-', color='orange', label='babel')

    ax.legend(loc='lower right', prop={'size':11})

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    ax.set_title(r'ROC curve')
    ax.set_xlabel(r'$False\/positive\/rate$', fontsize=15)
    ax.set_ylabel(r'$True\/positive\/rate$', fontsize=15)

    plt.savefig(rocFigName, format='pdf', bbox_inches='tight')

if __name__ == '__main__':

    simFileName1 = '../exp/Sim/Sim.Ribo.Rep3.G2K.Diff1K.Sh1.5.Sc0.5.cnt.txt'
    simFileName2 = '../exp/Sim/Sim.Rna.Rep3.G2K.Diff1K.Sh1.5.Sc0.5.cnt.txt'
    simFileName3 = '../exp/Sim/Sim.Ribo.Rep3.G2K.Diff1K.Sh0.8.Sc0.5.cnt.txt'
    simFileName4 = '../exp/Sim/Sim.Rna.Rep3.G2K.Diff1K.Sh0.6.Sc0.5.cnt.txt'
    simFileName5 = '../exp/Sim/Sim.Ribo.Rep3.G2K.Diff1K.Sh0.8.Sc0.5.cnt.txt'

    resTEtest1 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRb.Sh1.5.Sc0.5.res0.pkl'
    resTEtest2 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRna.Sh1.5.Sc0.5.res0.pkl'
    resTEtest3 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRb.Sh0.8.Sc0.5.res0.pkl'
    resTEtest4 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRna.Sh0.6.Sc0.5.res0.pkl'
    resTEtest5 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.res0.pkl'

    snspFigName = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1K.SnSp.res0.pdf'

    with open(resTEtest1, 'rb') as FileIn:
        data1 = pickle.load(FileIn)
    with open(resTEtest2, 'rb') as FileIn:
        data2 = pickle.load(FileIn)
    with open(resTEtest3, 'rb') as FileIn:
        data3 = pickle.load(FileIn)
    with open(resTEtest4, 'rb') as FileIn:
        data4 = pickle.load(FileIn)
    with open(resTEtest5, 'rb') as FileIn:
        data5 = pickle.load(FileIn)

    print 'Start SN1, SP1'
    SN1, SP1, P1 = cal_SN_SP(simFileName1, data1)
    print 'Start SN2, SP2'
    SN2, SP2, P2 = cal_SN_SP(simFileName2, data2)
    print 'Start SN3, SP3'
    SN3, SP3, P3 = cal_SN_SP(simFileName3, data3)
    print 'Start SN4, SP4'
    SN4, SP4, P4 = cal_SN_SP(simFileName4, data4)
    print 'Start SN5, SP5'
    SN5, SP5, P5 = cal_SN_SP(simFileName5, data5)
    plot_sn_sp(SN1, SP1, P1, SN2, SP2, P2, SN3, SP3, P3, SN4, SP4, P4, SN5, SP5, P5, snspFigName)
    print 'plot SN SP: Done!\n'

    babelFileName = '../exp/Sim/Sim.RnaRibo.Rep3.G2K.Diff1K.Sh0.8Rb.Sh0.6Rna.Sc0.5.res.babel.txt'
    rocFigName = '../exp/Sim/Sim.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.roc.pdf'
    resTEtest6 = '../exp/Sim/Sim.Merged.Rep3.G2K.Diff1KRbRna.Sh0.8Rb.Sh0.6Rna.Sc0.5.res1.pkl'
    with open(resTEtest6, 'rb') as FileIn:
        data6 = pickle.load(FileIn)

    print 'Start TPR1, FPR1'
    TPR1, FPR1 = cal_TPR_FPR(simFileName5, data5)
    print 'Start TPR2, FPR2'
    TPR2, FPR2 = cal_TPR_FPR(simFileName5, data6)
    print 'Start TPR3, FPR3'
    TPR3, FPR3 = cal_TPR_FPR_babel(simFileName5, babelFileName)
    plot_roc(TPR1, FPR1, TPR2, FPR2, TPR3, FPR3, rocFigName)
    print 'plot ROC: Done!'
