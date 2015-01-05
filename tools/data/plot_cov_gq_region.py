#!/usr/bin/env python

import sys
import cPickle as pickle
import gzip
import numpy as np
import pysam
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
import pdb

def get_cov(bam, obj, chr):

    bamFile = pysam.AlignmentFile(bam, 'rb')
    covRegion = []
    for j in range(obj.shape[0]):
        left = obj[j, 0]
        right = obj[j, 1]
        cov = []
        pos = []
        covSegment = np.zeros(right-left+1)
        for eachBase in bamFile.pileup(chr, left-1, right, truncate=True):
            cov.append(eachBase.nsegments)
            pos.append(eachBase.reference_pos)

            POS = np.asarray(pos)
            posRegion = np.arange(left-1, right)
            idx = np.in1d(posRegion, POS)
            covSegment[idx] = cov

        covRegion = np.hstack([covRegion, covSegment])

    bamFile.close()

    return covRegion

def get_cov_wrapper(npArray, bamFiles, sizeFactors, chr, strand):

    covRegion = np.zeros(1)
    for k in range(len(bamFiles)):
        covRegion_1Bam = get_cov(bamFiles[k], npArray, chr)
        covRegion_1Bam = covRegion_1Bam / sizeFactors[k]
        covRegion = covRegion + covRegion_1Bam
    covRegion = covRegion / len(bamFiles)
    if strand == '-':
        covRegion = covRegion[::-1]

    return covRegion

def normalize_length(npArray, destLength):

    ratio = npArray.size / float(destLength)
    points = np.arange(0, destLength) * ratio
    idx = np.floor(points).astype(int)
    npArrayNormLen = npArray[idx]

    return npArrayNormLen

def do_plot(COVDIFFdown, COVDIFFup, destLengthGQ, LengthFLK, FileOutName):

    fig, ax = plt.subplots()

    X = np.arange(1, COVDIFFdown.size + 1)
    ax.plot(X, COVDIFFdown, color='orange', linestyle='-', label='TE Down')
    ax.plot(X, COVDIFFup, color='skyblue', linestyle='-', label='TE Up')

    ax.set_xlim(0, COVDIFFdown.size + 1)
    ax.set_ylim(0.6, 1.8)
    ax.legend(loc='upper right', handlelength=3, prop={'size':9})

    ymin, ymax = plt.ylim()
    ax.plot((LengthFLK + 1, LengthFLK + 1), (ymax * 0.02, ymax * 0.98), linestyle='--', color='limegreen', lw=1)
    ax.plot((LengthFLK + destLengthGQ, LengthFLK + destLengthGQ), (ymax * 0.02, ymax * 0.98), linestyle='--', color='limegreen', lw=1)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlabel(r'$GQ\/region$', fontsize=15)
    ax.set_ylabel(r'$Fold Change\/(Treated/Control)$', fontsize=15)
    ax.set_title('Footprint at GQ in 5\'UTR of TE down & up genes')

    plt.savefig(FileOutName, format='pdf')

def get_cov_diff(geneList, gqPklz, destLengthGQ, LengthFLK):

    bamDir = '/cbio/grlab/projects/GuidoWendel/Silvestrol/results/alignments/bam/'

    bamFilesCondA = [bamDir + 'all_Control_1_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam', bamDir + 'all_Control_2_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam']
    bamFilesCondB = [bamDir + 'all_Silvestiol_1_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam', bamDir + 'all_Silvestiol_3_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam']

    sizeFactorsA = [1.134124, 0.802650]
    sizeFactorsB = [1.098537, 0.958124]

    FileIn = gzip.GzipFile(gqPklz, 'rb')
    (GENE_ID, CHR, GENE_STRAND, GQ_ALL_POS, GQ_COV_POS, GQ_INC_POS, FLK5P_COV_POS, FLK3P_COV_POS) = pickle.load(FileIn)
    FileIn.close()
    print 'Load GQ data: done.'

    COVDIFF_ALL = np.zeros(LengthFLK + destLengthGQ + LengthFLK)
    #FREQUNENCY  = np.zeros(LengthFLK + destLengthGQ + LengthFLK)

    FileIn = open(geneList, 'r')

    const = 1.0
    cutoff = 4.0
    cnt = 1
    for line in FileIn:
        if line.strip().startswith('ENSG'):
            geneID = line.strip()

            if geneID not in GQ_COV_POS:
                print '\n%s dose not exists in annotation file.' % geneID
                continue

            IDX = np.nonzero(np.asarray(GQ_COV_POS[geneID]) == 1)[0]
            for idx in IDX:
                expGQ = GQ_ALL_POS[geneID][idx]
                expGQcov_A = get_cov_wrapper(expGQ, bamFilesCondA, sizeFactorsA, CHR[geneID], GENE_STRAND[geneID])
                expGQcov_B = get_cov_wrapper(expGQ, bamFilesCondB, sizeFactorsB, CHR[geneID], GENE_STRAND[geneID])
                expGQcovDiff = (expGQcov_B + const) / (expGQcov_A + const)
                idxMask = np.logical_and(expGQcov_A < cutoff, expGQcov_B < cutoff)
                expGQcovDiff[idxMask] = 0.0

                FLK5P = FLK5P_COV_POS[geneID][idx]
                if type(FLK5P).__module__ == np.__name__:
                    FLK5Pcov_A = get_cov_wrapper(FLK5P, bamFilesCondA, sizeFactorsA, CHR[geneID], GENE_STRAND[geneID])
                    FLK5Pcov_B = get_cov_wrapper(FLK5P, bamFilesCondB, sizeFactorsB, CHR[geneID], GENE_STRAND[geneID])
                    FLK5PcovDiff = (FLK5Pcov_B + const) / (FLK5Pcov_A + const)
                    idxMask = np.logical_and(FLK5Pcov_A < cutoff, FLK5Pcov_B < cutoff)
                    FLK5PcovDiff[idxMask] = 0.0
                else:
                    FLK5PcovDiff = np.zeros(0)

                FLK3P = FLK3P_COV_POS[geneID][idx]
                if type(FLK3P).__module__ == np.__name__:
                    FLK3Pcov_A = get_cov_wrapper(FLK3P, bamFilesCondA, sizeFactorsA, CHR[geneID], GENE_STRAND[geneID])
                    FLK3Pcov_B = get_cov_wrapper(FLK3P, bamFilesCondB, sizeFactorsB, CHR[geneID], GENE_STRAND[geneID])
                    FLK3PcovDiff = (FLK3Pcov_B + const) / (FLK3Pcov_A + const)
                    idxMask = np.logical_and(FLK3Pcov_A < cutoff, FLK3Pcov_B < cutoff)
                    FLK3PcovDiff[idxMask] = 0.0
                else:
                    FLK3PcovDiff = np.zeros(0)

                expGQcovDiffNormLen = normalize_length(expGQcovDiff, destLengthGQ)

                COVDIFF_GQ_ONE = expGQcovDiffNormLen.copy()
                COVDIFF_5P_ONE = np.hstack([np.zeros(LengthFLK - FLK5PcovDiff.size), FLK5PcovDiff])
                COVDIFF_3P_ONE = np.hstack([FLK3PcovDiff, np.zeros(LengthFLK - FLK3PcovDiff.size)])
                COVDIFF_ONE = np.hstack([COVDIFF_5P_ONE, COVDIFF_GQ_ONE, COVDIFF_3P_ONE])
                #COVDIFF_ALL = COVDIFF_ALL + COVDIFF_ONE
                COVDIFF_ALL = np.vstack([COVDIFF_ALL, COVDIFF_ONE])
                #idxEff = np.nonzero(COVDIFF_ONE != 0.0)
                #FREQUNENCY[idxEff] = FREQUNENCY[idxEff] + 1.0

            sys.stdout.flush()
            print '\r%i genes finished ...' % cnt ,
            cnt += 1

    #COVDIFF = COVDIFF_ALL / FREQUNENCY
    COVDIFF = gmean(COVDIFF_ALL, axis=0).data

    FileIn.close()

    return COVDIFF

if __name__ == '__main__':

    outputTxt = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTE_DownUp.prot.5p.utr5.all.cutoff4.txt'
    outputPdf = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTE_DownUp.prot.5p.utr5.all.cutoff4.pdf'
    destLengthGQ = 30.0
    LengthFLK = 40.0

    #geneList = '/cbio/grlab/projects/RibosomeFootprint/Lists/SilTEdown.prot.5p.txt'
    #gqPklz = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTEdown.prot.5p.utr5.all.cutoff3.gq'
    #COVDIFF_Down = get_cov_diff(geneList, gqPklz, destLengthGQ, LengthFLK)

    #geneList = '/cbio/grlab/projects/RibosomeFootprint/Lists/SilTEup.prot.5p.txt'
    #gqPklz = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTEup.prot.5p.utr5.all.cutoff3.gq'
    #COVDIFF_Up = get_cov_diff(geneList, gqPklz, destLengthGQ, LengthFLK)

    #column1 = np.arange(1, COVDIFF_Down.size+1).reshape((COVDIFF_Down.size, 1)).astype(int).astype(str)
    #column2 = COVDIFF_Down.reshape((COVDIFF_Down.size, 1)).astype(str)
    #column3 = COVDIFF_Up.reshape((COVDIFF_Up.size, 1)).astype(str)
    #outNdarray = np.hstack([column1, column2, column3])
    #np.savetxt(outputTxt, outNdarray, fmt='%s', delimiter='\t', comments='')
    #print 'Save file: done.'

    COVDIFF_Down = np.loadtxt(outputTxt, dtype=float, skiprows=0, delimiter='\t', usecols=(1,))
    COVDIFF_Up   = np.loadtxt(outputTxt, dtype=float, skiprows=0, delimiter='\t', usecols=(2,))
    do_plot(COVDIFF_Down, COVDIFF_Up, destLengthGQ, LengthFLK, outputPdf)
    print 'Plot: done.'

