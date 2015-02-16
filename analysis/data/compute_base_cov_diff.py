#!/usr/bin/env python

import sys
import cPickle as pickle
import gzip
import string
import pysam
import numpy as np
import pdb

def read_genome(genomeFile):

    FileIn = open(genomeFile, 'r')

    genomeSeq = {}
    id = ''
    seq = ''
    n = 0
    for line in FileIn:
        sys.stdout.flush()
        if line.lstrip().startswith('>'):
            n += 1
            print '\rStart reading chromosome: %i...' % n,
            if id and seq:
                genomeSeq[id] = seq
            id = line.strip()[1:]
            seq = ''
        else:
            seq = seq + line.strip()

    FileIn.close()

    print '\nRead genome sequence done.'

    return genomeSeq

def find_len_extreme(Obj, seqFeature):

    """
    If there are two or more longest/shortest sequences, this function only returns the first one.
    The "shortest sequences" means the shortest one that longer than 0 bp. If all the isoforms are 
    0 bp, it returns idx as [0].
    """

    objLength = []
    extremeIdx = []
    for i in range(len(Obj)):
        if Obj[i].shape[0] == 0:
            eachLength = 0
        else:
            eachLength = sum(Obj[i][:,1] - Obj[i][:,0] + 1)
        objLength.extend([eachLength])

    if seqFeature == 'LONGEST':
        extremeIdx.append(objLength.index(max(objLength)))
    elif seqFeature == 'SHORTEST':
        if min(objLength) != 0:
            extremeIdx.append(objLength.index(min(objLength)))
        else:
            if sum(objLength) > 0:
                extremeIdx.append(objLength.index(sorted(each for each in objLength if each > 0)[0]))
            else:
                extremeIdx.append(objLength.index(0))
    else:
        pass

    return extremeIdx

def get_seq(genome, obj, Chr, Strand):

    seq = ''
    for j in range(obj.shape[0]):
        left = obj[j, 0]
        right = obj[j, 1]
        seqExon = genome[Chr][left-1:right]
        seq = seq + seqExon

    if Strand == '-':
        seq = seq[::-1].translate(string.maketrans('ATGCNatgcn', 'TACGNtacgn'))

    return seq

def get_cov(bam, obj, Chr):

    bamFile = pysam.AlignmentFile(bam, 'rb')
    covRegion = []
    for j in range(obj.shape[0]):
        left = obj[j, 0]
        right = obj[j, 1]
        cov = []
        pos = []
        covSegment = np.zeros(right-left+1)
        for eachBase in bamFile.pileup(Chr, left-1, right, truncate=True):
            cov.append(eachBase.nsegments)
            pos.append(eachBase.reference_pos)

            POS = np.asarray(pos)
            posRegion = np.arange(left-1, right)
            #idx = np.in1d(posRegion, POS).nonzero()[0]
            idx = np.in1d(posRegion, POS)
            covSegment[idx] = cov

        covRegion = np.hstack([covRegion, covSegment])

    bamFile.close()

    return covRegion

def get_diff(covRegion_A, covRegion_B, cutoff):

    const = 1.0
    idxUsed = np.invert(np.logical_and(covRegion_A < cutoff, covRegion_B < cutoff))
    covRegionFoldCh = np.zeros_like(covRegion_A, dtype=float)
    covRegionFoldCh[idxUsed] = (covRegion_B[idxUsed] + const) / (covRegion_A[idxUsed] + const)

    if np.any(idxUsed):
        covRegionDiff = np.empty_like(covRegionFoldCh, dtype=str)

        # Position with coverage below cutoff
        covRegionDiff[covRegionFoldCh == 0.0] = '.'

        # Position without coverage change
        covRegionDiff[covRegionFoldCh == 1.0] = ':'

        # Position with coverage increased
        covRegionDiff[np.logical_and(covRegionFoldCh > 1.0, covRegionFoldCh <= 6.0/5)] = '5'
        covRegionDiff[np.logical_and(covRegionFoldCh > 6.0/5, covRegionFoldCh <= 8.0/5)] = '6'
        covRegionDiff[np.logical_and(covRegionFoldCh > 8.0/5, covRegionFoldCh <= 10.0/5)] = '7'
        covRegionDiff[np.logical_and(covRegionFoldCh > 10.0/5, covRegionFoldCh <= 12.0/5)] = '8'
        covRegionDiff[covRegionFoldCh > 12.0/5] = '9'

        # Position with coverage deceased
        covRegionDiff[np.logical_and(covRegionFoldCh >= 5.0/6, covRegionFoldCh < 1.0)] = '4'
        covRegionDiff[np.logical_and(covRegionFoldCh >= 5.0/8, covRegionFoldCh < 5.0/6)] = '3'
        covRegionDiff[np.logical_and(covRegionFoldCh >= 5.0/10, covRegionFoldCh < 5.0/8)] = '2'
        covRegionDiff[np.logical_and(covRegionFoldCh >= 5.0/12, covRegionFoldCh < 5.0/10)] = '1'
        covRegionDiff[np.logical_and(covRegionFoldCh > 0.0, covRegionFoldCh < 5.0/12)] = '0'

        covRegionDiff = ''.join(covRegionDiff.tolist())
    else:
        covRegionDiff = '.' * covRegion_A.size

    return covRegionDiff

def main():

    genomeFile = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.STAR.fa'
    annoPklz = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/annotation/Homo_sapiens.GRCh37.75.primary_assembly.anno.pklz'
    bamDir = '/cbio/grlab/projects/GuidoWendel/Silvestrol/results/alignments/bam/'
    geneList = '/cbio/grlab/projects/RibosomeFootprint/Lists/SilTEup.prot.5p.txt'
    seqType = 'UTR5'
    seqFeature = 'ALL'
    outputFile = '/cbio/grlab/projects/RibosomeFootprint/CovChange/SilTEup.prot.5p.utr5.all.cutoff3.covdiff'

    cutoff = 3.0

    bamFilesCondA = [bamDir + 'all_Control_1_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam', bamDir + 'all_Control_2_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam']
    bamFilesCondB = [bamDir + 'all_Silvestiol_1_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam', bamDir + 'all_Silvestiol_3_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam']

    sizeFactorsA = [1.134124, 0.802650]
    sizeFactorsB = [1.098537, 0.958124]

    genome = read_genome(genomeFile)

    FileIn = gzip.GzipFile(annoPklz, 'rb')
    (gene_id, gene_name, gene_biotype, chr, gene_region, gene_strand, transc_id, transc_type, transc, cds, utr5, utr3) = pickle.load(FileIn)
    FileIn.close()
    print 'Load annotation: done.'

    FileIn = open(geneList, 'r')
    FileOut = open(outputFile, 'w')

    cnt = 0
    for line in FileIn:

        sys.stdout.flush()
        print '\r%i genes finished ...' % cnt ,

        if not line.lstrip().startswith('ENSG'):
            continue
        geneID = line.strip().split('\t')[0]

        if geneID not in chr:
            print '\n%s dose not exists in annotation file.' % geneID
            continue

        Chr = chr[geneID]
        Strand = gene_strand[geneID]

        if seqType == 'UTR5':
            Obj = utr5[geneID]
        elif seqType == 'UTR3':
            Obj = utr3[geneID]
        elif seqType == 'CDS':
            Obj = cds[geneID]
        elif seqType == 'TRANSCRIPT':
            Obj = transc[geneID]
        else:
            pass

        if seqFeature == 'ALL':
            idx = range(len(Obj))
        elif seqFeature == 'LONGEST':
            idx = find_len_extreme(Obj, seqFeature)
        else:
            pass

        for i in idx:
            if Obj[i].shape[0] == 0:
                continue

            transcID = transc_id[geneID][i]
            geneName = gene_name[geneID]

            seq = get_seq(genome, Obj[i], Chr, Strand)

            covRegion_A = np.zeros(1)
            for k in range(len(bamFilesCondA)):
                covRegion_A_oneBam = get_cov(bamFilesCondA[k], Obj[i], Chr)
                covRegion_A_oneBam = covRegion_A_oneBam / sizeFactorsA[k]
                covRegion_A = covRegion_A + covRegion_A_oneBam
            covRegion_A = covRegion_A / len(bamFilesCondA)
            if Strand == '-':
                covRegion_A = covRegion_A[::-1]

            covRegion_B = np.zeros(1)
            for k in range(len(bamFilesCondB)):
                covRegion_B_oneBam = get_cov(bamFilesCondB[k], Obj[i], Chr)
                covRegion_B_oneBam = covRegion_B_oneBam / sizeFactorsB[k]
                covRegion_B = covRegion_B + covRegion_B_oneBam
            covRegion_B = covRegion_B / len(bamFilesCondB)
            if Strand == '-':
                covRegion_B = covRegion_B[::-1]

            covRegionDiff = get_diff(covRegion_A, covRegion_B, cutoff)

            FileOut.write('>%s|%s|%s\n' % (geneID, transcID, geneName))
            FileOut.write('%s\n' % seq)
            FileOut.write('%s\n' % covRegionDiff)

        cnt += 1

    FileOut.close()
    FileIn.close()

if __name__ == '__main__':

   main()
