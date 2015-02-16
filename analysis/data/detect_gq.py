#!/usr/bin/env python

import sys
import numpy as np
import cPickle as pickle 
import gzip
import re
import pysam
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

def local_transc_pos_2_global(transcID, transc_id, transc, gene_strand, localPosStart, localPosEnd):

    idx = transc_id.index(transcID)
    transcript = transc[idx]

    if gene_strand == '+':
        lengthBef = 0
        for i in range(transcript.shape[0]):
            lengthNow = lengthBef + transcript[i, 1] - transcript[i, 0] + 1
            exonLeft  = lengthBef + 1
            exonRight = lengthNow
            if localPosStart >= exonLeft and localPosStart <= exonRight and localPosEnd >= exonLeft and localPosEnd <= exonRight:
                globalPos = np.array([[transcript[i, 0] + localPosStart - lengthBef - 1, transcript[i, 0] + localPosEnd - lengthBef - 1]])
                break
            elif localPosStart >= exonLeft and localPosStart <= exonRight and localPosEnd > exonRight:
                globalPos = np.array([[transcript[i, 0] + localPosStart - lengthBef - 1, transcript[i, 1]]])
            elif localPosStart < exonLeft and localPosEnd > exonRight:
                globalPos = np.vstack([globalPos, [transcript[i, 0], transcript[i, 1]]])
            elif localPosStart < exonLeft and localPosEnd >= exonLeft and localPosEnd <= exonRight:
                globalPos = np.vstack([globalPos, [transcript[i, 0], transcript[i, 0] + localPosEnd - lengthBef - 1]])
                break
            else:
                pass
            lengthBef = lengthNow
    elif gene_strand == '-':
        lengthBef = 0
        for i in range(transcript.shape[0])[::-1]:
            lengthNow = lengthBef + transcript[i, 1] - transcript[i, 0] + 1
            exonRight = lengthBef + 1
            exonLeft  = lengthNow
            if localPosStart >= exonRight and localPosStart <= exonLeft and localPosEnd >= exonRight and localPosEnd <= exonLeft:
                globalPos = np.array([[transcript[i, 1] - localPosEnd + lengthBef + 1, transcript[i, 1] - localPosStart + lengthBef + 1]])
                break
            elif localPosStart >= exonRight and localPosStart <= exonLeft and localPosEnd > exonLeft:
                globalPos = np.array([[transcript[i, 0], transcript[i, 1] - localPosStart + lengthBef + 1]])
            elif localPosStart < exonRight and localPosEnd > exonLeft:
                globalPos = np.vstack([[transcript[i, 0], transcript[i, 1]], globalPos])
            elif localPosStart < exonRight and localPosEnd >= exonRight and localPosEnd <= exonLeft:
                globalPos = np.vstack([[transcript[i, 1] - localPosEnd + lengthBef + 1, transcript[i, 1]], globalPos])
                break
            else:
                pass
            lengthBef = lengthNow
    else:
        pass

    return globalPos

def check_and_trim_FLK(GQ_ALL_POS_LOCAL, GQ_COV_POS, FLK_COV_POS_LOCAL, WHICH):

    if FLK_COV_POS_LOCAL:
        for i in range(len(FLK_COV_POS_LOCAL)):
            if FLK_COV_POS_LOCAL[i] != 0 and FLK_COV_POS_LOCAL[i] != -1:
                eachFLKstart = FLK_COV_POS_LOCAL[i][0]
                eachFLKstop  = FLK_COV_POS_LOCAL[i][1]
                idx = 0
                for eachGQ in GQ_ALL_POS_LOCAL:
                    # Only check if FLK overlaps with an expressed GQ.
                    if GQ_COV_POS[idx] == 1:
                        if WHICH == 5 and eachFLKstop == eachGQ[1]:
                            FLK_COV_POS_LOCAL[i] = -1
                        if WHICH == 5 and eachFLKstart <= eachGQ[1] and eachFLKstop > eachGQ[1]:
                            FLK_COV_POS_LOCAL[i][0] = eachGQ[1] + 1
                        if WHICH == 3 and eachFLKstart == eachGQ[0]:
                            FLK_COV_POS_LOCAL[i] = -1
                        if WHICH == 3 and eachFLKstart < eachGQ[0] and eachFLKstop >= eachGQ[1]:
                            FLK_COV_POS_LOCAL[i][1] = eachGQ[0] - 1
                    idx += 1

    return FLK_COV_POS_LOCAL

def get_unique_global_pos(rdGQ_ALL_POS, rdGQ_ALL_SEQ, rdGQ_COV_POS, rdGQ_INC_POS, rdFLK5P_COV_POS, rdFLK3P_COV_POS, rdGQ_IS_FIRST, rdCDS_COV_MU):

    nrdGQ_ALL_POS = {}
    nrdGQ_ALL_SEQ = {}
    nrdGQ_COV_POS = {}
    nrdGQ_INC_POS = {}
    nrdFLK5P_COV_POS = {}
    nrdFLK3P_COV_POS = {}
    nrdGQ_IS_FIRST = {}
    nrdCDS_COV_MU = {}

    numGene = 0
    numGQ_ALL = 0
    numGQ_COV = 0
    numGQ_INC = 0
    for eachKey in rdGQ_ALL_POS:
        if len(rdGQ_ALL_POS[eachKey]) != 0:
            uniqueListLists = []
            uniqueListArray = []
            uniqueListSEQ = []
            uniqueListCOV = []
            uniqueListINC = []
            uniqueListFLK5P = []
            uniqueListFLK3P = []
            uniqueListFIRST = []
            uniqueListCDS = []

            idx = 0
            for eachRdPos in rdGQ_ALL_POS[eachKey]:
                if eachRdPos.tolist() not in uniqueListLists:
                    uniqueListLists.append(eachRdPos.tolist())
                    uniqueListArray.append(eachRdPos)
                    uniqueListSEQ.append(rdGQ_ALL_SEQ[eachKey][idx])
                    uniqueListCOV.append(rdGQ_COV_POS[eachKey][idx])
                    uniqueListINC.append(rdGQ_INC_POS[eachKey][idx])
                    uniqueListFLK5P.append(rdFLK5P_COV_POS[eachKey][idx])
                    uniqueListFLK3P.append(rdFLK3P_COV_POS[eachKey][idx])
                    uniqueListFIRST.append(rdGQ_IS_FIRST[eachKey][idx])
                    uniqueListCDS.append(rdCDS_COV_MU[eachKey][idx])
                else:
                    IDX = uniqueListLists.index(eachRdPos.tolist())
                    if rdGQ_IS_FIRST[eachKey][idx] == 1:
                        uniqueListFIRST[IDX] = 1
                    if rdCDS_COV_MU[eachKey][idx] > uniqueListCDS[IDX]:
                        uniqueListCDS[IDX] = rdCDS_COV_MU[eachKey][idx]

                idx += 1

            nrdGQ_ALL_POS[eachKey] = uniqueListArray
            nrdGQ_ALL_SEQ[eachKey] = uniqueListSEQ
            nrdGQ_COV_POS[eachKey] = uniqueListCOV
            nrdGQ_INC_POS[eachKey] = uniqueListINC
            nrdFLK5P_COV_POS[eachKey] = uniqueListFLK5P
            nrdFLK3P_COV_POS[eachKey] = uniqueListFLK3P
            nrdGQ_IS_FIRST[eachKey] = uniqueListFIRST
            nrdCDS_COV_MU[eachKey] = uniqueListCDS

            numGene += 1
            numGQ_ALL = numGQ_ALL + len(uniqueListArray)
            numGQ_COV = numGQ_COV + len([e for e in uniqueListCOV if e == 1])
            numGQ_INC = numGQ_INC + len([e for e in uniqueListINC if e == 1])
        else:
            nrdGQ_ALL_POS[eachKey] = []
            nrdGQ_ALL_SEQ[eachKey] = []
            nrdGQ_COV_POS[eachKey] = []
            nrdGQ_INC_POS[eachKey] = []
            nrdFLK5P_COV_POS[eachKey] = []
            nrdFLK3P_COV_POS[eachKey] = []
            nrdGQ_IS_FIRST[eachKey] = []
            nrdCDS_COV_MU[eachKey] = []

    return (nrdGQ_ALL_POS, nrdGQ_ALL_SEQ, nrdGQ_COV_POS, nrdGQ_INC_POS, nrdFLK5P_COV_POS, nrdFLK3P_COV_POS, nrdGQ_IS_FIRST, nrdCDS_COV_MU, numGene, numGQ_ALL, numGQ_COV, numGQ_INC)

def main():

    annoPklz = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/annotation/Homo_sapiens.GRCh37.75.primary_assembly.anno.pklz'
    covDiffFile = '/cbio/grlab/projects/RibosomeFootprint/CovChange/SilTEdown.prot.5p.utr5.all.cutoff3.covdiff'
    gqPklz = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTEdown.prot.5p.utr5.all.cutoff3.gq'
    FLK = 40.0

    FileIn = gzip.GzipFile(annoPklz, 'rb')
    (gene_id, gene_name, gene_biotype, chr, gene_region, gene_strand, transc_id, transc_type, transc, cds, utr5, utr3) = pickle.load(FileIn)
    FileIn.close()
    print 'Load annotation done.'

    regSeq = re.compile(r'(G{2,})([TCGA]{0,12}?\1){3,5}', re.IGNORECASE)
    regSplit = re.compile(r'[TCA]{1,}', re.IGNORECASE)
    regDiff1 = re.compile(r'[0-9:]{5,}')
    regDiff2 = re.compile(r'[6-9]{5,}')
    regDiff3 = re.compile(r'[5-9]')
    regDiff4 = re.compile(r'[0-4]')

    GENE_ID = []
    CHR = {}
    GENE_STRAND = {}
    GQ_ALL_POS = {}
    GQ_ALL_SEQ = {}
    GQ_COV_POS = {}
    GQ_INC_POS = {}
    FLK5P_COV_POS = {}
    FLK3P_COV_POS = {}
    GQ_IS_FIRST = {}
    CDS_COV_MU = {}

    numBaseINC = 0
    numBaseDEC = 0
    numBase = 0

    FileIn = open(covDiffFile, 'r')
    cnt = 1
    for line in FileIn:
        if line.lstrip().startswith('>'):
            geneID = line.strip().split('|')[0].replace('>', '')
            transcID = line.strip().split('|')[1]
            geneName = line.strip().split('|')[2]

            if geneID not in GQ_ALL_POS:

                sys.stdout.flush()
                print '\rStart to process gene %i ...' % cnt ,
                cnt += 1

                GENE_ID.append(geneID)
                CHR[geneID] = chr[geneID]
                GENE_STRAND[geneID] = gene_strand[geneID]
                GQ_ALL_POS[geneID] = []
                GQ_ALL_SEQ[geneID] = []
                GQ_COV_POS[geneID] = []
                GQ_INC_POS[geneID] = []
                FLK5P_COV_POS[geneID] = []
                FLK3P_COV_POS[geneID] = []
                GQ_IS_FIRST[geneID] = []
                CDS_COV_MU[geneID] = []

            bamDir = '/cbio/grlab/projects/GuidoWendel/Silvestrol/results/alignments/bam/'
            bamFilesCondB = [bamDir + 'all_Silvestiol_1_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam', bamDir + 'all_Silvestiol_3_Ribo.uq.sorted.rRNAfiltered.rtrimfilterd.25-35.final.bam']
            sizeFactorsB = [1.098537, 0.958124]
            idx = np.in1d(transc_id[geneID], transcID).nonzero()[0]
            cdsCov = get_cov_wrapper(cds[geneID][idx], bamFilesCondB, sizeFactorsB, chr[geneID], gene_strand[geneID])
            if cdsCov.size == 0:
                meanCdsCov = 0.0
            else:
                meanCdsCov = np.mean(cdsCov)

            seq = next(FileIn).strip()
            cov = next(FileIn).strip()

            GQ_ALL_POS_LOCAL = []
            FLK5P_COV_POS_LOCAL = []
            FLK3P_COV_POS_LOCAL = []
            IS_FIRST = 1
            iterSeq = regSeq.finditer(seq)
            for matchSeq in iterSeq:
                if len(regSplit.split(matchSeq.group())) >= 3:
                    globalPosGQ = local_transc_pos_2_global(transcID, transc_id[geneID], transc[geneID], gene_strand[geneID], matchSeq.start()+1, matchSeq.end())
                    GQ_ALL_POS[geneID].append(globalPosGQ)
                    GQ_ALL_SEQ[geneID].append(matchSeq.group())
                    GQ_ALL_POS_LOCAL.append([matchSeq.start()+1, matchSeq.end()])
                    GQ_IS_FIRST[geneID].append(IS_FIRST)
                    CDS_COV_MU[geneID].append(meanCdsCov)
                    IS_FIRST = 0

                    covDiff = cov[matchSeq.start() : matchSeq.end()]
                    if regDiff1.search(covDiff):
                        GQ_COV_POS[geneID].append(1)
                        if matchSeq.start() == 0:
                            FLK5P_COV_POS_LOCAL.append(-1)
                        else:
                            PosFLK5Pstart = max(1, matchSeq.start() - FLK + 1)
                            PosFLK5Pend = matchSeq.start()
                            FLK5P_COV_POS_LOCAL.append([PosFLK5Pstart, PosFLK5Pend])

                        if matchSeq.end() == len(seq):
                            FLK3P_COV_POS_LOCAL.append(-1)
                        else:
                            PosFLK3Pstart = matchSeq.end() + 1
                            PosFLK3Pend = min(len(seq), matchSeq.end() + FLK)
                            FLK3P_COV_POS_LOCAL.append([PosFLK3Pstart, PosFLK3Pend])
                    else:
                        GQ_COV_POS[geneID].append(0)
                        FLK5P_COV_POS_LOCAL.append(0)
                        FLK3P_COV_POS_LOCAL.append(0)

                    covDiffFirstHalf = cov[matchSeq.start() : matchSeq.start() + (matchSeq.end() - matchSeq.start()) / 2 + 1]
                    if regDiff2.search(covDiffFirstHalf):
                        GQ_INC_POS[geneID].append(1)
                    else:
                        GQ_INC_POS[geneID].append(0)

            FLK5P_COV_POS_LOCAL = check_and_trim_FLK(GQ_ALL_POS_LOCAL, GQ_COV_POS[geneID], FLK5P_COV_POS_LOCAL, 5)
            FLK3P_COV_POS_LOCAL = check_and_trim_FLK(GQ_ALL_POS_LOCAL, GQ_COV_POS[geneID], FLK3P_COV_POS_LOCAL, 3)

            for i in range(len(FLK5P_COV_POS_LOCAL)):
                if FLK5P_COV_POS_LOCAL[i] != 0 and FLK5P_COV_POS_LOCAL[i] !=  -1:
                    globalPosFLK5P = local_transc_pos_2_global(transcID, transc_id[geneID], transc[geneID], gene_strand[geneID], FLK5P_COV_POS_LOCAL[i][0], FLK5P_COV_POS_LOCAL[i][1])
                    FLK5P_COV_POS[geneID].append(globalPosFLK5P)
                else:
                    FLK5P_COV_POS[geneID].append(FLK5P_COV_POS_LOCAL[i])

            for i in range(len(FLK3P_COV_POS_LOCAL)):
                if FLK3P_COV_POS_LOCAL[i] != 0 and FLK3P_COV_POS_LOCAL[i] !=  -1:
                    globalPosFLK3P = local_transc_pos_2_global(transcID, transc_id[geneID], transc[geneID], gene_strand[geneID], FLK3P_COV_POS_LOCAL[i][0], FLK3P_COV_POS_LOCAL[i][1])
                    FLK3P_COV_POS[geneID].append(globalPosFLK3P)
                else:
                    FLK3P_COV_POS[geneID].append(FLK3P_COV_POS_LOCAL[i])

            #iterDiff = regDiff.finditer(cov)
            #for matchDiff in iterDiff:
            #    for eachPosLocSeq in posLocSeq:
            #        if abs(matchDiff.start() - eachPosLocSeq[0]) <= 10: 
            #            globalPos = local_transc_pos_2_global(transcID, transc_id[geneID], transc[geneID], gene_strand[geneID], eachPosLocSeq[0]+1, eachPosLocSeq[1])
            #            GQ_POS[geneID].append(globalPos)
            #            break

            numBaseINC = numBaseINC + len(regDiff3.findall(cov))
            numBaseDEC = numBaseDEC + len(regDiff4.findall(cov))
            numBase = numBase + len(cov)

    FileIn.close()

    print '\nStart to remove repeated global positions ...'
    GQ_ALL_POS, GQ_ALL_SEQ, GQ_COV_POS, GQ_INC_POS, FLK5P_COV_POS, FLK3P_COV_POS, GQ_IS_FIRST, CDS_COV_MU, numGene, numGQ_ALL, numGQ_COV, numGQ_INC = get_unique_global_pos(GQ_ALL_POS, GQ_ALL_SEQ, GQ_COV_POS, GQ_INC_POS, FLK5P_COV_POS, FLK3P_COV_POS, GQ_IS_FIRST, CDS_COV_MU)

    #for eachGene in GQ_INC_POS:
    #    IDX = np.nonzero(np.asarray(GQ_INC_POS[eachGene]) == 1)[0]
    #    for idx in IDX:
    #        print GQ_ALL_SEQ[eachGene][idx]

    genes = (GENE_ID, CHR, GENE_STRAND, GQ_ALL_POS, GQ_ALL_SEQ, GQ_COV_POS, GQ_INC_POS, FLK5P_COV_POS, FLK3P_COV_POS, GQ_IS_FIRST, CDS_COV_MU)
    FileOut = gzip.GzipFile(gqPklz, 'wb')
    pickle.dump(genes, FileOut, pickle.HIGHEST_PROTOCOL)
    FileOut.close()
    print 'Saving the result: Done!'

    print '%i genes contain GQ.' % numGene
    print '%i GQs in total.' % numGQ_ALL
    print '%i GQs having RF reads.' % numGQ_COV
    print '%i GQs showing increased coverage.' % numGQ_INC

    print '%i bp showing increased coverage.' % numBaseINC
    print '%i bp showing decreased coverage.' % numBaseDEC
    print '%i bp in total.' % numBase

if __name__ == '__main__':

    main()

