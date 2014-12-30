#!/usr/bin/env python

import sys
import numpy as np
import cPickle as pickle 
import gzip
import re
import pdb

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
                globalPos = np.array([transcript[i, 0] + localPosStart - lengthBef - 1, transcript[i, 0] + localPosEnd - lengthBef - 1])
                break
            elif localPosStart >= exonLeft and localPosStart <= exonRight and localPosEnd > exonRight:
                globalPos = np.array([transcript[i, 0] + localPosStart - lengthBef - 1, transcript[i, 1]])
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
                globalPos = np.array([transcript[i, 1] - localPosEnd + lengthBef + 1, transcript[i, 1] - localPosStart + lengthBef + 1])
                break
            elif localPosStart >= exonRight and localPosStart <= exonLeft and localPosEnd > exonLeft:
                globalPos = np.array([transcript[i, 0], transcript[i, 1] - localPosStart + lengthBef + 1])
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

def get_unique_global_pos(redundantPos):

    uniquePos = {}

    numGQ = 0
    for eachKey in redundantPos:
        if len(redundantPos[eachKey]) != 0:
            uniqueListLists = []
            uniqueListArray = []
            for eachRdPos in redundantPos[eachKey]:
                if eachRdPos.tolist() not in uniqueListLists:
                    uniqueListLists.append(eachRdPos.tolist())
                    uniqueListArray.append(eachRdPos)

            uniquePos[eachKey] = uniqueListArray
            numGQ = numGQ + len(uniqueListArray)
        else:
            uniquePos[eachKey] = []

    numGene = len(uniquePos)

    return (uniquePos, numGene, numGQ)

def main():

    annoPklz = '/cbio/grlab/share/databases/genomes/H_sapiens/Homo_sapiens.GRCh37.ENSEMBL75/annotation/Homo_sapiens.GRCh37.75.primary_assembly.anno.pklz'
    covDiffFile = '/cbio/grlab/projects/RibosomeFootprint/CovChange/SilTEdown.prot.5p.utr5.all.cutoff3.covdiff'
    gqPklz = '/cbio/grlab/projects/RibosomeFootprint/GQ/SilTEdown.prot.5p.utr5.all.cutoff3.gq'

    FileIn = gzip.GzipFile(annoPklz, 'rb')
    (gene_id, gene_name, gene_biotype, chr, gene_region, gene_strand, transc_id, transc_type, transc, cds, utr5, utr3) = pickle.load(FileIn)
    FileIn.close()
    print 'Load annotation done.'

    regSeq = re.compile(r'(G{2,})([TCGA]{0,10}?\1){3}', re.IGNORECASE)
    regSplit = re.compile(r'[TCA]{1,}', re.IGNORECASE)
    regDiff = re.compile(r'[6-9]{10,}')

    GENE_ID = []
    CHR = {}
    GENE_STRAND = {}
    GQ_POS = {}

    FileIn = open(covDiffFile, 'r')
    cnt = 1
    for line in FileIn:
        if line.lstrip().startswith('>'):
            geneID = line.strip().split('|')[0].replace('>', '')
            transcID = line.strip().split('|')[1]
            geneName = line.strip().split('|')[2]
            if geneID not in GQ_POS:

                sys.stdout.flush()
                print '\rStart to process gene %i ...' % cnt ,
                cnt += 1

                GQ_POS[geneID] = []
                GENE_ID.append(geneID)
                CHR[geneID] = chr[geneID]
                GENE_STRAND[geneID] = gene_strand[geneID]

            seq = next(FileIn)
            cov = next(FileIn)

            iterSeq = regSeq.finditer(seq)
            posLocSeq = []
            for matchSeq in iterSeq:
                if len(regSplit.split(matchSeq.group())) >= 3 and len(matchSeq.group()) <= 50:
                    posLocSeq.append(matchSeq.span())

            iterDiff = regDiff.finditer(cov)
            for matchDiff in iterDiff:
                for eachPosLocSeq in posLocSeq:
                    if abs(matchDiff.start() - eachPosLocSeq[0]) <= 10: 
                        globalPos = local_transc_pos_2_global(transcID, transc_id[geneID], transc[geneID], gene_strand[geneID], eachPosLocSeq[0]+1, eachPosLocSeq[1])
                        GQ_POS[geneID].append(globalPos)
                        break

    FileIn.close()

    print '\nStart to remove repeated global positions ...'
    GQ_POS, numGene, numGQ = get_unique_global_pos(GQ_POS)

    genes = (GENE_ID, CHR, GENE_STRAND, GQ_POS)
    FileOut = gzip.GzipFile(gqPklz, 'wb')
    pickle.dump(genes, FileOut, pickle.HIGHEST_PROTOCOL)
    FileOut.close()
    print 'Saving the result: Done!'

    print '%i Genes\n%i GQs' % (numGene, numGQ)

if __name__ == '__main__':

    main()
