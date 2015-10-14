#!/usr/bin/env python

import sys
import os
import cPickle as pickle
import gzip
import pysam
import string

def remove_rRNA(fileInNameBam, fileInNamePickle, fileOutName):

    if os.path.exists(fileOutName):
        print '%s already exists!' % fileOutName
        return

    if os.path.exists(fileInNamePickle):
        FileIn = gzip.GzipFile(fileInNamePickle, 'rb')
        readIDs = pickle.load(FileIn)
        FileIn.close()

        bamFile = pysam.Samfile(fileInNameBam, 'rb')
        FileOut = pysam.Samfile(fileOutName, 'wb', template=bamFile)

        for eachRead in bamFile:
            if not eachRead.qname in readIDs:
                FileOut.write(eachRead)

        FileOut.close()
        bamFile.close()

        print 'Remove rRNA, %s: Done.' % fileOutName
    else:
        print '%s does not exit!' % fileInNamePickle
        return

def remove_ambiguous_read(fileInName, fileOutName):

    if os.path.exists(fileOutName):
        print '%s already exists!' % fileOutName
        return

    bamFile = pysam.Samfile(fileInName, 'rb')
    FileOut = pysam.Samfile(fileOutName, 'wb', template=bamFile)

    adapter = 'CTGTAGGC'

    for eachRead in bamFile:

        if eachRead.is_unmapped:
            continue

        alignedLen = eachRead.qlen

        if not eachRead.is_reverse and eachRead.cigar[-1][0] == 4:
            clip3EndLen = eachRead.cigar[-1][1]
            clip3EndSeq = eachRead.seq[-clip3EndLen:]
        elif eachRead.is_reverse and eachRead.cigar[0][0] == 4:
            clip3EndLen = eachRead.cigar[0][1]
            clip3EndSeq = eachRead.seq[:clip3EndLen]
            clip3EndSeq = clip3EndSeq[::-1].translate(string.maketrans('ATGCNatgcn', 'TACGNtacgn'))
        else:
            clip3EndLen = 0
            clip3EndSeq = ''

        if alignedLen < 15:
            continue
        elif alignedLen >= 15 and alignedLen <= 29 and alignedLen + clip3EndLen < 47:
            continue
        elif alignedLen >= 30 and alignedLen <= 50 and alignedLen + clip3EndLen < 42:
            continue
        else:
            pass

        if clip3EndLen != 0:
            for i in range(2):
                winSeq = clip3EndSeq[i:i+8]
                mismatch = 0
                for j in range(len(winSeq)):
                    if not winSeq[j] == adapter[j]:
                        mismatch += 1
                    if mismatch > 1:
                        break
                if mismatch <= 1:
                    FileOut.write(eachRead)
                    break
        else:
            FileOut.write(eachRead)

    FileOut.close()
    bamFile.close()

    print 'Remove ambiguous alignment, %s: Done.' % fileOutName

def remove_extreme_length(fileInName, fileOutName):

    if os.path.exists(fileOutName):
        print '%s already exists!' % fileOutName
        return

    bamFile = pysam.Samfile(fileInName, 'rb')
    FileOut = pysam.Samfile(fileOutName, 'wb', template=bamFile)

    for eachRead in bamFile:
        if eachRead.is_unmapped:
            continue

        alignedLen = eachRead.qlen

        if alignedLen >= 20 and alignedLen <= 40:
            FileOut.write(eachRead)

    FileOut.close()
    bamFile.close()

    print 'Select read by length, %s: Done.' % fileOutName 

if __name__ == '__main__':

    originalBAM = sys.argv[1]
    pickleFile = sys.argv[2]
    fileOutName1 = originalBAM.replace('.bam', '.NOrRNA.bam'
    remove_rRNA(originalBAM, pickleFile, fileOutName1)

    fileInName = fileOutName1
    fileOutName2 = originalBAM.replace('.bam', '.NOrRNA.NoAmbig.bam')
    remove_ambiguous_read(fileInName, fileOutName2)

    fileInName = fileOutName2
    fileOutName3 = originalBAM.replace('.bam', '.NOrRNA.NoAmbig.NoExt.bam')
    remove_extreme_length(fileInName, fileOutName3)
