#!/usr/bin/env python

import sys
import pysam
import os
import cPickle as pickle
import gzip

def get_rRNA_reads(fileInName, fileOutName):

    '''
    This script filters out the ambiguous read alignment to refine the rRNA read.
    '''

    bamFile = pysam.Samfile(fileInName, 'rb')

    readIDs = set()

    for eachRead in bamFile:

        if eachRead.is_unmapped:
            continue

        alignedLen = eachRead.qlen
        if alignedLen < 15:
            continue
        
        if not eachRead.is_reverse and eachRead.cigar[-1][0] == 4:
            clip3EndLen = eachRead.cigar[-1][1]
        elif eachRead.is_reverse and eachRead.cigar[0][0] == 4:
            clip3EndLen = eachRead.cigar[0][1]
        else:
            clip3EndLen = 0

        if alignedLen >= 15 and alignedLen <= 29 and alignedLen + clip3EndLen >= 47:
            readIDs.add(eachRead.qname)
        elif alignedLen >= 30 and alignedLen <= 50 and alignedLen + clip3EndLen >= 42:
            readIDs.add(eachRead.qname)
        else:
            pass

    bamFile.close()

    print '%i rRNA reads' % len(readIDs)

    FileOut = gzip.GzipFile(fileOutName, 'wb')
    pickle.dump(readIDs, FileOut, pickle.HIGHEST_PROTOCOL)
    FileOut.close()

if __name__ == '__main__':

    if len(sys.argv) != 2:
        sys.stderr.write('Usage:' + '\n' + 'python get_rRNA_reads.py <input Bam> <output file>' + '\n')
        sys.exit()

    get_rRNA_reads(sys.argv[1], sys.argv[2])
