#!/usr/bin/env python

import sys
import os
import cPickle as pickle
import gzip
import pysam

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
        print '%s does not exist!' % fileInNamePickle
        return

if __name__ == '__main__':

    originalBAM = sys.argv[1]
    pickleFile = sys.argv[2]
    fileOutName = originalBAM.replace('.bam', '.NOrRNA.bam')
    remove_rRNA(originalBAM, pickleFile, fileOutName)
