#!/usr/bin/env python
"""
Loading and preparing the data.
"""

import re
import sys
import numpy as np

class LoadInputs(object):
    """ Read the experiment description file, and use it to guide reading gene count file. """

    def __init__(self, opts):
        self.fileNameExper = opts.exptOutline
        self.fileNameCount = opts.cntFile
        self.experiment= np.empty([1, 1], dtype='str')
        self.exper     = np.empty([1, 1], dtype='str')
        self.experRibo = np.empty([1, 1], dtype='str')
        self.experRna  = np.empty([1, 1], dtype='str')
        self.experCtl  = np.empty([1, 1], dtype='str')
        self.experTrt  = np.empty([1, 1], dtype='str')
        self.nameConA  = np.empty([1, 1], dtype='str')
        self.nameConB  = np.empty([1, 1], dtype='str')
        self.idxRibo = np.empty([1, 1], dtype='int')
        self.idxRna  = np.empty([1, 1], dtype='int')
        self.idxCtl  = np.empty([1, 1], dtype='int')
        self.idxTrt  = np.empty([1, 1], dtype='int')
        self.geneIDs = np.empty([1, 1], dtype='str')
        self.countRibo = np.empty([1, 1], dtype='int')
        self.countRna  = np.empty([1, 1], dtype='int')
        self.headerRibo = np.empty([1, 1], dtype='str')
        self.headerRna  = np.empty([1, 1], dtype='str')
        self.libSizesRibo = np.empty([1, 1], dtype='float')
        self.libSizesRna  = np.empty([1, 1], dtype='float')
        self.matrix = np.empty([1, 1], dtype='int')
        self.muRaw     = np.empty([1, 1], dtype='float')
        self.muRawRibo = np.empty([1, 1], dtype='float')
        self.muRawRna  = np.empty([1, 1], dtype='float')
        self.dispRaw     = np.empty([1, 1], dtype='float')
        self.dispRawRibo = np.empty([1, 1], dtype='float')
        self.dispRawRna  = np.empty([1, 1], dtype='float')
        self.dispRawConv = np.empty([1, 1], dtype='bool')
        self.dispRawMthd = np.empty([1, 1], dtype='str')
        self.dispFitted     = np.empty([1, 1], dtype='float')
        self.dispFittedRibo = np.empty([1, 1], dtype='float')
        self.dispFittedRna  = np.empty([1, 1], dtype='float')
        self.dispFittedIdx  = np.empty([1, 1], dtype='int')
        self.dispFittedRiboIdx = np.empty([1, 1], dtype='int')
        self.dispFittedRnaIdx  = np.empty([1, 1], dtype='int')
        self.Lambda     = np.empty([1, 1], dtype='float')
        self.LambdaRibo = np.empty([1, 1], dtype='float')
        self.LambdaRna  = np.empty([1, 1], dtype='float')
        self.muAdj     = np.empty([1, 1], dtype='float')
        self.muAdjRibo = np.empty([1, 1], dtype='float')
        self.muAdjRna  = np.empty([1, 1], dtype='float')
        self.dispAdj     = np.empty([1, 1], dtype='float')
        self.dispAdjRibo = np.empty([1, 1], dtype='float')
        self.dispAdjRna  = np.empty([1, 1], dtype='float')
        self.dispAdjConv = np.empty([1, 1], dtype='bool')
        self.dispAdjMthd = np.empty([1, 1], dtype='str')
        self.pval  = np.empty([1, 1], dtype='float')
        self.padj  = np.empty([1, 1], dtype='float')
        self.TEctl = np.empty([1, 1], dtype='float')
        self.TEtrt = np.empty([1, 1], dtype='float')
        self.logFoldChangeTE = np.empty([1, 1], dtype='float')
        self.dispDiff = opts.dispDiff

    def parse_expt(self):
        """ Read the experiment description file """

        self.experiment = np.loadtxt(self.fileNameExper, dtype=str, delimiter=',', skiprows=1)
        self.exper = self.experiment.copy()

        seqType = np.unique(self.exper[:, 1])
        if seqType.size != 2:
            sys.stderr.write('Error: only two types of sequencing keyword are allowed. Please check the second column in experimental outline file.\n')
            sys.exit()
        for eachSeqType in seqType.tolist():
            if re.search(r'^Ribo-Seq$', eachSeqType, re.I):
                idxRibo = self.exper[:, 1] == eachSeqType
            if re.search(r'^RNA-Seq$', eachSeqType, re.I):
                idxRna = self.exper[:, 1] == eachSeqType
        try:
            idxRibo
        except NameError:
            sys.stderr.write('Error: in the second column in experimental outline file, please use \'Ribo-Seq\' to indicate which replicates are ribosome profiling data.\n')
            sys.exit()

        try:
            idxRna
        except NameError:
            sys.stderr.write('Error: in the second column in experimental outline file, please use \'RNA-Seq\' to indicate which replicates are RNA-Seq data.\n')
            sys.exit()

        condition, condIdx = np.unique(self.exper[:, 2], return_index=True)
        if condition.size != 2:
            sys.stderr.write('Error: only two conditions are allowed. Please check the last column in experimental outline file.\n')
            sys.exit()
        if condIdx[0] == 0:
            idxCtl = self.exper[:, 2] == condition[0]
            idxTrt = self.exper[:, 2] == condition[1]
            self.nameCondA = condition[0]
            self.nameCondB = condition[1]
        else:
            idxCtl = self.exper[:, 2] == condition[1]
            idxTrt = self.exper[:, 2] == condition[0]
            self.nameCondA = condition[1]
            self.nameCondB = condition[0]

        self.exper[idxRibo,1] = 'Ribo'
        self.exper[idxRna, 1] = 'mRna'
        self.exper[idxCtl, 2] = 'ConditionA'
        self.exper[idxTrt, 2] = 'ConditionB'

        self.experRibo = self.exper[idxRibo,0]
        self.experRna  = self.exper[idxRna, 0]
        self.experCtl  = self.exper[idxCtl, 0]
        self.experTrt  = self.exper[idxTrt, 0]

        return self

    def read_count(self):
        """ Load the discrete count file """

        with open(self.fileNameCount, 'r') as FileIn:
            header = np.array(FileIn.readline().strip().split('\t'), dtype=str)

        idxRibo = np.in1d(header, self.experRibo).nonzero()[0]
        idxRna  = np.in1d(header, self.experRna ).nonzero()[0]
        idxCtl  = np.in1d(header, self.experCtl ).nonzero()[0]
        idxTrt  = np.in1d(header, self.experTrt ).nonzero()[0]

        if idxRibo.size != self.experRibo.size or idxRna.size != self.experRna.size or idxCtl.size != self.experCtl.size or idxTrt.size != self.experTrt.size:
            sys.stderr.write('Error: At least one sample or replicate\'s name in count file and experimental outline file does not match.\n')
            sys.exit()

        idxRiboCtl = np.intersect1d(idxRibo, idxCtl)
        idxRiboTrt = np.intersect1d(idxRibo, idxTrt)
        idxRnaCtl  = np.intersect1d(idxRna,  idxCtl)
        idxRnaTrt  = np.intersect1d(idxRna,  idxTrt)

        self.headerRibo = header[np.hstack([idxRiboCtl, idxRiboTrt])]
        self.headerRna  = header[np.hstack([idxRnaCtl,  idxRnaTrt ])]

        geneIDs = np.loadtxt(self.fileNameCount, dtype=str, skiprows=1, usecols=(0,))
        self.geneIDs = geneIDs.reshape(geneIDs.size, 1)

        countRiboCtl = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRiboCtl)
        countRiboTrt = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRiboTrt)
        if idxRiboCtl.size == 1:
            countRiboCtl = countRiboCtl.reshape(countRiboCtl.size, 1)
        if idxRiboTrt.size == 1:
            countRiboTrt = countRiboTrt.reshape(countRiboTrt.size, 1)
        self.countRibo = np.hstack([countRiboCtl, countRiboTrt])

        countRnaCtl = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRnaCtl)
        countRnaTrt = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRnaTrt)
        if idxRnaCtl.size == 1:
            countRnaCtl = countRnaCtl.reshape(countRnaCtl.size, 1)
        if idxRnaTrt.size == 1:
            countRnaTrt = countRnaTrt.reshape(countRnaTrt.size, 1)
        self.countRna = np.hstack([countRnaCtl, countRnaTrt])

        self.idxRibo = np.arange(self.countRibo.shape[1])
        self.idxRna  = np.arange(self.countRna.shape[1]) + self.idxRibo.size
        self.idxCtl  = np.hstack([np.arange(idxRiboCtl.size), np.arange(idxRnaCtl.size) + self.idxRibo.size])
        self.idxTrt  = np.hstack([np.arange(idxRiboTrt.size) + idxRiboCtl.size, np.arange(idxRnaTrt.size) + self.idxRibo.size + idxRnaCtl.size])

        return self

