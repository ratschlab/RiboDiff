#!/usr/bin/env python

import sys
import numpy as np

class LoadInputs(object):

    """ Read the experiment description file, and use it to guide reading gene count file. 
        Store all data as a class object containing relative attributes """

    def __init__(self, opts):
        self.fileNameExper = opts.exptOutline
        self.fileNameCount = opts.cntFile
        self.experiment= np.empty([1, 1], dtype='str')
        self.exper     = np.empty([1, 1], dtype='str')
        self.experRibo = np.empty([1, 1], dtype='str')
        self.experRna  = np.empty([1, 1], dtype='str')
        self.experCtl  = np.empty([1, 1], dtype='str')
        self.experTrt  = np.empty([1, 1], dtype='str')
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
        self.dispFittedConv     = np.empty([1, 1], dtype='bool')
        self.dispFittedRiboConv = np.empty([1, 1], dtype='bool')
        self.dispFittedRnaConv  = np.empty([1, 1], dtype='bool')
        self.beta     = np.empty([1, 1], dtype='float')
        self.betaRibo = np.empty([1, 1], dtype='float')
        self.betaRna  = np.empty([1, 1], dtype='float')
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

        # add codes to make it case insensative
        idxRibo = self.exper[:, 1] == 'Ribo-Seq'
        idxRna  = self.exper[:, 1] == 'RNA-Seq'
        idxCtl  = self.exper[:, 2] == 'ConditionA'
        idxTrt  = self.exper[:, 2] == 'ConditionB'

        self.exper[idxRibo,1] = 'Ribo'
        self.exper[idxRna, 1] = 'mRna'
        self.exper[idxCtl, 2] = 'Control'
        self.exper[idxTrt, 2] = 'Treated'
        #

        self.experRibo = self.exper[idxRibo,0]
        self.experRna  = self.exper[idxRna, 0]
        self.experCtl  = self.exper[idxCtl, 0]
        self.experTrt  = self.exper[idxTrt, 0]

        return self

    def read_count(self):

        """ Load the discrete count file """

        with open(self.fileNameCount, 'r') as FileIn:
            header = np.array(FileIn.readline().strip().split('\t'), dtype=str)

        # add codes to dectect if the column name in the header of the count file agree with that in Experiment Outline File 
        idxRibo = np.in1d(header, self.experRibo).nonzero()[0]
        idxRna  = np.in1d(header, self.experRna ).nonzero()[0]
        idxCtl  = np.in1d(header, self.experCtl ).nonzero()[0]
        idxTrt  = np.in1d(header, self.experTrt ).nonzero()[0]
        #

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
        self.countRibo = np.hstack([countRiboCtl, countRiboTrt])

        countRnaCtl = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRnaCtl)
        countRnaTrt = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRnaTrt)
        self.countRna = np.hstack([countRnaCtl, countRnaTrt])

        self.idxRibo = np.arange(self.countRibo.shape[1])
        self.idxRna  = np.arange(self.countRna.shape[1]) + self.idxRibo.size
        self.idxCtl  = np.hstack([np.arange(idxRiboCtl.size), np.arange(idxRnaCtl.size) + self.idxRibo.size])
        self.idxTrt  = np.hstack([np.arange(idxRiboTrt.size) + idxRiboCtl.size, np.arange(idxRnaTrt.size) + self.idxRibo.size + idxRnaCtl.size])

        return self

