#!/usr/bin/env python

import sys
import numpy as np

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python loadinput.py Experiment_Outline_File Count_File' + '\n')

class LoadInputs(object):

    """ Read the experiment description file, and use it to guide to read gene count file 
        and store all data as a class object containing relative attributes """

    def __init__(self, opts):
        self.fileNameExper = opts.exptOutline
        self.fileNameCount = opts.cntFile
        self.experiment = np.empty([1, 1], dtype='str')
        self.exper = np.empty([1, 1], dtype='str')
        self.experRibo = np.empty([1, 1], dtype='str')
        self.experRNA  = np.empty([1, 1], dtype='str')
        self.experCtl = np.empty([1, 1], dtype='str')
        self.experTrt = np.empty([1, 1], dtype='str')
        self.idxRibo = np.empty([1, 1], dtype='int')
        self.idxRNA  = np.empty([1, 1], dtype='int')
        self.idxCtl = np.empty([1, 1], dtype='int')
        self.idxTrt = np.empty([1, 1], dtype='int')
        self.geneIDs = np.empty([1, 1], dtype='str')
        self.countRibo = np.empty([1, 1], dtype='int')
        self.countRNA  = np.empty([1, 1], dtype='int')
        self.headerRibo = np.empty([1, 1], dtype='str')
        self.headerRNA  = np.empty([1, 1], dtype='str')
        self.libSizesRibo = np.empty([1, 1], dtype='float')
        self.libSizesRNA  = np.empty([1, 1], dtype='float')
        self.matrix = np.empty([1, 1], dtype='int')
        self.disperRaw  = np.empty([1, 1], dtype='float')
        self.disperRawRibo = np.empty([1, 1], dtype='float')
        self.disperRawRNA  = np.empty([1, 1], dtype='float')
        self.disperRawConv = np.empty([1, 1], dtype='bool')
        self.disperRawMthd = np.empty([1, 1], dtype='str')
        self.disperFitted = np.empty([1, 1], dtype='float')
        self.disperFittedRibo = np.empty([1, 1], dtype='float')
        self.disperFittedRNA  = np.empty([1, 1], dtype='float')
        self.disperFittedIdx = np.empty([1, 1], dtype='int')
        self.disperFittedRiboIdx = np.empty([1, 1], dtype='int')
        self.disperFittedRNAIdx  = np.empty([1, 1], dtype='int')
        self.disperFittedConv = np.empty([1, 1], dtype='bool')
        self.disperFittedRiboConv = np.empty([1, 1], dtype='bool')
        self.disperFittedRNAConv  = np.empty([1, 1], dtype='bool')
        self.beta = np.empty([1, 1], dtype='float')
        self.betaRibo = np.empty([1, 1], dtype='float')
        self.betaRNA  = np.empty([1, 1], dtype='float')
        self.disperAdj = np.empty([1, 1], dtype='float')
        self.disperAdjRibo = np.empty([1, 1], dtype='float')
        self.disperAdjRNA  = np.empty([1, 1], dtype='float')
        self.disperAdjConv = np.empty([1, 1], dtype='bool')
        self.disperAdjMthd = np.empty([1, 1], dtype='str')
        self.pval = np.empty([1, 1], dtype='float')
        self.padj = np.empty([1, 1], dtype='float')
        self.TEctl = np.empty([1, 1], dtype='float')
        self.TEtrt = np.empty([1, 1], dtype='float')
        self.logFoldChangeTE = np.empty([1, 1], dtype='float')
        self.dispDiff = opts.dispDiff
        
    def parse_expt(self):

        """ Read the experiment description file """

        self.experiment = np.loadtxt(self.fileNameExper, dtype=str, delimiter=',', skiprows=1)
        self.exper = self.experiment.copy()

        # add codes to make it case insensative
        idxRibo = self.exper[:, 1] == 'Ribosome_Footprint'
        idxRNA  = self.exper[:, 1] == 'RNA_seq'
        idxCtl  = self.exper[:, 2] == 'Control'
        idxTrt  = self.exper[:, 2] == 'Treated'

        self.exper[idxRibo,1] = 'Ribo'
        self.exper[idxRNA, 1] = 'mRna'
        self.exper[idxCtl, 2] = 'Control'
        self.exper[idxTrt, 2] = 'Treated'
        #

        self.experRibo = self.exper[idxRibo,0]
        self.experRNA  = self.exper[idxRNA, 0]
        self.experCtl  = self.exper[idxCtl, 0]
        self.experTrt  = self.exper[idxTrt, 0]

        return self

    def read_count(self):

        """ Read the discrete reads count file """

        with open(self.fileNameCount, 'r') as FileIn:
            header = np.array(FileIn.readline().rstrip().split('\t'), dtype=str)

        # add codes to dectect if the column name in the header of the count file agree with that in Experiment Outline File 
        idxRibo = np.in1d(header, self.experRibo).nonzero()[0]
        idxRNA  = np.in1d(header, self.experRNA ).nonzero()[0]
        idxCtl  = np.in1d(header, self.experCtl ).nonzero()[0]
        idxTrt  = np.in1d(header, self.experTrt ).nonzero()[0]
        #

        idxRiboCtl = np.intersect1d(idxRibo, idxCtl)
        idxRiboTrt = np.intersect1d(idxRibo, idxTrt)
        idxRNACtl  = np.intersect1d(idxRNA,  idxCtl)
        idxRNATrt  = np.intersect1d(idxRNA,  idxTrt)

        self.headerRibo = header[np.hstack([idxRiboCtl, idxRiboTrt])]
        self.headerRNA  = header[np.hstack([idxRNACtl,  idxRNATrt ])]

        geneIDs = np.loadtxt(self.fileNameCount, dtype=str, skiprows=1, usecols=(0,))
        self.geneIDs = geneIDs.reshape(len(geneIDs), 1)

        countRiboCtl = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRiboCtl)
        countRiboTrt = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRiboTrt)
        self.countRibo = np.hstack([countRiboCtl, countRiboTrt])

        countRNACtl = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRNACtl)
        countRNATrt = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRNATrt)
        self.countRNA = np.hstack([countRNACtl, countRNATrt])

        self.idxRibo = np.arange(self.countRibo.shape[1])
        self.idxRNA  = np.arange(self.countRNA.shape[1]) + len(self.idxRibo)
        self.idxCtl  = np.hstack([np.arange(len(idxRiboCtl)), np.arange(len(idxRNACtl)) + len(self.idxRibo)])
        self.idxTrt  = np.hstack([np.arange(len(idxRiboTrt)) + len(idxRiboCtl), np.arange(len(idxRNATrt)) + len(self.idxRibo) + len(idxRNACtl)])

        return self

if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
    else:
        print '*'*25
        FileIn = LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_expt()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
