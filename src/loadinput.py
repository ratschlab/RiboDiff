import sys
import numpy as np

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python loadinput.py Experiment_Outline_File Gene_Count_File' + '\n')

class LoadInputs(object):

    """ Read the experiment description file, and use it to guide to read gene count file 
        and store all data as a class object containing relative attributes """

    def __init__(self, fileName1, fileName2):
        self.fileNameExper = fileName1
        self.fileNameCount = fileName2
        self.experiment = np.empty([1, 1], dtype=str)
        self.exper = np.empty([1, 1], dtype=str)
        self.experRF = np.empty([1, 1], dtype=str)
        self.experRNA = np.empty([1, 1], dtype=str)
        self.geneIDs = np.empty([1, 1], dtype=str)
        self.countRibo = np.empty([1, 1], dtype=int)
        self.countRNA = np.empty([1, 1], dtype=int)
        self.libSizesRibo = np.empty([1, 1], dtype=float)
        self.libSizesRNA = np.empty([1, 1], dtype=float)
        self.matrix = np.empty([1, 1], dtype=int)
        self.disperRaw = np.empty([1, 1], dtype=float)
        self.disperFitted = np.empty([1, 1], dtype=float)
        self.disperAdj = np.empty([1, 1], dtype=float)
        self.pval = np.empty([1, 1], dtype=float)

    def parse_exper(self):

        """ Read the experiment description file """

        self.experiment = np.loadtxt(self.fileNameExper, dtype=str, delimiter=',', skiprows=1)
        self.exper = self.experiment.copy()

        # add codes to make it case insensative
        idxRF = self.exper[:, 1] == 'Ribosome_Footprint'
        idxRNA = self.exper[:, 1] == 'RNA_seq'
        idxCtrl = self.exper[:, 1] == 'Control'
        idxTrt = self.exper[:, 1] == 'Treated'

        self.exper[idxRF, 1] = 'Ribo'
        self.exper[idxRNA, 1] = 'Rna'
        self.exper[idxCtrl, 2] = 'Control'
        self.exper[idxTrt, 2] = 'Treated'
        # 

        self.experRF = self.exper[idxRF, 0]
        self.experRNA = self.exper[idxRNA, 0]

        return self

    def read_count(self):

        """ Read the discrete reads count file """

        with open(self.fileNameCount, 'r') as FileIn:
            header = np.array(FileIn.readline().rstrip().split('\t'), dtype=str)

        # add codes to dectect if the column name in the header of the count file agree with that in Experiment Outline File 
        idxRF = np.in1d(header, self.experRF).nonzero()[0]
        idxRNA = np.in1d(header, self.experRNA).nonzero()[0]
        #

        geneIDs = np.loadtxt(self.fileNameCount, dtype=str, skiprows=1, usecols=(0,))
        self.geneIDs = geneIDs.reshape(len(geneIDs), 1)
        self.countRibo = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRF)
        self.countRNA = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRNA)

        return self

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        print '*'*25
        FileIn = LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_exper()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
