import sys
import numpy as np

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

class readInputs(object):

    """ Read the experiment description file, and use it to guide to read gene count file 
        and store all data as a class object containing relative attributes """

    def __init__(self, fileName1, fileName2):
        self.fileNameExper = fileName1
        self.fileNameCount = fileName2
        self.exper = np.empty([1, 1], dtype=str)
        self.geneIDs = np.empty([1, 1], dtype=str)
        self.countRibo = np.empty([1, 1], dtype=int)
        self.countRNA = np.empty([1, 1], dtype=int)

    def ParseExper(self):

        """ Read the experiment description file """

        self.exper = np.loadtxt(self.fileNameExper, dtype=str, delimiter=',', skiprows=1)

        # add codes to make it case insensative
        idxRF = self.exper[:, 1] == 'Ribosome_Footprint'
        idxRNA = self.exper[:, 1] != 'Ribosome_Footprint'
        #

        experRF = self.exper[idxRF, 0]
        experRNA = self.exper[idxRNA, 0]
        return (experRF, experRNA)

    def ReadCount(self, experRF, experRNA):

        """ Read the discrete reads count file """

        with open(self.fileNameCount, 'r') as FileIn:
            header = np.array(FileIn.readline().rstrip().split('\t'), dtype=str)

        # add codes to dectect if the column name in the header of the count file agree with that in Experiment Outline File 
        idxRF = np.in1d(header, experRF).nonzero()[0]
        idxRNA = np.in1d(header, experRNA).nonzero()[0]
        #

        self.geneIDs = np.loadtxt(self.fileNameCount, dtype=str, skiprows=1, usecols=(0,))
        self.countRibo = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRF)
        self.countRNA = np.loadtxt(self.fileNameCount, dtype=int, skiprows=1, usecols=idxRNA)
        return self

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        FileIn = readInputs(sys.argv[1], sys.argv[2])
        experRF, experRNA = FileIn.ParseExper()
        data = FileIn.ReadCount(experRF, experRNA)
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size

