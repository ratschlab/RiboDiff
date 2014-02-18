import sys
import numpy as np
import readInputs as ri

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

class normalizeLibrarySize(object):

    def CalculateLibrarySize(countNdarray):
        geoMeans = np.exp(np.mean(np.log(countNdarray), axis=1))
        librarySizes = np.median(countNdarray / geoMeans, axis=0)
        return librarySizes

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        FileIn = ri.readInputs(sys.argv[1], sys.argv[2])
        experRF, experRNA = FileIn.ParseExper()
        data = FileIn.ReadCount(experRF, experRNA)
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
        nls = normalizeLibrarySize()
        print type(data.countRibo)
        data.libSizesRibo = nls.CalculateLibrarySize(data.countRibo)
        data.libSizesRNA = nls.CalculateLibrarySize(data.countRNA)
        print 'Library size:\nribo footprint: %f\nRNA seq: %f' % data.libSizesRibo, data.libSizesRNA
