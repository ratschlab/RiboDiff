#!/usr/bin/env python

import sys
import numpy as np
import loadinput as ld

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python normalize.py Experiment_Outline_File Count_File' + '\n')

def lib_size(countNdarray):

    """ Calculate library size """

    if np.any(countNdarray < 0, axis=None):
        sys.stderr.write('Error: read count smaller than one is detected, please check input file.\n')
        sys.exit()

    countNdarrayTmp = countNdarray.copy()
    countNdarrayTmp = countNdarrayTmp + 1
    geoMeans = np.exp(np.mean(np.log(countNdarrayTmp), axis=1))

    librarySizes = np.zeros(countNdarray.shape[1])

    for i in range(countNdarray.shape[1]):
        idx = countNdarray[:, i] > 0
        librarySizes[i] = np.median(countNdarray[idx, i] / geoMeans[idx])

    return librarySizes

if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
    else:
        print '*'*25
        FileIn = ld.LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_expt()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Genes.' % data.geneIDs.size
        print '*'*25
        data.libSizesRibo = lib_size(data.countRibo)
        data.libSizesRna  = lib_size(data.countRna)
        print 'Library size:'
        np.set_printoptions(precision=3)
        print data.experRibo
        print data.libSizesRibo
        print data.experRna
        print data.libSizesRna
