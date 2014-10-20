#!/usr/bin/env python

import sys
import numpy as np
import loadinput as ld

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python normalize.py Experiment_Outline_File Count_File' + '\n')

def lib_size(countNdarray):

    """ Calculate library size: median( gene[i] cnt / geometric mean[i] across all libraries ) """

    # add code to check if there are counts smaller than zero
    #

    countNdarrayTmp = countNdarray.copy()
    countNdarrayTmp = countNdarrayTmp + 1
    geoMeans = np.exp(np.mean(np.log(countNdarrayTmp), axis=1))

    librarySizes = np.zeros(countNdarray.shape[1])

    for i in range(countNdarray.shape[1]):
        idx = countNdarray[:, i] > 0
        librarySizes[i] = np.median(countNdarray[idx, i] / geoMeans[idx])

    ### Plot histogram ############################
    #matrix = countNdarray / geoMeans[:, np.newaxis]
    #plt.hist(matrix[:,0], 40, range=(0, 2), color='blue')
    #plt.savefig('./libsizeHistgram.pdf')
    ###############################################

    return librarySizes

if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
    else:
        print '*'*25
        FileIn = ld.LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_expt()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
        print '*'*25
        data.libSizesRibo = lib_size(data.countRibo)
        data.libSizesRNA  = lib_size(data.countRNA)
        print 'Library size:'
        np.set_printoptions(precision=3)
        print data.experRibo
        print data.libSizesRibo
        print data.experRNA
        print data.libSizesRNA
