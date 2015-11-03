#!/usr/bin/env python
"""
Calculating the sequencing library size.
""" 

import sys
import numpy as np

def usage():
    """
    Help message.
    """
    sys.stderr.write('Usage:' + '\n' + 'python normalize.py <Count File>' + '\n')
    sys.stderr.write('Note: please calculate library size for RNA-Seq and ribosome footprint separately!' + '\n')

def lib_size(countNdarray):
    """ Calculating library size.
    
    @args countNdarray: read count data
    @type countNdarray: numpy array
    """

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

    if len(sys.argv) != 2:
        usage()
    else:
        print '*'*25

        with open(sys.argv[1], 'r') as FileIn:
            header = np.array(FileIn.readline().strip().split('\t'))

        count = np.loadtxt(sys.argv[1], dtype=int, delimiter='\t', skiprows=1, usecols=range(1, header.size))

        print 'Read input files: Done.\n%i Genes.' % count[:, 0].size
        print '*'*25
        libSizes = lib_size(count)
        print 'Library size:'
        np.set_printoptions(precision=3)
        print header[1:]
        print libSizes
