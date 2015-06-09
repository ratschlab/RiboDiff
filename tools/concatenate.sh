#!/usr/bin/env python

"""
This script concatenates simulated 'Ribo' and 'RNA' count files.
Usage: python concatenate.py <Ribo count file> <RNA count file> <output file>
"""

import sys
import numpy as np

def main():

    cntFileRibo = sys.argv[1]
    cntFileRna  = sys.argv[2]
    cntFileAll  = sys.argv[3]

    idsRibo = np.loadtxt(cntFileRibo, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))
    idsRna  = np.loadtxt(cntFileRna, dtype=str, delimiter='\t', skiprows=1, usecols=(0,))

    if not np.array_equal(idsRibo, idsRna):
        sys.stderr.write('\nError: ID in two input files are not the same!\n\n')
        sys.exit()
    else:
        ids = idsRibo.copy()

    with open(cntFileRibo, 'r') as FileInRibo:
        headerRibo = FileInRibo.readline().strip().split('\t')
        sampNumRibo = headerRibo.index('Dispersion')
    with open(cntFileRna, 'r') as FileInRna:
        headerRna = FileInRna.readline().strip().split('\t')
        sampNumRna = headerRna.index('Dispersion')

    cntRibo = np.loadtxt(cntFileRibo, dtype=str, delimiter='\t', skiprows=1, usecols=tuple(range(1, sampNumRibo)))
    cntRna  = np.loadtxt(cntFileRna, dtype=str, delimiter='\t', skiprows=1, usecols=tuple(range(1, sampNumRna)))

    ids = np.reshape(ids, (ids.size, 1))
    cntAll  = np.hstack([ids, cntRibo, cntRna])

    header = ['Entry']

    nameRibo = np.array(headerRibo[1 : sampNumRibo])
    idxRibo = np.hstack([np.arange(np.where(nameRibo == 'conditionA')[0].size) + 1, np.arange(np.where(nameRibo == 'conditionB')[0].size) + 1])
    for i in range(nameRibo.size):
        col  = 'Rb' + nameRibo[i] + str(idxRibo[i])
        header.append(col)

    nameRna  = np.array(headerRna[1 : sampNumRna])
    idxRna = np.hstack([np.arange(np.where(nameRna == 'conditionA')[0].size) + 1, np.arange(np.where(nameRna == 'conditionB')[0].size) + 1])
    for i in range(nameRna.size):
        col  = 'Rna' + nameRna[i] + str(idxRna[i])
        header.append(col)

    header = '\t'.join(header)

    np.savetxt(cntFileAll, cntAll, fmt='%s', delimiter='\t', header=header, comments='')

if __name__ == '__main__':

    main()
