import sys
import numpy as np
import loadinput as ld

def usage():
    sys.stderr.write('Usage:' + '\n' + 'python normalize.py Experiment_Outline_File Gene_Count_File' + '\n')

def lib_size(countNdarray):

    """ Calculate library size: median( gene[i] cnt / geometric mean[i] across all libraries ) """

    # add code to check if there are counts smaller than zero
    #

    countNdarrayTmp = countNdarray.copy()
    countNdarrayTmp = countNdarrayTmp + 1    #### better than countNdarrayTmp[countNdarrayTmp==0] = 1 ###
    geoMeans = np.exp(np.mean(np.log(countNdarrayTmp), axis=1))
    librarySizes = np.median(countNdarray / geoMeans[:, np.newaxis], axis=0)

    return librarySizes

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
    else:
        print '*'*25
        FileIn = ld.LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_exper()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
        print '*'*25
        data.libSizesRibo = lib_size(data.countRibo)
        data.libSizesRNA = lib_size(data.countRNA)
        print 'Library size:'
        np.set_printoptions(precision=3)
        print data.experRF
        print data.libSizesRibo
        print data.experRNA
        print data.libSizesRNA
