import sys
import numpy as np
import loadinput as ld
import normalize as nm
import statsmodels.api as sm

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

def create_matrix_0(experNdarray, model='H0'):

    width = len(experNdarray[:, 0])
    if model == 'H0':
        explanatory = np.zeros((width, 2))
        explanatory[:, 0] = 1
        explanatory[experNdarray[:,1]=='Ribo', 1] = 1
    elif model == 'H1':
        explanatory = np.zeros((width, 3))
        explanatory[:, 0] = 1
        explanatory[experNdarray[:,1]=='Ribo', 1] = 1
        explanatory[np.logical_and(experNdarray[:,1]=='Ribo', experNdarray[:,2]=='Treated'), 1] = 1
    else:
        sys.stderr.write('ERROR: The parameter \'model\' in create_sparse() can only accept \'H0\' or \'H1\' as its input.\n')
    return explanatory

def create_matrix_1(experNdarray, model='H0'):

    width = len(experNdarray[:, 0])
    idxRF = np.nonzero(experNdarray[:,1] == 'Ribo')[0]
    col_num = np.sum(experNdarray[:,1] == 'Ribo')
    explanatory = np.zeros((width, col_num))
    for i in range(col_num):
        idx = np.logical_or(np.arange(width) == idxRF[i], np.logical_and(experNdarray[:, 1] == 'Rna', experNdarray[:, 2] == experNdarray[idxRF[i], 2]))
        explanatory[idx, i] = 1

    is_Rna = np.zeros((width, 1))
    is_Rna[experNdarray[:, 1] == 'Rna'] = 1
    explanatory = np.hstack((is_Rna, explanatory))

    if model == 'H1':
        interact = np.zeros((width, 1))
        interact[np.logical_and(experNdarray[:,1]=='Ribo', experNdarray[:,2]=='Treated')] = 1
        explanatory = np.hstack((explanatory, interact))

    return explanatory

def estimate_disp(data):

    #explanatory = create_matrix_0(data.exper, model='H1')
    explanatory = create_matrix_1(data.exper, model='H1')
    #explanatory = sm.add_constant(explanatory, prepend=False)

    num = len(data.geneIDs)
#   for i in range(num):
    i = 0
    response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])
    disper = 0.1
    modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
    #modNB = sm.NegativeBinomial(response, explanatory, loglike_method='nb2')
    results = modNB.fit()
    print results.summary()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    else:
        print '*'*25
        FileIn = ld.LoadInputs(sys.argv[1], sys.argv[2])
        data = FileIn.parse_exper()
        data = FileIn.read_count()
        print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size
        print '*'*25
        data.libSizesRibo = nm.lib_size(data.countRibo)
        data.libSizesRNA = nm.lib_size(data.countRNA)
        print 'Library size:'
        #np.set_printoptions(precision=3)
        print data.experRF
        print data.libSizesRibo
        print data.experRNA
        print data.libSizesRNA
        print '*'*25
        data = estimate_disp(data)
