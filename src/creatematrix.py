import sys
import numpy as np

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
        explanatory[np.logical_and(experNdarray[:,1]=='Ribo', experNdarray[:,2]=='Treated'), 2] = 1
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

def create_matrix_2(experNdarray, model='H0'):

    width = len(experNdarray[:, 0])
    if model == 'H0':
        explanatory = np.zeros((width, 3))
        explanatory[experNdarray[:,2]=='Control', 0] = 1
        explanatory[experNdarray[:,2]=='Treated', 1] = 1
        explanatory[experNdarray[:,1]=='Rna', 2] = 1
    elif model == 'H1':
        explanatory = np.zeros((width, 4))
        explanatory[experNdarray[:,2]=='Control', 0] = 1
        explanatory[experNdarray[:,2]=='Treated', 1] = 1
        explanatory[experNdarray[:,1]=='Rna', 2] = 1
        explanatory[np.logical_and(experNdarray[:,1]=='Ribo', experNdarray[:,2]=='Treated'), 3] = 1
    else:
        sys.stderr.write('ERROR: The parameter \'model\' in create_sparse() can only accept \'H0\' or \'H1\' as its input.\n')

    return explanatory    
