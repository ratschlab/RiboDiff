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
        explanatory[experNdarray[:,1]=='Ribosome_Footprint', 1] = 1
    elif model == 'H1':
        explanatory = np.zeros((width, 3))
        explanatory[:, 0] = 1
        explanatory[experNdarray[:,1]=='Ribosome_Footprint', 1] = 1
        explanatory[np.logical_and(experNdarray[:,1]=='Ribosome_Footprint' and experNdarray[:,2]=='Treated'), 1] = 1
    else:
        sys.stderr.write('ERROR: The parameter \'model\' in create_sparse() can only accept \'H0\' or \'H1\' as its input.\n')
    return explanatory

def create_matrix_1(experNdarray, model='H0'):



def estimate_disp(data):

    explanatory = create_matrix_0(data.exper)
    explanatory = sm.add_constant(explanatory, prepend=False)

    num = len(data.geneIDs)
#	for i in range(num):
    i = 1
    response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
    modNB = sm.GLM(response, explanatory, family=sm.families.NegativeBinomial())
    #modNB = sm.NegativeBinomial(response, explanatory, loglike_method='nb2')
    results = modNB.fit()
    print results.summary()
    print response
    print explanatory

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
        np.set_printoptions(precision=3)
        print data.experRF
        print data.libSizesRibo
        print data.experRNA
        print data.libSizesRNA
        data = estimate_disp(data)
