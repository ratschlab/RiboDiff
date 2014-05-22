import sys
import loadinput as ld
import normalize as nm
import estimatedisp as ed
import numpy as np
import creatematrix as cm
import statsmodels.api as sm
from scipy.stats import chi2

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python testcount.py Experiment_Outline_File Gene_Count_File Output_File' + '\n')

def test_count(data):

    print 'Start the statistical test.'

    num = len(data.geneIDs)
    pval = np.empty((num, 1))
    pval.fill(np.nan)

    explanatory0 = cm.create_matrix_2(data.exper, model='H0')
    explanatory1 = cm.create_matrix_2(data.exper, model='H1')
    librarySizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    for i in range(num):

        if i % 50 == 0:
            print '%i genes finished...' % i
        if i+1 == num:
            print '%i genes finished...' % num

        response = np.hstack([data.countRibo[i, :], data.countRNA[i, :]])
        disper = data.disperFitted[i]
        if not np.isnan(disper):
            modNB0 = sm.GLM(response, explanatory0, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            modNB1 = sm.GLM(response, explanatory1, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            result0 = modNB0.fit()
            result1 = modNB1.fit()
            pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, explanatory1.shape[1] - explanatory0.shape[1])

    data.pval = pval

    return data

def write_output(data, fileName):

    geneIDs = data.geneIDs
    disperRaw = data.disperRaw.astype(str)
    disperFitted = data.disperFitted.astype(str)
    pval = data.pval.astype(str)
    outNdarray = np.hstack([geneIDs, disperRaw, disperFitted, pval])
    np.savetxt(fileName, outNdarray, fmt='%s', delimiter='\t')

if __name__ == '__main__':
    if len(sys.argv) < 4:
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
        print '*'*25
        data = ed.estimate_disp(data, method='simp')
        print 'Estimate dispersion: Done.'
        print '*'*25
        data = test_count(data)
        print 'Statistical test: Done.'
        print '*'*25
        write_output(data, sys.argv[3])
        print 'Write output file: Done.'
        print '*'*25
