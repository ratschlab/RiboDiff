import sys
import loadinput as ld
import normalize as nm
import estimatedisp as ed
import numpy as np
import creatematrix as cm
import statsmodels.api as sm
from scipy.stats import chi2
import statsmodels.sandbox as sms
import cPickle as pickle

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
        disper = data.disperAdj[i]
        if not np.isnan(disper):
            modNB0 = sm.GLM(response, explanatory0, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            modNB1 = sm.GLM(response, explanatory1, family=sm.families.NegativeBinomial(alpha=disper), offset=np.log(librarySizes))
            result0 = modNB0.fit()
            result1 = modNB1.fit()
            pval[i] = 1 - chi2.cdf(result0.deviance - result1.deviance, explanatory1.shape[1] - explanatory0.shape[1])

    data.pval = pval

    return data

def adj_pval(data, mthd='BH'):

    pval = data.pval.copy()
    idx = ~np.isnan(pval)

    if mthd == 'BH':
        method = 'fdr_bh'
    elif mthd == 'Bonferroni':
        method = 'bonferroni'
    elif mthd == 'Holm':
        method = 'holm'
    elif mthd == 'Hochberg':
        method = 'simes-hochberg'
    elif mthd == 'Hommel':
        method = 'hommel'
    elif mthd == 'BY':
        method = 'fdr_by'
    elif mthd == 'TSBH':
        method = 'tsbh'
    else:
        sys.stderr.write('ERROR: The methods for multiple test correction can only accept \'Bonferroni\', \'Holm\', \'Hochberg\', \'Hommel\', \'BH\', \'BY\' or \'TSBH\' as its input.\n')

    mtc = sms.stats.multicomp.multipletests(pval[idx], alpha=0.05, method=method, returnsorted=False)

    padj = pval.copy()
    padj[idx] = mtc[1]
    data.padj = padj

    return data

def cal_TEchange(data):

    const = 1e-5

    count = np.hstack([data.countRibo, data.countRNA])
    libSizes = np.hstack([data.libSizesRibo, data.libSizesRNA])

    idxRFCtl  = np.intersect1d(data.idxRF,  data.idxCtl) - 1
    idxRNACtl = np.intersect1d(data.idxRNA, data.idxCtl) - 1
    TEctl = (const + np.mean((count/libSizes)[:, idxRFCtl], axis=1)) / (const + np.mean((count/libSizes)[:, idxRNACtl], axis=1))

    idxRFTrt  = np.intersect1d(data.idxRF,  data.idxTrt) - 1
    idxRNATrt = np.intersect1d(data.idxRNA, data.idxTrt) - 1
    TEtrt = (const + np.mean((count/libSizes)[:, idxRFTrt], axis=1)) / (const + np.mean((count/libSizes)[:, idxRNATrt], axis=1))

    logFoldChange = np.log2(TEtrt) - np.log2(TEctl)
    logFoldChange = np.reshape(logFoldChange, (len(logFoldChange), 1))

    data.logFoldChange = logFoldChange

    return data

def write_output(data, fileName):

    geneIDs = data.geneIDs
    disperRaw = data.disperRaw.astype(str)
    disperFitted = data.disperFitted.astype(str)
    disperAdj = data.disperAdj.astype(str)
    pval = data.pval.astype(str)
    padj = data.padj.astype(str)
    logFoldChange = data.logFoldChange
    outNdarray = np.hstack([geneIDs, disperRaw, disperFitted, disperAdj, pval, padj, logFoldChange])
    header = "geneIDs\tdisperRaw\tdisperFitted\tdisperAdj\tpval\tpadj\tlogFoldChange"
    np.savetxt(fileName, outNdarray, fmt='%s', delimiter='\t', header=header, comments='')

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
        data = ed.estimate_disp(data, mthd='simp')
        print 'Estimate dispersion: Done.'
        print '*'*25
        data = test_count(data)
        data = adj_pval(data, mthd='BH')
        print 'Statistical test: Done.'
        print '*'*25
        data = cal_TEchange(data)
        print 'Calculate TE fold change: Done.'
        print '*'*25
        write_output(data, sys.argv[3])
        print 'Write output file: Done.'
        print '*'*25
        pos = sys.argv[3].rfind('/') + 1
        pathOutput = sys.argv[3][:pos]
        pklFile = pathOutput + 'data.pkl'
        with open(pklFile, 'wb') as FileOut:
            pickle.dump(data, FileOut)
        print 'Save data: Done.'
        print '*'*25