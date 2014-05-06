import sys
import loadinput as ld
import normalize as nm
import creatematrix as cm
import dispersion as dp

def usage():

    sys.stderr.write('Usage:' + '\n' + 'python ReadInputs.py Experiment_Outline_File Gene_Count_File' + '\n')

def estimate_disp(data, method='simp'):

    #explanatory = cm.create_matrix_0(data.exper, model='H1')
    #explanatory = cm.create_matrix_1(data.exper, model='H1')
    explanatory = cm.create_matrix_2(data.exper, model='H1')
    #explanatory = sm.add_constant(explanatory, prepend=False)
    data.matrix = explanatory

    if method == 'simp':
        data = dp.disper_simp(data)
        data = dp.fit_disper_gamma(data)

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
        data = estimate_disp(data, method='simp')
        print 'Estimate dispersion: Done.'
        print '*'*25
