#!/usr/bin/env python

import sys
from parseopts import *
import pdb

def main():

    opts = parse_options(sys.argv)

    import os
    import numpy as np
    import loadinput as ld
    import normalize as nm
    import estimatedisp as ed
    import testcnt as tc
    import writeres as wr
    import plot as pl

    print '*'*25
    FileIn = ld.LoadInputs(opts)
    data = FileIn.parse_expt()
    data = FileIn.read_count()
    print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size

    print '*'*25
    data.libSizesRibo = nm.lib_size(data.countRibo)
    data.libSizesRNA  = nm.lib_size(data.countRNA)
    print 'Library size:'
    np.set_printoptions(precision=3)
    print data.experRibo
    print data.libSizesRibo
    print data.experRNA
    print data.libSizesRNA

    print '*'*25
    data = ed.estimate_disp(data, opts)
    print 'Estimate dispersion: Done.'

    print '*'*25
    data = tc.test_count(data, opts)
    data = tc.adj_pval(data, opts)
    print 'Statistical test: Done.'

    print '*'*25
    data = tc.cal_TEchange(data)
    print 'Calculate TE and fold change: Done.'

    print '*'*25
    wr.write_result(data, opts)
    print 'Write output file: Done.'

    print '*'*25
    wr.save_data(data, opts)
    try:
        os.remove(opts.resPath + 'TmpData.pkl')
    except OSError:
        pass
    print 'Save data: Done.'

    if opts.plots:
        print '*'*25
        pl.make_plots(data, opts)
        print 'Make plots: Done.'

    print '*'*25

if __name__ == '__main__':

    main()
