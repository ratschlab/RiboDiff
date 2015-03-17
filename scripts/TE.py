#!/usr/bin/env python
"""
Master script running RiboDiff program.

Usage: 
    python TE.py -e <exp_outline.txt> -c <cnt_table.txt> -o <result.txt> -d [0] -r [1] -p [1] -q [0.1]

Details: 
    Check the MANUAL in RibiDiff directory. Type "python TE.py -h" to get the help information.
"""

import sys
import os
from ribodiff.parseopts import *

def main():
    """
    Master function definition to call different functionalities of RiboDiff.
    """

    opts = parse_options(sys.argv)

    import numpy as np
    import ribodiff.loadinput as ld
    import ribodiff.normalize as nm
    import ribodiff.estimatedisp as ed
    import ribodiff.testcnt as tc
    import ribodiff.writeres as wr
    import ribodiff.plot as pl

    print '*'*25
    FileIn = ld.LoadInputs(opts)
    data = FileIn.parse_expt()
    data = FileIn.read_count()
    print 'Read input files: Done.\n%i Gene(s) to be tested.' % data.geneIDs.size

    print '*'*25
    data.libSizesRibo = nm.lib_size(data.countRibo)
    data.libSizesRna  = nm.lib_size(data.countRna)
    print 'Library size:'
    np.set_printoptions(precision=3)
    print data.experRibo
    print data.libSizesRibo
    print data.experRna
    print data.libSizesRna

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
    args = sys.argv
    main()
