#!/usr/bin/env python 
"""
test run script 

python setup.py test  
"""

import os 
import sys 
import unittest
import subprocess 

class TransEfficiency(unittest.TestCase):

    def testTransEfficiency(self):
        """
        test function 
        """

        setup_file = sys.modules['__main__'].__file__
        setup_dir = os.path.abspath(os.path.dirname(setup_file))

        test_cmd = "python src/ribodiff/functional_test_te.py -e tests/exp_outline.txt -c tests/test_gene_count_table.txt -o tests/test_result.txt -d 0 -r 1 -p 1 -q 0.1"

        try:
            test_process = subprocess.Popen(test_cmd, shell=True)
            returncode = test_process.wait()

            if returncode != 0:
                raise Exception, "Exit status return code = %i" % returncode

            print "Test run finished. Please check the results in the %s/tests/ directory." % setup_dir

        except Exception, e:
            print "Error running test cases.\n%s" % str(e) 
            sys.exit(0)


if __name__=="__main__":
    unittest.main()
