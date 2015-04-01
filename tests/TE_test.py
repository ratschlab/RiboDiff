#!/usr/bin/env python 
"""
the script integrated with setuptools test command

to execute this script: 
    python setup.py test  
"""

import os 
import sys 
import unittest
import subprocess 

class TransEfficiency(unittest.TestCase):

    def testTransEfficiency(self):
        """
        """

        setup_file = sys.modules['__main__'].__file__
        setup_dir = os.path.abspath(os.path.dirname(setup_file))

        test_cmd = "python src/ribodiff/functional_test_te.py -e tests/experimental_design.csv -c tests/test_gene_count_table.tab -o tests/test_result.txt -d 0 -r 1 -p 1 -q 0.1"

        try:
            test_process = subprocess.Popen(test_cmd, shell=True)
            returncode = test_process.wait()

            if returncode != 0:
                raise Exception, "Exit status return code = %i" % returncode

            print "Test run finished, results are here %s/tests/" % setup_dir

        except Exception, e:
            print "Error running test cases.\n%s" % str(e) 
            sys.exit(0)


if __name__=="__main__":
    unittest.main()
