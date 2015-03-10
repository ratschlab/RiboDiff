#!/usr/bin/env python 
"""
"""

import os 
import sys 
import subprocess 

setup_file = sys.modules['__main__'].__file__
print setup_file
setup_dir = os.path.abspath(os.path.dirname(setup_file))
print setup_dir

test_cmd = "python src/ribodiff/TE.py -e test/exp_outline.txt -c test/test_gene_count_table.txt -o test/test_result.txt -d 0 -r 1 -p 1 -q 0.1"

try:
    os.chdir(setup_dir)
    test_process = subprocess.Popen(test_cmd, shell=True)
    test_process.wait()
except Exception, e:
    print "Error running test cases.\n%s" % str(e) 
    sys.exit(0)

print "Test run finished. Please check the results in the ./test/ directory."
