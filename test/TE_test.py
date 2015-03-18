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
sub_module_dir = os.path.join(setup_dir, "src")
print sub_module_dir

test_cmd = "python %s/scripts/TE.py -e %s/test/exp_outline.txt -c %s/test/test_gene_count_table.txt -o %s/test/test_result.txt -d 0 -r 1 -p 1 -q 0.1" % (setup_dir, setup_dir, setup_dir, setup_dir)

try:
    os.chdir(sub_module_dir)
    test_process = subprocess.Popen(test_cmd, shell=True)
    returncode = test_process.wait()

    if returncode != 0:
        raise Exception, "Exit status return code = %i" % returncode

    print "Test run finished. Please check the results in the %s/test/ directory." % setup_dir

except Exception, e:
    print "Error running test cases.\n%s" % str(e) 
    sys.exit(0)
