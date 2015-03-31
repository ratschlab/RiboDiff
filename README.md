SOFTWARE: RiboDiff
----------

###VERSION

0.1

###Authors

Yi Zhong, Theofanis Karaletsos, Philipp Drewe, Vipin Sreedharan and Gunnar Raetsch.

###DESCRIPTION

RiboDiff is a statistical tool to detect the protein translation 
efficiency change from ribosome footprint profiling data and RNA-Seq
data under two different experiment conditions.

###URL

http://bioweb.me/ribo

###REQUIREMENTS
* Python2 >= 2.6.6
* Numpy >= 1.8.0
* Scipy >= 0.13.3
* Matplotlib >= 1.3.0 
* Statsmodels >= 0.5.0

These requirements can either be installed individually or as a Python distribution
that includes all the required packages. Please find more details at
http://www.scipy.org/install.html

###INSTALLATION

For a global install (you need root permissions to do so):

    python setup.py build
    python setup.py test (optional step)
    python setup.py install

For a local install:

    python setup.py build
    python setup.py test (optional step)
    python setup.py install --prefix=/PATH/PREFIX

Please replace /PATH/PREFIX with your prefered path where RiboDiff will be installed. 
Make sure this path exists in your PATHONPATH environment variable. If not, you 
can edit PATHONPATH in your ~/.bashrc if you are using Linux or ~/.bash_profile if 
Mac OS X by adding the folowing line:

	export PYTHONPATH=/PATH/PREFIX/lib/python2.7/site-packages/:$PYTHONPATH

###CONTENTS

All relevant scripts for RiboDiff are located in the subdirectory src. 

* src   - main codebase for RiboDiff;
* tests  - dataset and script for functional test;
* test-data  - test dataset for Galaxy system;
        (move to your galaxy_root_folder/test-data/)
        You may need to move the test files into your test-data directory so galaxy can find them. 
        If you want to run the functional tests example as: 

        sh run_functional_tests.sh -id ribodiff

* tools - util functions for simulating negative binomial count data.
* scripts - TE.py - the main script to start RiboDiff.

###GALAXY

https://galaxy.cbio.mskcc.org/tool_runner?tool_id=ribodiff

###DOCUMENTATION

To use RiboDiff, please refer to the instructions in MANUAL in this directory.

###LICENSE

RiboDiff is licensed under the GPL version 3 or any later version
(cf. LICENSE).

###CITE US

If you use RiboDiff in your research you are kindly asked to cite the
following publication:

Zhong Y, Karaletsos T, Drewe P, et al. RiboDiff: Detecting Changes of Translation 
Efficiency from Ribosome Footprints. bioRxiv. doi: 10.1101/017111.
