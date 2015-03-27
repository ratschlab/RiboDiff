SOFTWARE: RiboDiff
----------

###VERSION

0.1

###Authors

Yi Zhong, Philipp Drewe, Theofanis Karaletsos, Vipin Sreedharan and Gunnar Raetsch.

###DESCRIPTION

RiboDiff is a statistical framework to detect the protein translation 
efficiency change from ribosome footprint profiling data and RNA-Seq
data in two different experiment conditions.

###URL

http://bioweb.me/ribo

###REQUIREMENTS
* Python >= 2.6.6
* Numpy >= 1.8.0
* Scipy >= 0.13.3
* Statsmodels >= 0.6.0
* Matplotlib >= 1.3.0 

###INSTALLATION

For a global install, for which you need root permissions::

    python setup.py build
    python setup.py test [optional]
    python setup.py install 

For a local install::

    python setup.py build
    python setup.py test [optional]
    python setup.py install --prefix=$HOME

###CONTENTS

All relevant scripts for RiboDiff are located in the subdirectory src. 

* src   - main codebase for RiboDiff;
* tests  - dataset and script for functional test;
* tools - util functions for simulating negative binomial count data.
* scripts - TE.py - the main script to start RiboDiff.

###GALAXY

https://galaxy.cbio.mskcc.org

###DOCUMENTATION

To use RiboDiff, please refer to the instructions in MANUAL in this directory.

###LICENSE

RiboDiff is licensed under the GPL version 3 or any later version
(cf. LICENSE).

###CITE US

If you use RiboDiff in your research you are kindly asked to cite the
following publication:

XXX
