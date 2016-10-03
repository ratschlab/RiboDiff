#!/usr/bin/env python
"""
Installation scripts based on these:
http://www.scotttorborg.com/python-packaging/minimal.html
http://packages.python.org/distribute/setuptools.html
"""

from setuptools import setup

def readme():
    with open('README.md') as f:
       return f.read()


setup(name='RiboDiff',
      version='0.2.1',
      description='Detect protein translational efficiency difference using deep sequencing data.',
      long_description=readme(),
      keywords='ribo-seq rna-seq differential-expression-testing translational-efficiency',
      url='http://galaxy.cbio.mskcc.org',
      author=['Yi Zhong', 'Vipin T. Sreedharan'],
      author_email=['zhongy@cbio.mskcc.org', 'vipin@cbio.mskcc.org'], 
      license='GPLv3',
      packages=['ribodiff'],
      package_dir={'ribodiff': 'src/ribodiff'},
      install_requires=['numpy', 'scipy', 'statsmodels', 'matplotlib'],
      test_suite="tests",
      zip_safe=False)
