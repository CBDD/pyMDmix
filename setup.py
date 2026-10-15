#!/usr/bin/env python

import sys, os
from setuptools import setup, Extension
from setuptools.command.install import install

# Make sure I have the right Python version.
if sys.version_info[:2] < (3, 9):
    print("pyMDMix requires Python 3.9 or later. Python %d.%d detected" % \
        sys.version_info[:2])
    print("Please upgrade your version of Python.")
    sys.exit(1)

# Make sure AMBERHOME environ variable is set
if not os.environ.get('AMBERHOME'):
    print("AMBERHOME env variable not set! Please set this variable pointing to AMBER package installation directory.")
    sys.exit(1)

script_list = ['src/mdmix']

def getRequirements():
    requirements=[]
    with open('requirements.txt', 'r') as reqfile:
        for line in reqfile.readlines():
            requirements.append(line.strip())
    return requirements

def getVersionFromInit():
    # Single source of truth: __version__ in pyMDMix/__init__.py
    import re
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pyMDMix', '__init__.py')) as f:
        return re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', f.read()).group(1)

setup(
    name='pyMDMix',
    zip_safe=False,
    version= getVersionFromInit(),
    description='Molecular Dynamics with organic solvent mixtures setup and analysis',
    author='Daniel Alvarez-Garcia',
    author_email='algarcia.daniel@gmail.com',
    url='',
    packages=['pyMDMix','pyMDMix.Actions', 'pyMDMix.Commands'],
    include_package_data=True,
    package_data={'pyMDMix.data':['pyMDMix/data/*']},
    scripts=script_list,
    python_requires='>=3.9',
    classifiers=['Programming Language :: Python :: 3', 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'],
    install_requires=getRequirements()
)
