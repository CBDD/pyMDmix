#!/usr/bin/env python

import sys, os
from setuptools import setup, Extension
from setuptools.command.install import install

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 5):
    print "pyMDMix requires Python 2.5 or later. Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of Python."
    sys.exit(1)

# Make sure AMBERHOME environ variable is set
if not os.environ.get('AMBERHOME'):
    print "AMBERHOME env variable not set! Please set this variable pointing to AMBER package installation directory."
    sys.exit(1)

script_list = ['src/mdmix']

def getRequirements():
    requirements=[]
    with open('requirements.txt', 'r') as reqfile:
        for line in reqfile.readlines():
            requirements.append(line.strip())
    return requirements

def getVersionFromInit():
    return '0.2.6'

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
    dependency_links=['https://sourcesup.renater.fr/frs/?group_id=180&release_id=2467#stable-releases-_2.8.1-title-content'],
    install_requires=getRequirements()
)
