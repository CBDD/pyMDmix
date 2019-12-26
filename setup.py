#!/usr/bin/env python
"""Setuptools-based setup script for pyMDMix

REQUIREMENTS:
	- Working installation of NumPy
	- Working installation of Biskit

OPTIONAL:
	- pycuda for accelerated computing
	- matplotlib and scipy for visualization and running of clustering in HotSpotsManager functions

For a basic installation just type the command::
	python setup.py install

By default we use setuptools <http://pypi.python.org/pypi/setuptools>.
"""

import sys, os
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension
from setuptools.command.install import install

class CustomInstall(install):
    def createSolventDB(self, path):
		from pyMDMix.Solvents import SolventManager
    	import glob
		print "Building solvent database..."
		dbpath = os.path.join(path, 'SOLVENTS.db')
    	configfilelist = glob.glob(os.path.join(path, '*.config'))
    	maker = SolventManager()
    	for conffile in configfilelist:
			print 'Adding solvent from %s' % conffile
			maker.saveSolvent(maker.createSolvent(conffile), db = dbpath, createEmpty = True)
		if os.path.exists(dbpath): print "DONE creating solvent DB"

    def run(self):
		install.run(self)
		# Run solvent library building
		solventdbpath = os.path.join(self.install_lib, 'pyMDMix','data','solventlib')
		self.createSolventDB(solventdbpath)

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

scriptlist = ['src/mdmix']

def getRequirements():
	requirements=[]
	with open('requirements.txt', 'r') as reqfile:
		for line in reqfile.readlines():
			requirements.append(line.strip())
	return requirements

def getVersionFromInit():
	import pyMDMix
	return pyMDMix.__version__

setup(name='pyMDMix',
		cmdclass={'install': CustomInstall},
		zip_safe=False,
		version= getVersionFromInit(),
		description='Molecular Dynamics with organic solvent mixtures setup and analysis',
		author='Daniel Alvarez-Garcia',
		author_email='algarcia.daniel@gmail.com',
		url='',
		packages=['pyMDMix','pyMDMix.Actions'],
		include_package_data=True,
		package_data={'pyMDMix.data':['pyMDMix/data/*']},
		scripts=scriptlist,
		dependency_links=['https://sourcesup.renater.fr/frs/?group_id=180&release_id=2467#stable-releases-_2.8.1-title-content'],
		install_requires=getRequirements())
