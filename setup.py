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

By default we use setuptools <http://pypi.python.org/pypi/setuptools>.  The
details of such an "EasyInstall" installation procedure are shown on

  http://peak.telecommunity.com/DevCenter/EasyInstall

By changing the code below you can also switch to a standard distutils
installation.
"""
#------------------------------------------------------------
# selection of the installation system
#------------------------------------------------------------
#
# Standard distutils-based installation:
#
##from distutils.core import setup, Extension

# setuptools ("EasyInstall") installation:
#
# If you want EasyInstall features then enable the next three lines and comment
# out the preceding line 'from distutils.core import ...'
#
import sys, os
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension
from setuptools.command.install import install
import numpy

class CustomInstall(install):
    def createSolventDB(self, path):
	from pyMDMix.Solvents import SolventManager
    	import glob
	print "Building solvent database..."
	dbpath = os.path.join(path, 'SOLVENTS.db')
    	configfilelist = glob.glob(os.path.join(path, '*.config'))
    	maker=SolventManager()
    	for conffile in configfilelist:
        	print 'Adding solvent from %s'%conffile
        	maker.saveSolvent(maker.createSolvent(conffile), db=dbpath, createEmpty=True)
	if os.path.exists(dbpath): print "DONE creating solvent DB"

    def run(self):
        install.run(self)
	# Run solvent library building
	solventdbpath = os.path.join(self.install_lib, 'pyMDMix','data','solventlib')
	self.createSolventDB(solventdbpath)
	#import glob
	#import subprocess
	##for f in glob.glob(solventdbpath+os.sep+'*.config'):
	##	print "Adding solvent in %s file"%f
	#p = subprocess.Popen(['python',solventdbpath+os.sep+'createSolventDB.py'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#print p.stdout.read()
	#print p.stderr.read()
	#print solventdbpath+os.sep+"SOLVENT.db"

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 5):
    print "pyMDMix requires Python 2.5 or later. Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of Python."
    sys.exit(-1)

# Make sure AMBERHOME environ variable is set
if not os.environ.get('AMBERHOME'):
    print "AMBERHOME env variable not set! Please set this variable pointing to AMBER package installation directory."

#scriptlist = ['scripts/prepareMDMixProject.py','scripts/runCenteringAndRawEnergyCalculations.py','scripts/createMinotauroQueueInput.py','scripts/createPaintersQueueInput.py', 'scripts/runReplicaCentering.py','scripts/extendReplica.py', 'scripts/printReplicaInfo.py','scripts/createReplicaCenteringInput.py','scripts/runReplicaDensityAndRawCalculation.py', 'scripts/mdmix']
scriptlist = ['src/mdmix']

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = [numpy_include]


setup(name='pyMDMix',
      cmdclass={'install': CustomInstall},
      zip_safe=False,
      version='0.1.3',
      description='Molecular Dynamics with organic solvent mixtures setup and analysis',
      author='Daniel Alvarez-Garcia',
      author_email='algarcia.daniel@gmail.com',
      url='',
      packages=['pyMDMix','pyMDMix.Actions'],# 'MDMix.KDTree'],
      include_package_data=True,
      package_data={'pyMDMix.data':['pyMDMix/data/*']},
      scripts=scriptlist,
	dependency_links=['https://sourcesup.renater.fr/frs/?group_id=180&release_id=2467#stable-releases-_2.8.1-title-content'],
      install_requires=["numpy","ScientificPython","Biskit","scipy","matplotlib","mechanize"])

      #extras_require={'HotSpotsManager':["scipy","matplotlib"],
#			"CUDA":["pycuda"]})
