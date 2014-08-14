#!/usr/bin/env python
##
## pyMDMix --- http://mdmix.sourceforge.net
## Software for preparation, analysis and quality control
## of solvent mixtures molecular dynamics
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
## Please cite your use of pyMDMix in published work:
##
##    TOBEPUBLISHED
##
"""
Module adapted from Biskit.test
===============================
Inherited functions are modified to make testing work on pyMDMix tree structure.
All rights to Raik.
"""

from Biskit.test import *
import pyMDMix.tools as T
import os.path as osp

# TEST CASE FOR FUTURE
FUTURE = 6

class PyMDMixTestLoader( BiskitTestLoader ):
    """
    Wrap for BiskitTestLoader to match pyMDMix structure
    """
    def modulesFromPath( self, path=osp.split(T.projectRoot())[0], module='pyMDMix' ):
	import glob
	module_folder = module.replace('.', os.path.sep)
        files = glob.glob( os.path.join( path, module_folder,'*.py' ) )
        files = map( T.stripFilename, files )
        files = [ f for f in files if f[0] != '_' ]
        r = []
        for f in files:
            try:
                r += [ __import__( '.'.join([module, f]), globals(),
                                   None, [module]) ]
            except:
                pass  ## temporary // remove after testing
        return r
    def collectTests( self, path=osp.split(T.projectRoot())[0], module='pyMDMix' ):
	modules = self.modulesFromPath( path=path, module=module )
        self.addTestsFromModules( modules )


############################################
### Script functions adapted from Biskit ###
############################################

def _use( defaults ):
    print """
Run unittest tests for pyMDMix.

    test.py [-i |include tag1 tag2..| -e |exclude tag1 tag2..|
             -p |package1 package2..|
             -v |verbosity| -log |log-file| -nox ]

    i    - include tags, only run tests with at least one of these tags   [All]
    e    - exclude tags, do not run tests labeled with one of these tags  [old]
         valid tags are:
             long   - long running test case
             pvm    - depends on PVM
             exe    - depends on external application
             old    - is obsolete
         (If no tags are given to -i this means all tests are included)

    p    - packages to test, e.g. pyMDMix.Actions pyMDMix.Parsers         [All]
    v    - int, verbosity level, 3 switches on several graphical plots      [2]
    log  - path to logfile (overriden); empty -log means STDOUT        [STDOUT]
    nox  - suppress test plots                                          [False]
    dry  - do not actually run the test but just collect tests          [False]

Examples:

    * Run all but long or obsolete tests from pyMDMix and pyMDMix.Actions:
    test.py -e old long  -p pyMDMix pyMDMix.Actions


Default options:
"""
    for key, value in defaults.items():
        print "\t-",key, "\t",value

    sys.exit(0)


def _str2tags( s ):
    """convert list of string options to list of valid TAGS"""
    try:
        r = [ x.upper() for x in s if x ] ## to list of uppercase str
        r = [ eval( x ) for x in r ]      ## to list of int
    except:
        EHandler.error('unrecognized tags: %r'%s)

    return r

def _convertOptions( o ):
    o['i'] = _str2tags( T.toList( o['i'] ) )
    o['e'] = _str2tags( T.toList( o['e'] ) )
    o['v'] = int( o['v'] )
    o['nox'] = ('nox' in o)
    o['dry'] = ('dry' in o)
    o['debug'] = ('debug' in o)
    if o['log']:
        o['log'] = LogFile( o['log'] )
    else:
        o['log'] = StdLog()
    o['p'] = T.toList( o['p'] )


if __name__ == '__main__':
    from Biskit import EHandler
    import sys


    defaults = {'i':'',
                'e':'future',
                'p':['pyMDMix', 'pyMDMix.Actions'],
                'v':'2',
                'log': '', ##T.testRoot()+'/test.log',
                }

    o = T.cmdDict( defaults )

    if len( sys.argv ) == 1 and 'test.py' in sys.argv[0]:
        _use( defaults )
        
    _convertOptions( o )

    BiskitTest.VERBOSITY = o['v']
    BiskitTest.DEBUG = o['debug']
    
    l = PyMDMixTestLoader( allowed=o['i'], forbidden=o['e'],
                          verbosity=o['v'], log=o['log'], debug=o['debug'])


    for package in o['p']:
        print 'collecting ', repr( package )
        l.collectTests( module=package )

    l.run( dry=o['dry'] )
    l.report()

    print "DONE"
