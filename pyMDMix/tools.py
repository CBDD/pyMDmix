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
Extension of Biskit.tools module
================================

Incorporates all functions of Biskit.tools
Overwrites some of them and adds new ones

"""

__author__="dalvarez"
__date__ ="$16-ene-2014 17:49:02$"
import sys, os
import os.path as osp
import logging

# Disable Biskit Loading warning
DISABLE_BISKIT_LOADING_WARNS = True
if DISABLE_BISKIT_LOADING_WARNS:
    olderr = sys.stderr
    sys.stderr = open(os.devnull, 'wb')
    import Biskit.tools
    sys.stderr.close()
    sys.stderr = olderr

from Biskit.tools import *

class InvalidPath( ToolsError ):
    pass

class InvalidBinary( ToolsError ):
    pass


def projectRoot(file=None):
    """
    Root of pyMDMix project.

    @return: absolute path of the root of current project::
    @rtype: string
    """
    ## import this module
    import pyMDMix
    ## get location of this module
    f = absfile(pyMDMix.__file__)
    f = osp.dirname(f)

    if file: f = osp.join(f, file)

    return f

def dataRoot(subfolder=''):
    """
    Root of pyMDMix data directory.

    @param subfolder: str, optional sub-folder of data folder

    @return: absolute path
    @rtype: string
    """
    return osp.join(projectRoot(), 'data', subfolder)

def testRoot(*subpath):
    """
    Root of pyMDMix test directory.

    @return: absolute path
    @rtype: string
    """
    f = dataRoot('test')
    f = osp.join(f, *subpath)
    return f

def templatesRoot(file=None):
    """
    Root of pyMDMix test directory.

    @return: absolute path
    @rtype: string
    """
    f = dataRoot('templates')
    if file: f = osp.join(f, file)
    return f

def solventsRoot(file=None):
    """
    Root of pyMDMix test directory.

    @return: absolute path
    @rtype: string
    """
    f = dataRoot('solventlib')
    if file: f = osp.join(f, file)
    return f

def simplifyNestedList(inlist, l=[]):
    "Return a plain list from nested list"
    for el in inlist:
        if isinstance(el, list):
            simplifyNestedList(el, l)
#            [l.append(i) for i in el]
        else:
            l.append(el)
    return l

def filePermission(path):
    dic = {}
    p = path
    if not os.access(p, os.F_OK):
        dic["NOEXISTS"] = True
    else:
        dic["NOEXISTS"] = False
    if(os.access(p,os.R_OK)):
        dic["READ"] = True
    else:
        dic["READ"] = False
    if(os.access(p, os.W_OK)):
        dic["WRITE"] = True
    else:
        dic["WRITE"] = False
    if (os.access(p, os.X_OK)):
        dic["EXECUTE"] = True
    else:
        dic["EXECUTE"] = False

    return dic


def validPath(v):
    """
    @param v: potential path name
    @type  v: str

    @return: validated absolute Path
    @rtype : str

    @raise InvalidPath: if path is not found
    """
    try:
        v = absfile( v )
        if not v or not os.path.exists( v ):
            raise InvalidPath, 'invalid path %r' % v

        return v

    except InvalidPath, e:
        raise
    except Exception, e:
        raise InvalidPath, 'error during path validation: %r' % str(e)


def validBinary(v):
    """
    @param v: potential binary path
    @type  v: str

    @return: validated absolute path to existing binary
    @rtype : str

    @raise InvalidBinary: if path is not found
    """
    try:
        if not v:
            raise IOError, 'empty path'

        return absbinary( v )

    except IOError, msg:
        raise InvalidBinary( str(msg) )

def amberMaskToDict(maskstring):
    "Split string of format RESNAME@ATOMNAME,ATOMNAME;RESNAME@ATOMNAMES and return a dictionary {RES:[ATOMNAME, ATOMNAME..], RES:[ATOMNAMES]}"
    if ':' == maskstring[0]: maskstring = maskstring[1:]
    splited = maskstring.split(';')
    out = {}
    for el in splited:
        if '@' in el:
            split1 = el.split('@')
            if not out.has_key(split1[0]):
                out[split1[0]] =[k.strip() for k in split1[1].split(',')]
            else:
                out[split1[0]] += [k.strip() for k in split1[1].split(',')]
        else:
            out[el] = ['all'] # All atoms if just residue is given
    return out


def parseNumMask(mask):
    "Transform mask of form 1,2,5:10 to a list of ints"
    import numpy as npy
    s = mask.split(',')
    out = []
    for el in s:
        if ':' in el:
            a, b = map(int, el.split(':'))
            out += npy.arange(a,b+1).tolist()
        else:
            out += [int(el)]
    out.sort()
    return out

def numListToMask(numlist):
    "Transform a list of integers into mask eg:1-5,7,8-20"
    if len(numlist) <= 1: return numlist
    from operator import itemgetter
    from itertools import groupby
    ranges = []
    numlist.sort()
    for key, group in groupby(enumerate(numlist), lambda (index, item): index - item):
        group = map(itemgetter(1), group)
        if len(group) > 1:
            ranges.append('{0}-{1}'.format(group[0], group[-1]))
        else:
            ranges.append('{}'.format(group[0]))
    return ','.join(ranges)

##### CONSTRUCT A BROWSER
from Browser import Browser
BROWSER = Browser()

##### CONSTRUCT EXECUTOR
from Executor import Executor
EXECUTOR = Executor()

##### LOGGING TOOLS
class LogFormatter(logging.Formatter):
    "Custom formatter for logging INFO records in different format to warnings and errors"
    def format(self, record):
        #compute s according to record.levelno
        #for example, by setting self._fmt
        #according to the levelno, then calling
        #the superclass to do the actual formatting
        if record.levelno == logging.INFO:
            self._fmt = "%(levelname)s\t%(message)s"
        else:
            self._fmt = "###\t%(levelname)s (%(name)s)\t%(message)s"
        return super(LogFormatter, self).format(record)

import traceback,sys

def traceback_plus():
    """
    Print the usual traceback information, followed by a listing of all the
    local variables in each frame.
    """
    tb = sys.exc_info()[2]
    while 1:
        if not tb.tb_next:
            break
        tb = tb.tb_next
    stack = []
    f = tb.tb_frame
    while f:
        stack.append(f)
        f = f.f_back
    stack.reverse()
    traceback.print_exc()
    print "Locals by frame, innermost last"
    for frame in stack:
        print
        print "Frame %s in %s at line %s" % (frame.f_code.co_name,
                                             frame.f_code.co_filename,
                                             frame.f_lineno)
        for key, value in frame.f_locals.items():
            print "\t%20s = " % key,
            #We have to be careful not to cause a new error in our error
            #printer! Calling str() on an unknown object could cause an
            #error we don't want.
            try:                   
                print value
            except:
                print "<ERROR WHILE PRINTING VALUE>"

#### CLEAN NAMESPACE
del Browser, Executor


if __name__ == "__main__":
    print projectRoot()
    print solventsRoot()
    l = [[1,2,3],4,[[5,6],[[7],[8],[9,10,11]],[[12,13],14]]]
    print simplifyNestedList(l)

    BROWSER = Browser()

    class A(object):
        def __init__(self, name):
            self.name=name
            #global BROWSER
            self.browser = BROWSER

        def chdir(self, d):
            print 'PREV:', self.name, self.browser.getcwd()
            self.browser.chdir(d)
            print 'NEW:', self.name, self.browser.getcwd()

    testA = A('a')
    testB = A('b')

    testA.chdir('..')
    testB.chdir('..')
    testA.browser.goup()
    print testA.name, testA.browser.cwd
    print testB.name, testB.browser.cwd
    testB.browser.goback()
    print testA.name, testA.browser.cwd