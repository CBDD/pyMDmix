#! /usr/bin/env python
##  ------------------------ pyMDMix -----------------------------------
##                  http://mdmix.sourceforge.net
##  --------------------------------------------------------------------
## 
##  Software for preparation, analysis and quality control
##  of solvent mixtures molecular dynamics.
## 
##  Copyright (C) 2014 dalvarez
## 
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
##  Please cite your use of pyMDMix in published work:
## 
##              TOBEPUBLISHED
## 
##  --------------------------------------------------------------------

__author__="dalvarez"
__date__ ="$11-mar-2014 17:33:33$"

import argparse
import sys
import pyMDMix
import numpy as npy

from pyMDMix import MDMixError

class MDMixClient(object):
    def __init__(self):
        self.actions = {}
        self.parser = self.createParser()

    def createParser(self):
        from pyMDMix.Commands import Create, Info, Add, Remove, Queue, Plot, Analyze, Tools
        parser = argparse.ArgumentParser("mdmix")
        parser.add_argument("--log", action="store", dest="logfile", help="Logging file. Default: output to stdout")
        parser.add_argument("--debug", action="store_true", dest="debug", default=False, help="Print debugging info")
        subparsers = parser.add_subparsers(help='commands', dest='command')
        for action in [
            Create.Create(),
            Info.Info(),
            Add.Add(),
            Remove.Remove(),
            Queue.Queue(),
            Plot.Plot(),
            Analyze.Analyze(),
            Tools.Tools()
        ]:
            self.actions[action.cmdstring] = action
            action.create_parser(subparsers)
        return parser

    def header(self):
        
        return """
        ==========================================================
        ||              pyMDMix User Interface                  ||
        ==========================================================
        ||  Author: Daniel Alvarez-Garcia                       ||
        ||  Version : %s                                     
        ==========================================================
        
        """%pyMDMix.__version__

    def createRootLogger(self, level, logFile=None):
        pyMDMix.setLogger(level, logFile)

    def run(self):
        print self.header()
        import time
        import sys
        t0 = time.time()
        parserargs =  self.parser.parse_args()
        command = parserargs.command
        
        # If logging file, redirect stdout and stderr to file
        if parserargs.logfile:
            sys.stderr = open(parserargs.logfile, 'w+')
            sys.stdout = sys.stderr
            print ' '.join(sys.argv)

        if parserargs.debug: level = 'DEBUG'
        else: level = 'INFO'
        self.createRootLogger(level,parserargs.logfile)

        #Study the different actions
        self.actions[command].action(parserargs)

        sys.stdout.flush()
        print "Total execution time: %.3fs"%(time.time()-t0)

def main():
    client = MDMixClient()
    client.run()

if __name__ == "__main__":
    try: 
        main()
    except KeyboardInterrupt: print "Forcing MDMix UI exit!"
    sys.exit(0)
