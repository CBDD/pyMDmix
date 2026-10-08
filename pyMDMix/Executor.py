#  ------------------------ pyMDMix ----------------------------------- 
#                  http://mdmix.sourceforge.net
#  -------------------------------------------------------------------- 
# 
#  Software for preparation, analysis and quality control
#  of solvent mixtures molecular dynamics.
# 
#  Copyright (C) 2014 daniel
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Please cite your use of pyMDMix in published work:
# 
#              TOBEPUBLISHED
# 
#  -------------------------------------------------------------------- 

__author__="daniel"
__date__ ="$Mar 17, 2014 3:12:43 PM$"

import sys
import os
import time
import logging
import subprocess as sub
import multiprocessing

class ExecutorThread(multiprocessing.Process):
    def __init__(self, queue):
        multiprocessing.Process.__init__(self)
        self.queue = queue
#        self.sid = sid #sessionid
        
    def run(self):
        while True:
            # Get next task or exit if None
            exe = self.queue.get()
            if exe is None: break
            if isinstance(exe[0], str): 
                # Shell comands
                cmd, path = exe
                exit = sub.Popen(cmd, shell=True, cwd=path).wait()
                if exit: print("cmd exited with non-zero code: %s"%cmd, file=sys.stderr)
            else:
                # Methods and functions
                fx, args, kwargs = exe
                fx(*args,**kwargs)

class Executor(object):
    def __init__(self, nthreads=1, twait=1):
        self.log = logging.getLogger("Executor")
        self.nthreads = nthreads
        self.procList = []
        self.twait = twait
        self.threads = []
        self.queue = multiprocessing.Queue()
        self.ncalls = multiprocessing.Manager().Value('i',0)
        self.log.debug("Init executor nthreads: %d"%nthreads)
    
    def __del__(self):
        [self.queue.put(None) for _ in range(self.nthreads)]
    
    def start(self):
        self.ncalls.value += 1
        self.log.debug("Registering new caller. ncalls: %d"%self.ncalls.value)
        if self.ncalls.value == 1: self.__setThreads() # init threads if first call
                
    def terminate(self):
        self.ncalls.value -= 1
        if self.ncalls.value <= 0: # Terminate threads if last call
            self.log.debug("Terminating threads")
            [self.queue.put(None) for _ in range(self.nthreads)]
            self.threads = []
#            [t.terminate() for t in self.threads]
#            [t.join() for t in self.threads]
#            self.__terminateThreads()
        
    def __setThreads(self):
        "Start threads acording to self.maxproc"
        self.log.debug("Starting threads")
        [self.threads.append(ExecutorThread(self.queue)) for _ in range(self.nthreads)]
        [t.start() for t in self.threads]
#        [t.join() for t in self.threads]
            
    def changeNthreads(self, nthreads):
        "Change de number of processes to run in parallel"
        self.nthreads = nthreads or 1
        self.log.debug("Changed executor nthreads: %d"%nthreads)

    def submitCmd(self, fx=None, cmd=None, path=None, *args, **kwargs):
        """
        Submit command. Will wait for any slot in the processes to empty to launch current.
        
        :arg str cmd: String with full execution command. 
        :arg str path: Path where to execute the program. Path will be passed to ``subprocess.Popen`` as *cwd* argument.
        """
        if self.ncalls.value == 0: self.start()
        if cmd:
            self.log.debug("Submitting to queue: %s (%s)"%(cmd,path or 'cwd'))
            self.queue.put_nowait((cmd,path))
        elif fx:
            self.queue.put_nowait((fx, args, kwargs))
        
    def waitJobCompletion(self):
        "Wait until all jobs in the queue are done"
        while not self.queue.empty(): 
#            self.log.debug("Wating completion")
            time.sleep(self.twait)
    

if __name__ == "__main__":
    print("Hello World")
