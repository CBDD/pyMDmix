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
    def __init__(self, queue, stats=None):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        self.stats = stats      # (submitted, done, t0) shared counters used to report progress
#        self.sid = sid #sessionid

    def _report(self, label):
        "Count a finished job and print the progress line."
        if not self.stats: return
        submitted, done, t0 = self.stats
        with done.get_lock():
            done.value += 1
            n = done.value
        total = max(submitted.value, n)
        elapsed = time.time() - t0.value
        print("[%d/%d, %3d%%, %s elapsed] %s finished" % (n, total, 100 * n // total, _fmt_time(elapsed), label or 'job'), flush=True)

    def run(self):
        while True:
            # Get next task or exit if None
            exe = self.queue.get()
            if exe is None: break
            if isinstance(exe[0], str): 
                # Shell comands (cmd, path[, label])
                cmd, path = exe[0], exe[1]
                label = exe[2] if len(exe) > 2 else None
                exit = sub.Popen(cmd, shell=True, cwd=path).wait()
                if exit: print("cmd exited with non-zero code: %s"%cmd, file=sys.stderr)
                self._report(label or cmd)
            else:
                # Methods and functions
                fx, args, kwargs = exe
                fx(*args,**kwargs)
                self._report(getattr(fx, '__name__', 'function'))

def _fmt_time(seconds):
    "Format seconds as H:MM:SS."
    seconds = int(seconds)
    return "%d:%02d:%02d" % (seconds // 3600, (seconds % 3600) // 60, seconds % 60)

class Executor(object):
    def __init__(self, nthreads=1, twait=1):
        self.log = logging.getLogger("Executor")
        self.nthreads = nthreads
        self.procList = []
        self.twait = twait
        self.threads = []
        self.queue = multiprocessing.Queue()
        self.ncalls = multiprocessing.Manager().Value('i',0)
        # Shared job counters (submitted, finished, start time) for progress reporting and for waiting
        # until the jobs have really finished (an empty queue only means they were picked up)
        self.submitted = multiprocessing.Value('i', 0)
        self.done = multiprocessing.Value('i', 0)
        self.t0 = multiprocessing.Value('d', time.time())
        self.log.debug("Init executor nthreads: %d"%nthreads)
    
    def __del__(self):
        try:
            [self.queue.put(None) for _ in range(self.nthreads)]
        except Exception:   # i.e. interpreter shutting down: Python >= 3.13 does not allow starting threads then
            pass
    
    def start(self):
        self.ncalls.value += 1
        self.log.debug("Registering new caller. ncalls: %d"%self.ncalls.value)
        if self.ncalls.value == 1:
            # first call: reset the progress counters and init threads
            self.submitted.value = 0
            self.done.value = 0
            self.t0.value = time.time()
            self.__setThreads()
                
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
        [self.threads.append(ExecutorThread(self.queue, (self.submitted, self.done, self.t0))) for _ in range(self.nthreads)]
        [t.start() for t in self.threads]
#        [t.join() for t in self.threads]
            
    def changeNthreads(self, nthreads):
        "Change de number of processes to run in parallel"
        self.nthreads = nthreads or 1
        self.log.debug("Changed executor nthreads: %d"%nthreads)

    def submitCmd(self, fx=None, cmd=None, path=None, label=None, *args, **kwargs):
        """
        Submit command. Will wait for any slot in the processes to empty to launch current.
        
        :arg str cmd: String with full execution command. 
        :arg str path: Path where to execute the program. Path will be passed to ``subprocess.Popen`` as *cwd* argument.
        :arg str label: Short name of the job (for example "ETAO_2 md7") shown in the progress line when it finishes.
        """
        if self.ncalls.value == 0: self.start()
        if cmd:
            self.log.debug("Submitting to queue: %s (%s)"%(cmd,path or 'cwd'))
            with self.submitted.get_lock(): self.submitted.value += 1
            self.queue.put_nowait((cmd,path,label))
        elif fx:
            with self.submitted.get_lock(): self.submitted.value += 1
            self.queue.put_nowait((fx, args, kwargs))
        
    def waitJobCompletion(self):
        "Wait until all submitted jobs have finished (not only until they have left the queue)."
        while self.done.value < self.submitted.value:
            # Do not wait forever if every worker started by this process died
            if self.threads and not any(t.is_alive() for t in self.threads):
                self.log.warning("Executor workers are not running; %d job(s) did not finish."%(self.submitted.value - self.done.value))
                break
            time.sleep(self.twait)
    

if __name__ == "__main__":
    print("Hello World")
