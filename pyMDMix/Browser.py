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
__date__ ="$Mar 17, 2014 3:11:29 PM$"

##### DIRECTORY BROWSING TOOL
import os
import os.path as osp
import logging

class BrowserError( Exception):
    pass

class WrongDirectory(BrowserError):
    pass

class Browser(object):
    "Wrap os functions to work with project and recognize folders"
    def __init__(self, project=None):
        self.log = logging.getLogger("Browser")
        self.cwd = os.getcwd()
        self.prevDir = os.getcwd()
        self.home = os.getcwd() # Project Home path. Until set, init with cwd
        self.project = project
        if self.project: self.setHome(self.project.projectPath)

    def setProject(self, project):
        self.project = project
        self.setHome(project.projectPath)

    def chdir(self, newdir):
        if not osp.exists(newdir): raise WrongDirectory, "Directory %s does not exists"%(newdir)
        os.chdir(newdir)
        self.update()

    def update(self):
        self.prevDir = self.cwd
        self.cwd = os.getcwd()
        self.log.debug("Directory change from %s to %s"%(self.prevDir, self.cwd))
        return self.cwd

    def goMD(self):
        self.goHome()
        self.chdir('MD')

    def goup(self):
        "Like cd.."
        self.chdir(os.pardir)

    def setHome(self, homedir=None):
        "Move to project home folder and do setHome there"
        if homedir and os.getcwd() != homedir:
            if not osp.exists(homedir): raise WrongDirectory, "Directory %s does not exists"%(homedir)
            os.chdir(homedir)
        self.home = self.update()

    def gotoReplica(self, replica, subfolder=''):
#        if self.cwd != self.home: self.goHome()
        path = replica.path
        if subfolder: path = osp.join(path, subfolder)
        self.chdir(path)

    def goHome(self):
        self.chdir(self.home)

    def goback(self):
        "Like cd -"
        self.chdir(self.prevDir)

    def getcwd(self):
        cwd = os.getcwd()
        if self.cwd != cwd: self.update()
        return cwd


if __name__ == "__main__":
    print "Hello World"
