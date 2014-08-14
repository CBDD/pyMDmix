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
__date__ ="$Mar 17, 2014 5:16:35 PM$"

__all__ = []

import pkgutil
import inspect

for loader, name, is_pkg in pkgutil.walk_packages(__path__):
    module = loader.find_module(name).load_module(name)
#    globals()[name] = module
#    __all__.append(name)
    for name, value in inspect.getmembers(module):
        if name.startswith('__'):
            continue

        globals()[name] = value
        __all__.append(name)

del inspect, pkgutil, name, loader, is_pkg