#! /usr/bin/python

# Given a configuration file, add a solvent to the solvent database

__author__="daniel"
__date__ ="$Oct 11, 2012 4:43:37 PM$"

if __name__ == "__main__":
    import sys
    from pyMDMix.Solvents import SolventManager
    import glob
    configfilelist = glob.glob('*.config')
    maker=SolventManager()
    for conffile in configfilelist:
        print('Adding solvent from %s'%conffile)
        maker.saveSolvent(maker.createSolvent(conffile), db='SOLVENTS.db', createEmpty=True)
    print("DONE")
