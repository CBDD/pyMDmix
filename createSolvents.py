import os
def createSolventDB(path):
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

createSolventDB(os.path.join(os.curdir, "pyMDMix/data/solventlib" ))