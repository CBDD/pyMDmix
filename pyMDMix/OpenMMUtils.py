import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u

import sys

def getHeavyAtomIndexExceptWat(modeller):
    """ Return the index of all heavy atoms, except solvent, from a modeller
     object.
    """
    # Get ligand index
    HeavyAtomIndex = []
    for residue in modeller.topology.residues():
        if (residue.name != "HOH"):
            # print residue.name
            for atom in residue.atoms():
                if (atom.element != app.element.hydrogen):
                    # print atom.name, atom.index, atom.element
                    HeavyAtomIndex.append(atom.index)

    return HeavyAtomIndex


def applyHarmonicPositionalRestraints(system, forceConstant, inpcrd,
                                      indexOfAtomsToBeModified):

    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")

    force.addGlobalParameter("k", forceConstant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for i in indexOfAtomsToBeModified:
        force.addParticle(i, inpcrd.getPositions()[i])

    system.addForce(force)

def setContextFromRst(simulation, inpcrd):
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setVelocities(inpcrd.velocities)

    x = mm.Vec3(inpcrd.getBoxVectors()[0][0].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[0][1].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[0][2].value_in_unit(u.nanometers))

    y = mm.Vec3(inpcrd.getBoxVectors()[1][0].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[1][1].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[1][2].value_in_unit(u.nanometers))

    z = mm.Vec3(inpcrd.getBoxVectors()[2][0].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[2][1].value_in_unit(u.nanometers),
            inpcrd.getBoxVectors()[2][2].value_in_unit(u.nanometers))

    simulation.context.setPeriodicBoxVectors(x,y,z)


