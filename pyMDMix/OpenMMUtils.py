import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u

import sys

def applyHarmonicPositionalRestraints(system, forceConstantInKcalPerMolePerAngSquared,
                                      inpcrd, indexOfAtomsToBeModified):
    """ This is essentially mimicking AMBER's restraint_wt"""

    forceConstant = u.Quantity(value=forceConstantInKcalPerMolePerAngSquared,
               unit=u.kilocalorie/(u.mole * u.angstrom * u.angstrom))

    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")

    force.addGlobalParameter("k",
       forceConstant.in_units_of(u.kilojoule/(u.mole * u.nanometer * u.nanometer )))

    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for i in indexOfAtomsToBeModified:
        force.addParticle(i, inpcrd.getPositions()[i])

    system.addForce(force)

def setContextFromRst(simulation, inpcrd):
    """ Restore the simulation's context from all the information
    present in an AMBER rst7 file.
    """
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


