import argparse
from loguru import logger 
import parmed as pmd
import openmm as mm
from openmm import LangevinMiddleIntegrator, app, unit



def load_amber_files(prmtop, inpcrd):
  topology = app.AmberPrmtopFile(prmtop)
  coordinates = app.AmberInpcrdFile(crd_file)
  return topology, coordinates
  

def build_system(top):
  return topology.createSystem(nobondedMethod=app.PME,
                               nonbondedCutoff=1.0 * unit.nanometer,
                               constraints=app.HBonds
                              )

def apply_restraint(topology, positions, k, chain_indices=None, mode="all"):
    restraint = mm.CustomExternalForce("k * (periodincdistance(x, y, z, x0, y0, z0))^2")
    restraint.addGlobalParameter("k", k * unit.kilocalorie_per_mole / unit.angstrom**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    chains = list(topology.chains())
    if chain_indices is None:
        chain_indices = range(len(chains))

    n_restrained = 0
    for chain_index in chain_indices:
        if chain_index >= len(chains):
            logger.warning(f"Chain index {chain_index} out of range. Skipping.")
            continue
        for atom in chains[chain_index].atoms():
            if mode == "heavy" and atom.element.symbol == "H":
                continue
            if mode == "backbone" and atom.name not in ("CA", "C", "N"):
                continue
            restraint.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
            n_restrained += 1

    system.addForce(restraint)
    logger.info(f"Applied restraints to {n_restrained} atoms (mode={mode}) with k={k: .2f}")
    return restraint



