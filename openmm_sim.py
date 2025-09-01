import argparse
from loguru import logger 
import parmed as pmd
import openmm as mm
from openmm import LangevinMiddleIntegrator, app, unit



def load_amber_files(prmtop, inpcrd):
    topology = app.AmberPrmtopFile(prmtop)
    coordinates = app.AmberInpcrdFile(crd_file)
    return topology, coordinates


def build_system(topology):
    return topology.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer, constraints=app.HBonds)


def set_integrator(temperature, friction=1/unit.picosecond, timestep=0.002*unit.picoseconds):
    return mm.LangevinMiddleIntegrator(temperature*unit.kelvin, friction, timestep)


def create_simulation(topology, system, integrator, positions=None, velocities=None, platformName="CUDA", platformProps=None):

    # Set platform properties
    platform = mm.Platform.getPlatformByName(platformName)
    if platformProps is None:
        platformProps = {}

    # Create simulation 
    simulation = app.Simulation(topology, system, integrator, platform, platformProps)

    # Assign positions and velocities if provided
    if positions is not None:
        simulation.context.setPositions(positions)
    if velocities is not None:
        simulation.context.setVelocities(velocities)

    return simulation 


def apply_restraint(topology, positions, k, chain_indices=None, mode="all", logger=None):
    restraint = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k", k * unit.kilocalorie_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    chains = list(topology.topology.chains())
    if chain_indices is None:
      chain_indices = range(len(chains))
 
    n_restrained = 0
    for chain_index in chain_indices:
        if chain_index >= len(chains):
            msg = f"Chain index {chain_index} out of range. Skipping."
            (logger.warning(msg) if logger else print(msg))
            continue
        for atom in chains[chain_index].atoms():
            if mode == "heavy" and atom.element.symbol == "H":
                continue
            if mode == "backbone" and atom.name not in ("CA", "C", "N", "O"):
                continue
            restraint.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
            n_restrained += 1

    system.addForce(restraint)
    msg = f"Applied restraints to {n_restrained} atoms (mode={mode}) with k={k:.2f} kcal/mol/nm^2"
    (logger.info(msg) if logger else print(msg))
    return restraint


### Energy Minimisation ###

def minimise(topology, positions, temperature, min_prefix, restraint_k=0.0, restraint_mode="all", platformName="CPU", platformProps={}):
    system = build_system(topology)
  
    if restraint_k > 0:
        apply_restraint(system, topology, positions, restraint_k, mode=restraint_mode, logger=logger)

    integrator = set_integrator(temperature)

    platform = mm.Platform.getPlatformByName(platformName)
    simulation = app.Simulation(topology.topology, system, integrator, platform, platformProps)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()

    # Save energy minimised pdb file
    state = simulation.context.getState(getPositions=True)
    with open(f"{min_prefix}_min.pdb", "w") as f:
        app.PDBFile.writeFile(topology.topology, state.getPositions(), f)
    return simulation




### NVT ###

def run_nvt_eq(simulation, positions, velocity=None, temp=300.0, stepsPerStage=5000, restraint_modes=None, k_values=None):

    system = build_system(topology)

    if restraint_modes:
        pos_restraints = apply_restraint(topology, positions, k, chain_indices=None, mode="all", logger=None)
        system.addForce(pos_restraints)

    thermostat = mm.AndersenThermostat(temperature * unit.kelvin, 1/unit.picosecond)
    system.addForce(thermostat)

    integrator = set_integrator(temperature)
    platform = mm.Platform.getPlatformByName(platformName)
    simulation = create_simulation(topology, 
                                   system, 
                                   integrator, 
                                   positions=None, 
                                   velocities=None, 
                                   platformName="CUDA", 
                                   platformProps=None
                                  )

    if velocity is None:
        simulation.context.setVelocitiesToTemperature(temp * unit.kelvin, 1)

    fileName = "{}-NVT-{}.log".format(jobname, int(temperature))
    reporterLog = app.StateDataReporter(fileName, logFrequency, step=True, temperature=True, kineticEnergy=True,
                                                                    potentialEnergy=True, volume=True, separator=' ')
    simulation.reporters.append(reporterLog)

    fileName = "{}-NVT-{}.dcd".format(jobname, int(temperature))
    simulation.reporters.append(app.DCDReporter(fileName, trajFrequency))
    simulation.step(steps)

    return simulation.context.getState(getPositions=True, getVelocities=True)
    


def run_nvt_eq(
    topology,
    positions,
    velocity=None,
    temp=300.0,
    steps=5000,
    restraint_modes=None,
    k_values=None,
    platformName="CUDA",
    platformProps=None,
    jobname="simulation",
    logFrequency=1000,
    trajFrequency=1000
):
    # Build system
    system = build_system(topology)

    # Apply restraints if requested
    if restraint_modes and k_values:
        pos_restraints = apply_restraint(
            topology,
            positions,
            k_values,
            chain_indices=None,
            mode=restraint_modes,
            logger=None
        )
        system.addForce(pos_restraints)

    # Add thermostat
    thermostat = mm.AndersenThermostat(temp * unit.kelvin, 1/unit.picosecond)
    system.addForce(thermostat)

    # Integrator & simulation
    integrator = set_integrator(temp)
    platform = mm.Platform.getPlatformByName(platformName)
    if platformProps is None:
        platformProps = {}
    simulation = app.Simulation(topology.topology, system, integrator, platform, platformProps)

    # Initial conditions
    if velocity is not None:
        simulation.context.setVelocities(velocity)
    else:
        simulation.context.setVelocitiesToTemperature(temp * unit.kelvin, 1)

    # Reporters
    fileName = f"{jobname}-NVT-{int(temp)}.log"
    reporterLog = app.StateDataReporter(
        fileName,
        logFrequency,
        step=True,
        temperature=True,
        kineticEnergy=True,
        potentialEnergy=True,
        volume=True,
        separator=' '
    )
    simulation.reporters.append(reporterLog)

    fileName = f"{jobname}-NVT-{int(temp)}.dcd"
    simulation.reporters.append(app.DCDReporter(fileName, trajFrequency))

    # Run simulation
    simulation.step(steps)

    return simulation.context.getState(getPositions=True, getVelocities=True)















