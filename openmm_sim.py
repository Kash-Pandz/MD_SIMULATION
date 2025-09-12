import argparse
from loguru import logger
import openmm as mm
from openmm import app, unit
from openmm.app import (
    Simulation,
    PDBFile,
    StateDataReporter,
    DCDReporter,
    AmberPrmtopFile,
    AmberInpcrdFile,
)


def parse_args():
    parser = argparse.ArgumentParser(description="MD Simulation with OpenMM")
    parser.add_argument("--prmtop", required=True, help="Amber prmtop file")
    parser.add_argument("--inpcrd", required=True, help="Amber inpcrd file")
    parser.add_argument("--temperature", type=float, default=300.0, help="Target temperature for production (K)")
    parser.add_argument("--platform", default="CPU", help="OpenMM platform: CPU, CUDA, OpenCL")
    parser.add_argument("--checkpoint_npt", default="npt.chk", help="Checkpoint file for NPT equilibration")
    parser.add_argument("--checkpoint_prod", default="prod.chk", help="Checkpoint file for production run")
    parser.add_argument("--nvt_steps", nargs=3, type=int, default=[100000, 100000, 100000],
                        help="Number of steps for each NVT stage (3 integers)")
    parser.add_argument("--npt_steps", type=int, default=500000, help="Number of steps for NPT equilibration")
    parser.add_argument("--prod_steps", type=int, default=50000000, help="Number of steps for production run")
    return parser.parse_args()


def load_amber_files(prmtop, inpcrd):
    """Parse AMBER topology and coordinate files."""
    topology = AmberPrmtopFile(prmtop)
    coordinates = AmberInpcrdFile(inpcrd)
    if coordinates.boxVectors is not None:
        coordinates = AmberInpcrdFile(inpcrd, periodicBoxVectors=coordinates.boxVectors)
    return topology, coordinates


def build_system(topology):
    """Build system."""
    return topology.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1*unit.nanometer,
        constraints=app.HBonds
    )


def set_integrator(temperature, friction=1/unit.picosecond, timestep=0.002):
    """Initialise integrator.""""
    return mm.LangevinMiddleIntegrator(temperature*unit.kelvin, friction, timestep*unit.picoseconds)


def apply_restraint(system, topology, positions, k=0.0, chain_indices=None, mode="all"):
    """A function to apply atom group restraints.""" 
    if k <= 0.0:
        return None

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
            logger.warning(f"Chain index {chain_index} out of range. Skipping...")
            continue
        for atom in chains[chain_index].atoms():
            if mode == "heavy" and atom.element.symbol == "H":
                continue
            if mode == "backbone" and atom.name not in ("CA", "C", "N", "O"):
                continue
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, pos.tolist())
            n_restrained += 1

    if n_restrained > 0:
        system.addForce(restraint)
        logger.info(f"Applied restraints to {n_restrained} atoms (mode={mode}) with k={k:.2f} kcal/mol/nm²")
    return restraint


### Energy Minimisation ###

def run_min(topology, positions, temperature, min_prefix, restraint_k=0.0, platform_name="CPU", platform_props=None):
    if platform_props is None:
        platform_props = {}
    system = build_system(topology)
    apply_restraint(system, topology, positions, k=restraint_k, mode="all")
    integrator = set_integrator(temperature)
    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = Simulation(topology.topology, system, integrator, platform=platform, platformProperties=platform_props)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True)
    with open(f"{min_prefix}_min.pdb", "w") as f:
        PDBFile.writeFile(topology.topology, state.getPositions(), f)
    return simulation


### NVT Equilibration ###

def run_nvt(topology, positions, velocities, nvt_stages, platform_name="CPU", platform_props=None):
    if platform_props is None:
        platform_props = {}
    state = None
    for i, (temp, steps, restraint_k) in enumerate(nvt_stages):
        system = build_system(topology)
        apply_restraint(system, topology, positions, k=restraint_k, mode="all")
        integrator = set_integrator(temp)
        platform = mm.Platform.getPlatformByName(platform_name)
        simulation = Simulation(topology.topology, system, integrator, platform=platform, platformProperties=platform_props)
        simulation.context.setPositions(positions)
        if velocities is not None:
            simulation.context.setVelocities(velocities)
        else:
            simulation.context.setVelocitiesToTemperature(temp*unit.kelvin)
        simulation.reporters.append(StateDataReporter(f"NVT_stage{i+1}_{temp}K.log", 500,
                                                      step=True, temperature=True, potentialEnergy=True, kineticEnergy=True))
        simulation.reporters.append(DCDReporter(f"NVT_stage{i+1}_{temp}K.dcd", 500))
        logger.info(f"NVT stage {i+1}: T={temp}K, restraint_k={restraint_k} kcal/mol/nm², steps={steps}")
        simulation.step(steps)
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
    return state


### NPT Equilibration ###

def run_npt(topology, positions, velocities, steps, temperature=300, restraint_k=0.2, timestep=0.002, platform_name="CPU", platform_props=None, checkpoint_file=None):
    if platform_props is None:
        platform_props = {}
    system = build_system(topology)
    apply_restraint(system, topology, positions, k=restraint_k, mode="all")
    system.addForce(mm.MonteCarloBarostat(1*unit.bar, temperature*unit.kelvin))
    integrator = set_integrator(temperature, timestep=timestep)
    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = Simulation(topology.topology, system, integrator, platform=platform, platformProperties=platform_props)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    simulation.reporters.append(StateDataReporter(f"NPT_{temperature}K.log", 500,
                                                  step=True, temperature=True, potentialEnergy=True,
                                                  kineticEnergy=True, volume=True, density=True, pressure=True))
    simulation.reporters.append(DCDReporter(f"NPT_{temperature}K.dcd", 500))
    if checkpoint_file:
        try:
            with open(checkpoint_file, "rb") as f:
                simulation.context.loadCheckpoint(f.read())
            logger.info(f"Loaded checkpoint {checkpoint_file}")
        except FileNotFoundError:
            logger.info(f"No checkpoint file found!.")
    simulation.step(steps)
    if checkpoint_file:
        with open(checkpoint_file, "wb") as f:
            f.write(simulation.context.createCheckpoint())
        logger.info(f"Saved checkpoint {checkpoint_file}")
    return simulation.context.getState(getPositions=True, getVelocities=True)


### Production NPT ###

def run_prod(topology, positions, velocities, steps, temperature=300, timestep=0.002, platform_name="CPU", platform_props=None, checkpoint_file=None):
    if platform_props is None:
        platform_props = {}
    system = build_system(topology)
    system.addForce(mm.MonteCarloBarostat(1*unit.bar, temperature*unit.kelvin))
    integrator = set_integrator(temperature, timestep=timestep)
    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = Simulation(topology.topology, system, integrator, platform=platform, platformProperties=platform_props)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    simulation.reporters.append(StateDataReporter(f"PROD_NPT_{temperature}K.log", 500,
                                                  step=True, temperature=True, potentialEnergy=True,
                                                  kineticEnergy=True, volume=True, density=True, pressure=True))
    simulation.reporters.append(DCDReporter(f"PROD_NPT_{temperature}K.dcd", 500))
    if checkpoint_file:
        try:
            with open(checkpoint_file, "rb") as f:
                simulation.context.loadCheckpoint(f.read())
            logger.info(f"Loaded checkpoint {checkpoint_file}")
        except FileNotFoundError:
            logger.info(f"No checkpoint file found!.")
    simulation.step(steps)
    if checkpoint_file:
        with open(checkpoint_file, "wb") as f:
            f.write(simulation.context.createCheckpoint())
        logger.info(f"Saved checkpoint {checkpoint_file}")
    return simulation.context.getState(getPositions=True, getVelocities=True)


def main():
    args = parse_args()
    topology, coordinates = load_amber_files(args.prmtop, args.inpcrd)
    positions = coordinates.positions
    velocities = None

    logger.info("Starting energy minimization...")
    sim = run_min(topology, positions, args.temperature, "min", restraint_k=5.0, platform_name=args.platform)
    state = sim.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions()
    velocities = state.getVelocities()

    # NVT equilibration
    nvt_stages = [
        (100, args.nvt_steps[0], 5.0),
        (200, args.nvt_steps[1], 2.0),
        (300, args.nvt_steps[2], 0.5),
    ]
    logger.info("Starting NVT equilibration...")
    state = run_nvt(topology, positions, velocities, nvt_stages, platform_name=args.platform)
    positions = state.getPositions()
    velocities = state.getVelocities()

    # NPT equilibration
    logger.info("Starting NPT equilibration...")
    state = run_npt(topology, positions, velocities, steps=args.npt_steps, temperature=args.temperature,
                    timestep=0.002, platform_name=args.platform, checkpoint_file=args.checkpoint_npt)
    positions = state.getPositions()
    velocities = state.getVelocities()

    # Production
    logger.info("Starting production NPT...")
    run_prod(topology, positions, velocities, temperature=args.temperature, steps=args.prod_steps,
             timestep=0.002, platform_name=args.platform, checkpoint_file=args.checkpoint_prod)
    logger.info("Production finished.")


if __name__ == "__main__":
    main()
