import argparse
from loguru import logger
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds, PDBFile
import parmed as pmd


def parse_args():
    parser = argparse.ArgumentParser(description="Convert AMBER files to GROMACS")
    parser.add_argument("--do_minimise", action="store_true", help="Run energy minimisation using OpenMM")
    parser.add_argument("--prmtop", type=str, default="SYSTEM.prmtop", help="AMBER prmtop file")
    parser.add_argument("--inpcrd", type=str, default="SYSTEM.inpcrd", help="AMBER inpcrd file")
    parser.add_argument("--gmx_gro", type=str, default="SYSTEM.gro", help="Output GROMACS .gro file")
    parser.add_argument("--gmx_top", type=str, default="SYSTEM.top", help="Output GROMACS .top file")
    parser.add_argument("--min_pdb", type=str, default=None, help="Optional energy minimised PDB file")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def configure_logging(verbose: bool):
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(lambda msg: print(msg, end=""), level=level)


def minimise_and_convert_to_gromacs(prmtop_file: str, inpcrd_file: str, gmx_gro: str, gmx_top: str, min_pdb: str = None):
    logger.info("Loading AMBER files for energy minimisation")
    inpcrd = AmberInpcrdFile(inpcrd_file)
    prmtop = AmberPrmtopFile(prmtop_file)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=Hbonds
    )

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 
        1 / unit.picosecond, 
        0.002 * unit.picoseconds
    )
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    logger.info("Minimizing energy...")
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()

    if min_pdb:
        logger.info(f"Writing minimized PDB to {min_pdb}")
        with open(min_pdb, 'w') as f:
            PDBFile.writeFile(simulation.topology, positions, f)

    logger.info(f"Writing GROMACS files: {gmx_gro}, {gmx_top}")
    structure = pmd.openmm.load_topology(simulation.topology, system, positions)
    structure.save(gmx_top, overwrite=True)
    structure.save(gmx_gro, overwrite=True)

    logger.info("AMBER to GROMACS complete.")

  
def main():

    args= parser.parser_args()
    configure_logging(args.verbose)

    minimise_and_convert_to_gromacs(
        prmtop_file=args.prmtop,
        inpcrd_file=args.inpcrd,
        gmx_gro=args.gmx_gro,
        gmx_top=args.gmx_top,
        min_pdb=args.min_pdb
    )

if __name__ == "__main__":
    main()
