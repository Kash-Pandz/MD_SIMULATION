import argparse
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import (
    PDBFile,
    ForceField, 
    Modeller, 
    PME, 
    HBonds, 
    Simulation
)
from openmm import LangevinMiddleIntegrator
from openmm.unit import kelvin, picosecond, femtosecond, nanometers, molar
import parmed as pmd
from loguru import logger


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare system for using OpenMM.")
    parser.add_argument("--input_pdb", required=True, help="Input .pdb file")
    parser.add_argument("--output_pdb", required=True, help="Output .pdb file")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for Gromacs files")
    parser.add_argument("--forcefield", default="amber14-all.xml", help="OpenMM force field")
    parser.add_argument("--ph", type=float, default=7.0, help="System pH")
    parser.add_argument("--water", default="amber14/tip3pfb.xml", help="OpenMM water model")
    parser.add_argument("--padding", type=float, default=1.0, help="Box padding (nm)")
    parser.add_argument("--ionic_strength", type=float, default=0.15, help="Ionic strength (mM)")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument("--timestep", type=float, default=2.0, help="Time Step (fs)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def configure_logging(verbose: bool):
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(lambda msg: print(msg, end=""), level=level)


def fix_pdb_structure(input_pdb: Path) -> PDBFixer:
    logger.info(f"Loading PDB file: {input_pdb.name}")
    fixer = PDBFixer(filename=str(input_pdb))
    logger.info("Fixing Input Structure...")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.addMissingAtoms()
    return fixer


def solvate_and_minimize_system(
        fixer: PDBFixer,
        output_pdb: Path,
        gmx_prefix: str,
        forcefield: str,
        ph: float, 
        water_model: str,
        padding: float, 
        ionic_strength: float, 
        temperature: float,
        timestep: float
):
    logger.info("Loading force field and water model...")
    forcefield = ForceField(forcefield, water_model)
    
    modeller = Modeller(fixer.topology, fixer.positions)

    logger.info(f"Adding hydrogens at pH {ph}...")
    modeller.addHydrogens(forcefield, pH=ph)

    logger.info("Adding solvent and ions...")
    modeller.addSolvent(
        forcefield,
        model=water_model,
        padding=padding * nanometers,
        ionicStrength=ionic_strength * molar,
    )

    logger.info("Creating system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometers,
        constraints=HBonds,
        rigidWater=False
    )
    
    integrator = LangevinMiddleIntegrator(
        temperature * kelvin,
        1 / picosecond,
        timestep * femtosecond
    )

    logger.info(f"Minimizing energy...")
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()

    logger.info(f"Writing minimized PDB to {output_pdb.name}")
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, positions, f)

    logger.info(f"Writing Gromacs files to {gmx_prefix}.gro and {gmx_prefix}.top")
    structure = pmd.openmm.load_topology(modeller.topology, system, positions)
    structure.save(f"{gmx_prefix}.top", overwrite=True)
    structure.save(f"{gmx_prefix}.gro", overwrite=True)

    logger.success(f"System prepared and Gromacs files written.")


def main():
    args = parse_args()
    configure_logging(args.verbose)

    pdb_path = Path(args.input_pdb)
    output_pdb = Path(args.output_pdb)

    fixer = fix_pdb_structure(pdb_path)
    solvate_and_minimize_system(
        fixer,
        output_pdb,
        args.output_prefix,
        args.forcefield,
        args.ph,
        args.water,
        args.padding,
        args.ionic_strength,
        args.temperature, 
        args.timestep,
    )


if __name__ == "__main__":
    main()
