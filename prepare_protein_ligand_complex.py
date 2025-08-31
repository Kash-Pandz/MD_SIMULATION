import argparse
import subprocess
import sys
from pathlib import Path
from loguru import logger
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds, PDBFile
import parmed as pmd


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare protein-ligand system.")
    parser.add_argument("--input_pdb", required=True, help="Input protein-ligand complex")
    parser.add_argument("--ligand_resname", default="LIG", help="Ligand resname")
    parser.add_argument("--ligand_charge", type=int, default=0, help="Net charge of ligand")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--do_minimise", action="store_true", help="Run energy minimisation with OpenMM")
    parser.add_argument("--gmx_gro", type=str, default="system.gro", help="Output GROMACS (.gro) file")
    parser.add_argument("--gmx_top", type=str, default="system.top", help="Output GROMACS (.top) file")
    parser.add_argument("--min_pdb", type=str, default="em.pdb", help="Optional energy minimised file")
    return parser.parse_args()


def configure_logging(verbose: bool):
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(sys.stderr,
               level=level,
               format = "<green>{time:HH:mm:ss}</green> | <cyan>{level}</cyan> | {message}"
    )
   

def run_cmd(cmd, cwd=None, shell=False):
    logger.debug(f"Running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    result = subprocess.run(cmd, cwd=cwd, shell=shell, check=True,
                            stdout=sys.stdout, stderr=sys.stderr)


def prepare_protein(input_pdb: str, output_pdb: str = "prot_p4a.pdb") -> str:
    logger.info("Cleaning protein with pdb4amber...")
    run_cmd([
        "pdb4amber",
        "-i", input_pdb,
        "-o", output_pdb,
        "--reduce"
    ])
    return output_pdb


def extract_ligand(complex_pdb: str, resn: str) -> str:
    output_pdb = f"{resn}_final.pdb"
    logger.info(f"Extracting ligand {resn} from complex...")
    run_cmd(f"awk '$4==\"{resn}\" {{print $0}}' {complex.pdb} > {output_pdb}", shell=True)
    return output_pdb
    

def prepare_ligand(ligand_pdb: str, charge: int = 0, resn: str = "LIG"):
    h_pdb = f"{resn}_H.pdb"
    logger.info("Protonating ligand...")
    run_cmd(["obabel", "-ipdb", ligand_pdb, "-opdb", "-O", "h_pdb", "-h"])

    sybl_mol2 = f"{resn}_sybl.mol2"
    logger.info("Running antechamber (SYBYL cleanup)...")
    run_cmd([ 
        "antechamber", "-i", h_pdb, "-fi", "pdb", "-o", sybl_mol2,
        "-fo", "mol2", "-j", "5", "-at", "sybyl", "-dr", "no", "-s", "2"
    ])
    
    final_mol2 = f"{resn}.mol2"
    logger.info("Running antechamber (GAFF2 + charges)...")
    run_cmd([
        "antechamber", "-i", sybl_mol2, "-fi", "mol2", "-o", final.mol2, "-fo", "mol2",
        "-at", "gaff2", "-c", "bcc", "-s", "2", "-rn", resn, "-nc", str(charge)
    ])

    frcmod = f"{resn}.frcmod"
    logger.info("Running parmchk2...")
    run_cmd([
        "parmchk2", "-i", final_mol2, "-f", "mol2", "-o", frcmod
    ])
    
    return final_mol2, frcmod


def write_tleap(protein_pdb: str, lig_mol2: str, lig_frcmod: str, resn: str, output: str = "tleap.in") -> str:
    tleap_script = f"""
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

loadamberparams {lig_frcmod}
protein = loadPDB {protein_pdb}
{resn} = loadmol2 {lig_mol2}

system = combine {{{resn} protein}}
savepdb system system.dry.pdb

check system

solvateBox system TIP3PBOX 10 iso
addions2 system Cl- 0
addions2 system Na+ 0

savePDB system system.pdb
saveAmberParm system system.prmtop system.inpcrd
quit
"""
    with open(output, "w") as f:
        f.write(tleap_script)
    return output


def minimize_and_convert_to_gromacs(prmtop_file: str, inpcrd_file: str, gmx_gro: str, gmx_top: str, pdb_out: str = None):
    logger.info("Loading Amber file for energy minimisation")
    inpcrd = AmberInpcrdFile(inpcrd_file)
    prmtop = AmberPrmtopFile(prmtop_file, periodicBoxVectors=inpcrd.boxVectors)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds
    )

    integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
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

    logger.info(f"Writing Gromacs files to {gmx_prefix}.gro and {gmx_prefix}.top")
    structure = pmd.openmm.load_topology(simulation.topology, system, positions)
    structure.save(gmx_top, overwrite=True)
    structure.save(gmx_gro, overwrite=True)

    logger.success(f"System prepared and Gromacs files written.")



def main():
    args = parse_args()
    configure_logging(args.verbose)

    complex = str(Path(args.input_pdb).resolve())
    
    logger.info("Step 1: Protein preparation...")
    protein_pdb = prepare_protein(complex)
    
    logger.info("Step 2: Ligand preparation...")
    ligand_pdb = extract_ligand(complex, resn=args.ligand_resname)
    lig_mol2, lig_frcmod = prepare_ligand(ligand_pdb, resn=args.ligand_resname, charge=args.ligand_charge)

    logger.info("Step 3: Building system with tleap...")
    tleap_in = write_tleap(protein_pdb, lig_mol2, lig_frcmod, resn=args.ligand_resname)
    run_cmd(["tleap", "-f", tleap_in])

    if args.do_minimise:
        logger.info("Step 4: Minimizing and converting to GROMACS...")
        minimize_and_convert_to_gromacs("system.prmtop", 
                                        "system.inpcrd",
                                        gmx_gro=args.gmx_gro, 
                                        gmx_top=args.gmx_top,
                                        min_pdb=args.min_pdb
                                       )

if __name__ == "__main__":
    main()
