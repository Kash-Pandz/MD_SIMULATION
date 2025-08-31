import argparse
import subprocess
import sys
from pathlib import Path
from loguru import logger


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare protein-ligand system.")
    parser.add_argument("--input_pdb", required=True, help="Input protein-ligand complex")
    parser.add_argument("--ph", type=float, default=7.0, help="pH for protein")
    parser.add_argument("--ligand_resname", default="LIG", help="Ligand resname")
    parser.add_argument("--ligand_charge", type=int, default=0, help="Net charge of ligand")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
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
    result = subprocess.run(cmd, cwd=cwd, check=True, shell=shell)


def protonate_protein(input_pdb: str, output_pdb: str = "prot_propka.pdb", pH: float = 7.0):
    logger.info(f"Protonating protein at pH {pH} using pdb2pqr (PropKa)...")
    run_cmd([
        "pdb2pqr",
        "--ff=amber",
        f"--with-ph={pH}",
        input_pdb,
        output_pdb
    ])
    return output_pdb


def prepare_protein(input_pdb: str, output_pdb: str = "prot_p4a.pdb") -> str:
    logger.info("Cleaning protein with pdb4amber...")
    run_cmd([
        "pdb4amber",
        "-i", input_pdb,
        "-o", output_pdb
    ])
    return output_pdb


def extract_ligand(complex_pdb: str, output_pdb: str = "LIG_final.pdb", resn: str = "LIG") -> str:
    logger.info(f"Extracting ligand {resn} from complex...")
    run_cmd(f"grep '^HETATM' {complex_pdb} | awk '{{if ($4==\"{resn}\") print $0}}' > {output_pdb}", shell=True)
    return output_pdb


def prepare_ligand(ligand_pdb: str, charge: int = 0, resn: str = "LIG"):
    logger.info("Protonating ligand...")
    h_pdb = f"{resid}_H.pdb"
    run_cmd(["obabel", "-ipdb", ligand_pdb, "-opdb", "-O", "h_pdb", "-h"])

    logger.info("Running antechamber (SYBYL cleanup)...")
    sybl_mol2 = f"{resid}_sybl.mol2"
    run_cmd([ 
        "antechamber", "-i", h_pdb, "-fi", "pdb", "-o", sybl_mol2,
        "-fo", "mol2", "-j", "5", "-at", "sybyl", "-dr", "no", "-s", "2"
    ])
    
    logger.info("Running antechamber (GAFF2 + charges)...")
    final_mol2 = f"{resid}.mol2"
    run_cmd([
        "antechamber", "-i", sybl_mol2, "-fi", "mol2", "-o", final.mol2, "-fo", "mol2",
        "-at", "gaff2", "-c", "bcc", "-s", "2", "-rn", resn, "-nc", str(charge)
    ])
    
    logger.info("Running parmchk2...")
    frcmod = f"{resid}.frcmod"
    run_cmd([
        "parmchk2", "-i", final.mol2, "-f", "mol2", "-o", frcmod
    ])
    
    return final_mol2, frcmod


def write_tleap(protein_pdb: str, lig_mol2: str, lig_frcmod: str, resn: str = "LIG", output: str = "tleap.in") -> str:
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



def main():
    args = parse_args()
    configure_logging(args.verbose)

    complex_file = str(Path(args.input_pdb).resolve())
    
    logger.info("Step 1: Protein preparation...")
    propka_pdb = protonate_protein(complex_file, pH=args.ph)
    protein_pdb = prepare_protein(propka_pdb)
    
    logger.info("Step 2: Ligand preparation...")
    ligand_pdb = extract_ligand(complex_file, resn=args.ligand_resname)
    lig_mol2, lig_frcmod = prepare_ligand(ligand_pdb, resn=args.ligand_resname, charge=args.ligand_charge)

    logger.info("Step 3: Building system with tleap...")
    tleap_in = write_tleap(protein_pdb, lig_mol2, lig_frcmod, resn=args.ligand_resname)
    run_cmd(["tleap", "-f", tleap_in])


if __name__ == "__main__":
    main()
