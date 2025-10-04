import os
import subprocess
from loguru import logger
import psi4
import resp


def calculate_resp_charges(ligand_pdb: str, resn: str = "LIG", mol_name="LIGNAME"):
    """
    RESP charges at HF/6-31G* with optional geometry optimisation at B3LYP/6-31G*. 
    Writes .mol2 file with RESP charges.
    """
    # Add/Check hydrogens using OpenBabel
    #lig_h = "ligand_H.pdb"
    #subprocess.run(["obabel", ligand_pdb, "-O", lig_h, "-h"], check=True)

    # Convert PDB to XYZ
    xyz_file = "ligand_tmp.xyz"
    subprocess.run(["obabel", "-ipdb", ligand_pdb, "-oxyz", "-O", xyz_file], check=True)

    with open(xyz_file, "r") as f:
        xyz_str = f.read()

    # Psi4 setup
    psi4.set_memory("4 GB")
    psi4.set_num_threads(2)
    psi_mol = psi4.geometry(xyz_str)

    # Optional geometry optimisation at B3LYP/6-31G*
    logger.info("Optimizing geometry (B3LYP/6-31G*)...")
    psi4.optimize("B3LYP/6-31G*", molecule=psi_mol)

    # ESP calculation at HF/6-31G*
    logger.info("Computing ESP (HF/6-31G(d))...")
    params = {
        'METHOD_ESP': 'HF',
        'BASIS_ESP': '6-31G*',
        'RESP_A': 0.0005,
        'RESP_B': 0.1,
        'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0],
        'VDW_POINT_DENSITY': 1,
    }
    resp_result = resp.resp([psi_mol], params)
    charges = resp_result[1]
    
    os.remove(xyz_file)

    # Load PDB into ParmEd
    parm = pmd.load_file(ligand_pdb)

    # Assign the RESP charges
    logger.info(f"Assigning {len(charges)} RESP charges tp {len(parm.atoms)} atoms...")
    assert len(charges) == len(parm.atoms), f"Mismatch between ({len(charges)}) and ({len(parm.atoms)})!"

    for atom, charge in zip(parm.atoms, charges):
        atom.charge = charge

    # Write .mol2 file with RESP charges
    mol2_file = f"{resn}_RESP.mol2"
    parm.save(mol2_file, format="mol2")
    logger.success(f"Saved mol2 file with RESP charges")
    
    return mol2_file
