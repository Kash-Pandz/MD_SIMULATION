import subprocess
import psi4
import resp


def calculate_resp_charges(ligand_pdb):
    """Geometry Optimisation at B3LYP/6-31G* then compute RESP charges at HF/6-31G*."""
    psi4.set_memory("4 GB")

    # Convert PDB to XYZ
    xyz_file = "ligand_tmp.xyz"
    subprocess.run(["obabel", "-ipdb", ligand_pdb, "-oxyz", "-O", xyz_file], check=True)

    # Load XYZ into Psi4
    with open(xyz_file, "r") as f:
        xyz_str = f.read()
    psi_mol = psi4.geometry(xyz_str)

    # Geometry optimisation (B3LYP/6-31G*)
    logger.info("Optimizing geometry (B3LYP/6-31G*)...")
    psi4.optimize("B3LYP/6-31G*", molecule=psi_mol)

    # ESP calculation (HF/6-31G*)
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

    os.remove(xyz_file)
    return resp_result[1]


def write_mol2_with_resp(ligand_pdb, charges, resn="LIG"):
    """Outputs a ligand mol2 file with RESP charges."""
    mol2_file = f"{resn}_RESP.mol2"
    # Convert PDB file to mol2
    subprocess.run(["obabel", "-ipdb", ligand_pdb, "-omol2", "-O", mol2_file], check=True)

    # Add RESP charges
    with open(mol2_file, "r") as f:
        lines = f.readlines()

    new_lines = []
    atom_idx = 0
    in_atoms = False
    for line in lines:
        if line.strip().startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue
        if line.strip().startswith("@<TRIPOS>BOND"):
            in_atoms = False
            new_lines.append(line)
            continue
        if in_atoms and len(line.split()) >= 9:
            parts = line.split()
            parts[8] = f"{charges[atom_idx]:.6f}"
            new_lines.append(" ".join(parts) + "\n")
            atom_idx += 1
        else:
            new_lines.append(line)

    with open(mol2_file, "w") as f:
        f.writelines(new_lines)

    logger.info(f"Saved ligand with RESP charges: {mol2_file}")
    return mol2_file
