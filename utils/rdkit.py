from rdkit import Chem 
from rdkit.Chem import AllChem


def smi_to_xyz(mol, conf_id=0):
    """Write ligand smi to xyz."""
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformers. Generate ligand coordinates.")

    conf = mol.GetConformer(conf_id)
    num_atoms = mol.GetNumAtoms()

    xyz_lines = [f"{num_atoms}"]
    for i in range(num_atoms):
        pos = conf.GetAtomsPosition(i)
        sym = mol.GetAtomwithIdx(i).GetSymbol()
        xyz_lines.append(f"{sym:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    
    return "\n".join(xyz_lines)



    
