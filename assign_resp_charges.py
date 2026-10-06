import os
import subprocess
import warnings
from loguru import logger
import psi4
import resp
import parmed as pmd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligand_from_smiles(smiles: str, pdb_file: str, output_prefix: str = "LIG"):
    """
    Fix bond orders on a PDB using SMILES template, add H.
    Heavy atoms stay in their docked positions.
    Returns (pdb_with_H, xyz_with_H).
    """
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")

    mol_pdb = Chem.MolFromPDBFile(pdb_file, removeHs=True, sanitize=False)
    if mol_pdb is None:
        raise ValueError(f"RDKit could not read PDB: {pdb_file}")

    mol_fixed = AllChem.AssignBondOrdersFromTemplate(template, mol_pdb)
    mol_h = Chem.AddHs(mol_fixed, addCoords=True)

    # Normalise residue name on ALL atoms to the short output_prefix (e.g. "LIG")
    # The input PDB may have long residue names (e.g. "LIG3C") that exceed the
    # 3-4 character PDB limit and cause antechamber to fail.
    # Also patch H atoms that AddHs created with missing/default monomer info.
    resn_short = output_prefix[:3]  # PDB residue name: max 3 chars
    for atom in mol_h.GetAtoms():
        mi = atom.GetMonomerInfo()
        if atom.GetAtomicNum() == 1:
            # H atom: ensure it has monomer info, inheriting from its parent heavy atom
            parent = atom.GetNeighbors()[0]
            parent_mi = parent.GetMonomerInfo()
            if mi is None:
                mi = Chem.AtomPDBResidueInfo()
                mi.SetName(f" H{atom.GetIdx():<2d}")
                mi.SetIsHeteroAtom(True)
                atom.SetMonomerInfo(mi)
            if parent_mi is not None:
                mi.SetResidueNumber(parent_mi.GetResidueNumber())
                mi.SetChainId(parent_mi.GetChainId())
        # Set the normalized residue name on every atom (heavy and H)
        if mi is not None:
            mi.SetResidueName(resn_short)

    # Sanity check to ensure heavy atoms have not moved
    conf_pdb = mol_pdb.GetConformer()
    conf_h = mol_h.GetConformer()

    n_heavy = mol_pdb.GetNumAtoms()
    coords_pdb = np.array(
        [conf_pdb.GetAtomPosition(i) for i in range(n_heavy)]
    )
    heavy_indices = [
        i for i, a in enumerate(mol_h.GetAtoms()) if a.GetAtomicNum() > 1
    ]
    coords_h = np.array(
        [conf_h.GetAtomPosition(i) for i in heavy_indices]
    )

    if len(coords_pdb) != len(coords_h):
        raise RuntimeError(
            f"Heavy-atom count mismatch: PDB has {len(coords_pdb)}, "
            f"mol_h has {len(coords_h)}"
        )

    rmsd = float(np.sqrt(np.mean(np.sum((coords_pdb - coords_h) ** 2, axis=1))))
    if rmsd > 0.001:
        raise RuntimeError(f"Heavy atoms moved during addition of Hs! RMSD = {rmsd:.3f} Å")
    logger.info(f"Heavy atom RMSD after H addition: {rmsd:.3f} Å")

    # Write outputs — PDB for downstream use, XYZ for Psi4, SDF for mol2 conversion
    pdb_out = f"{output_prefix}_H.pdb"
    xyz_out = f"{output_prefix}_H.xyz"
    sdf_out = f"{output_prefix}_H.sdf"
    names_out = f"{output_prefix}_atomnames.txt"
    Chem.MolToPDBFile(mol_h, pdb_out)
    Chem.MolToXYZFile(mol_h, xyz_out)
    Chem.MolToMolFile(mol_h, sdf_out)  # SDF preserves bond orders

    # Save atom names from the PDB (heavy atoms keep original names, H atoms
    # were named during patching above). These are needed so the final mol2
    # atom names match the complex PDB for tleap.
    atom_names = []
    for atom in mol_h.GetAtoms():
        mi = atom.GetMonomerInfo()
        if mi is not None:
            atom_names.append(mi.GetName().strip())
        else:
            atom_names.append(f"{atom.GetSymbol()}{atom.GetIdx()}")
    with open(names_out, "w") as f:
        for name in atom_names:
            f.write(name + "\n")

    n_total = mol_h.GetNumAtoms()
    n_h = n_total - mol_h.GetNumHeavyAtoms()
    logger.info(f"Wrote {n_total} atoms ({mol_h.GetNumHeavyAtoms()} heavy + {n_h} H) → {pdb_out}")

    return pdb_out, xyz_out, sdf_out, names_out


def calculate_resp_charges(ligand_pdb: str, resn: str = "LIG", charge: int = 0,
                           multiplicity: int = 1, smiles: str = None):
    """
    RESP charges at HF/6-31G*.
    Writes .mol2 file with GAFF2 atom types and RESP charges.
    Pipeline: RDKit SDF → obabel mol2 → antechamber GAFF2 retype → charge patch.
    """
    # Prepare structure with hydrogens
    if smiles:
        ligand_h_pdb, xyz_file, sdf_file, names_file = prepare_ligand_from_smiles(smiles, ligand_pdb, resn)
    else:
        ligand_h_pdb = ligand_pdb
        xyz_file = "ligand_tmp.xyz"
        sdf_file = None
        names_file = None
        subprocess.run(["obabel", "-ipdb", ligand_pdb, "-oxyz", "-O", xyz_file], check=True)

    # Read XYZ and prepend charge/multiplicity for Psi4
    with open(xyz_file, "r") as f:
        lines = f.readlines()
    atom_lines = [l.strip() for l in lines[2:] if len(l.strip().split()) >= 4]
    xyz_str = f"{charge} {multiplicity}\n" + "\n".join(atom_lines)

    # Psi4 setup
    psi4.core.clean()             # reset state from any previous run
    psi4.core.clean_timers()      # clear stale timers
    psi4.core.clean_options()     # reset options
    psi4.set_memory("2 GB")
    psi4.set_num_threads(2)
    psi4.core.set_output_file(f"{resn}_psi4.log", False)
    psi_mol = psi4.geometry(xyz_str)
    psi_mol.update_geometry()
    logger.info(f"Molecule: {psi_mol.natom()} atoms, charge={charge}, mult={multiplicity}")

    psi4.set_options({
        "basis": "6-31G*",
        "reference": "rhf" if multiplicity == 1 else "uhf",
        "scf_type": "df",
    })

    # Two-stage RESP fitting at HF/6-31G* (Bayly/Cieplak/Cornell protocol)
    # Stage 1: weak restraint (a=0.0005) on all atoms, no symmetry enforcement
    # Stage 2: stronger restraint (a=0.001), freeze polar/well-determined atoms,
    #          refit methyl/methylene H with equivalence constraints enforced
    logger.info("Computing ESP and RESP charges (HF/6-31G*)...")

    # --- Stage 1 ---
    logger.info("RESP Stage 1: weak restraint (a=0.0005), all atoms free...")
    options = {
        'METHOD_ESP': 'HF',
        'BASIS_ESP': '6-31G*',
        'RESP_A': 0.0005,
        'RESP_B': 0.1,
        'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0],
        'VDW_POINT_DENSITY': 1,
        'IHFREE': True,
        'TOLER': 1e-5,
        'MAX_IT': 25,
    }

    charges1 = resp.resp([psi_mol], options)
    logger.info(f"Stage 1 charges — sum: {np.sum(charges1[1]):.4f}, "
                f"range: [{np.min(charges1[1]):.4f}, {np.max(charges1[1]):.4f}]")

    # --- Stage 2 ---
    # Use resp.set_stage2_constraint() to automatically freeze non-methyl/methylene
    # atoms and add equivalence constraints for equivalent H groups
    logger.info("RESP Stage 2: stronger restraint (a=0.001), "
                "freezing polar atoms, equivalencing methyl/methylene H...")
    options['RESP_A'] = 0.001
    resp.set_stage2_constraint(psi_mol, charges1[1], options)

    # Reuse the saved grid and ESP files from stage 1 (avoids recomputing ESP)
    options['grid'] = ['1_%s_grid.dat' % psi_mol.name()]
    options['esp'] = ['1_%s_grid_esp.dat' % psi_mol.name()]
    psi_mol.set_name('stage2')

    charges2 = resp.resp([psi_mol], options)
    charges = charges2[1]
    logger.info(f"Stage 2 charges — sum: {np.sum(charges):.4f}, "
                f"range: [{np.min(charges):.4f}, {np.max(charges):.4f}]")

    # Clean up intermediate xyz file
    try:
        os.remove(xyz_file)
    except OSError:
        pass

    # Save plain-text charges
    charge_file = f"{resn}_resp_charges.txt"
    np.savetxt(charge_file, charges, fmt="%.6f")

    # Generate mol2 with correct bond orders and SYBYL atom types
    # Strategy: use the SDF from RDKit (has correct bond orders) and convert
    # via obabel to mol2 (assigns SYBYL types), then patch in RESP charges.
    mol2_file = f"{resn}_resp.mol2"
    mol2_tmp = f"{resn}_tmp.mol2"

    if sdf_file and os.path.isfile(sdf_file):
        # SDF preserves bond orders — obabel will assign SYBYL types correctly
        result = subprocess.run(
            ["obabel", "-isdf", sdf_file, "-omol2", "-O", mol2_tmp],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"obabel STDERR:\n{result.stderr}")
            raise RuntimeError("obabel SDF to mol2 conversion failed.")
    else:
        # Fallback: convert PDB via obabel (less reliable for bond orders)
        result = subprocess.run(
            ["obabel", "-ipdb", ligand_h_pdb, "-omol2", "-O", mol2_tmp],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"obabel STDERR:\n{result.stderr}")
            raise RuntimeError("obabel PDB to mol2 conversion failed.")

    # Retype to GAFF2 via antechamber 
    mol2_retyped = f"{resn}_gaff2_tmp.mol2"
    retype_result = subprocess.run([
        "antechamber",
        "-i", mol2_tmp,
        "-fi", "mol2",      # mol2 input - bonds already correct
        "-o", mol2_retyped,
        "-fo", "mol2",
        "-at", "gaff2",
        "-j", "4",          # keep existing bond types, run atom type assignment
        "-pf", "y",
        "-dr", "no",        # disable acdoctor for cleaner run
    ], capture_output=True, text=True)

    if retype_result.returncode != 0:
        logger.error(f"antechamber retype STDOUT:\n{retype_result.stdout}")
        logger.error(f"antechamber retype STDERR:\n{retype_result.stderr}")
        raise RuntimeError("antechamber GAFF2 retyping failed. Check logs above.")

    try:
        os.remove(mol2_tmp)
    except OSError:
        pass
    mol2_tmp = mol2_retyped
    logger.info("Atom types converted to GAFF2 via antechamber")

    # Load original atom names from PDB (for tleap compatibility)
    atom_names = None
    if names_file and os.path.isfile(names_file):
        with open(names_file, "r") as f:
            atom_names = [line.strip() for line in f if line.strip()]

    # Patch RESP charges, residue name, and atom names into the mol2 file
    _patch_mol2_charges(mol2_tmp, charges, resn, mol2_file, atom_names=atom_names)

    # Clean up temp file
    try:
        os.remove(mol2_tmp)
    except OSError:
        pass

    logger.success(f"Saved mol2 file with RESP charges: {mol2_file}")

    # Generate frcmod for missing GAFF2 parameters
    frcmod_file = f"{resn}.frcmod"
    frcmod_result = subprocess.run([
        "parmchk2",
        "-i", mol2_file,
        "-f", "mol2",
        "-o", frcmod_file,
        "-s", "gaff2",
    ], capture_output=True, text=True)

    if frcmod_result.returncode != 0:
        logger.error(f"parmchk2 STDERR:\n{frcmod_result.stderr}")
        raise RuntimeError("parmchk2 failed. Check logs above.")
    logger.success(f"Saved frcmod file: {frcmod_file}")

    return mol2_file, frcmod_file


def patch_mol2_charges(mol2_in: str, charges: np.ndarray, resn: str, mol2_out: str,
                       atom_names: list = None, net_charge: int = None,
                       charge_tol: float = 1e-3):
    """
    Read mol2 file, replace the partial charghes with RESP charges and set the
    residue name (ATOM + SUBSTRUCTURE blocks and the MOLECULE name line).
    Optionally replace atom names to match the complex PDB (for tleap).
    Write the patched file to mol2_out.
    """
    struct = pmd.load_file(mol2_in, structure=True)
    atoms = struct.atoms
 
    charges = np.asarray(charges, dtype=float)
    if charges.ndim != 1:
        raise ValueError(f"charges must be 1-D, got shape {charges.shape}")
    if len(atoms) != len(charges):
        raise RuntimeError(
            f"Charge count mismatch: {mol2_in} has {len(atoms)} atoms, "
            f"but {len(charges)} RESP charges were supplied."
        )
 
    if net_charge is not None:
        total = charges.sum()
        if abs(total - net_charge) > charge_tol:
            raise RuntimeError(
                f"RESP charges sum to {total:.6f}, expected {net_charge:+d} "
                f"(tolerance {charge_tol})."
            )
 
    if len(resn) > 4:
        raise ValueError(f"Residue name {resn!r} is longer than 4 characters.")
 
    if atom_names is not None:
        if len(atom_names) != len(atoms):
            raise RuntimeError(
                f"Atom name count mismatch: {mol2_in} has {len(atoms)} atoms, "
                f"but {len(atom_names)} names were supplied."
            )
        if len(set(atom_names)) != len(atom_names):
            dupes = sorted({n for n in atom_names if atom_names.count(n) > 1})
            raise ValueError(f"Atom names must be unique within a residue; "
                             f"duplicates: {dupes}")
        too_long = [n for n in atom_names if len(n) > 4]
        if too_long:
            raise ValueError(f"Atom names must be <= 4 characters; offenders: "
                             f"{too_long}")
 
    for i, atom in enumerate(atoms):
        atom.charge = float(charges[i])
        if atom_names is not None:
            atom.name = atom_names[i]
 
    for res in struct.residues:
        res.name = resn
 
    # Keeps the MOLECULE name line in sync with the residue name; some tools
    # (antechamber, parmchk2) read that rather than SUBSTRUCTURE.
    struct.name = resn
 
    with warnings.catch_warnings():
        warnings.simplefilter("error", pmd.exceptions.ParmedWarning)
        pmd.formats.Mol2File.write(struct, mol2_out)
