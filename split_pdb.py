import argparse
import MDAnalysis as mda


def split_pdb(pdb_file, ligand_name, protein_file="protein.pdb", ligand_file="ligand.pdb"):
    """Splits input PDB file into protein and ligand files."""

    # Load PDB file into MDAnalysis universe
    universe = mda.Universe(pdb_file)

    # Select protein and ligand atoms
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms(f"resname {ligand_name}")

    # Write protein and ligand pdb files
    protein.write(protein_file)
    ligand.write(ligand_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split PDB into protein and ligand files.")
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--ligand_name", required=True, help="Ligand residue name")
    parser.add_argument("--protein_out", default="protein.pdb", help="Output protein PDB file")
    parser.add_argument("--ligand_out", default="ligand.pdb", help="Output ligand PDB file")

    args = parser.parse_args()

    split_pdb(args.input, args.ligand_name, args.protein_out, args.ligand_out)
