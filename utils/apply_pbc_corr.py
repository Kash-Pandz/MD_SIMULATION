import argparse
import MDAnalysis as mda
import MDAnalysis.transformations as trans


def fix_pbc(
    topology: str,
    trajectory: str,
    output: str,
    selection: str = "protein",  # (protein or resname LIG)
    water_sel: str = "resname HOH WAT SOL"
) -> None:
    """
    Apply Periodic Boundary Condition corrections to a MD trajectory.
    """
    universe = mda.Universe(topology, trajectory)
    group = universe.select_atoms(selection)
    water = universe.select_atoms(water_sel)

    workflow = [
        trans.unwrap(universe.atoms),
        trans.center_in_box(group, center="geometry"),
        trans.wrap(water, compound="residues"),
    ]
    universe.trajectory.add_transformations(*workflow)

    print(f"Writing: {output}")
    with mda.Writer(output, n_atoms=universe.atoms.n_atoms) as writer:
        for ts in universe.trajectory:
            writer.write(universe)
    print(f"PBC fix for {trajectory} complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix PBC artifacts in MD trajectories")
    parser.add_argument("-s", "--topology", required=True, help="Topology file (.pdb, .prmtop, .gro, etc.)")
    parser.add_argument("-f", "--trajectory", required=True, help="Trajectory file (.dcd, .xtc, .nc, etc.)")
    parser.add_argument("-o", "--output", required=True, help="Output trajectory file")
    parser.add_argument("--selection", default="protein", help="Atom selection to center on (default: protein)")
    parser.add_argument("--water", default="resname HOH WAT SOL", help="Water selection string")
    args = parser.parse_args()

    fix_pbc(args.topology, args.trajectory, args.output, args.selection, args.water)
