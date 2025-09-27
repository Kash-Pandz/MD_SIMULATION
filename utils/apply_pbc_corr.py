import argparse
import MDAnalysis as mda
import MDAnalysis.transformations as trans

def fix_pbc(
    topology: str,
    trajectory: str,
    output: str,
    selection: str = "protein"  #"protein or resname LIG"
) -> None:
    """
    Apply Periodic Boundary Condition corrections to a MD trajectory.

    Args:
        topology: 
            Topology file (e.g. .gro, .tpr, .psf, .pdb & .prmtop).
        trajectory: 
            Trajectory file (e.g. .xtc, .dcd & .nc).
        output: 
            Output trajectory filename.
        selection: 
            Atom selection for PBC correction (e.g. protein & protein+ligand)
    """
    
    universe = mda.Universe(top_file, traj_file)

    group = universe.select_atoms(selection)
    water = universe.select_atoms("resname SOL")

    workflow = [
        trans.unwrap(universe.atoms),
        trans.center_in_box(group, center="geometry"),
        trans.wrap(water, compound="residues"),
    ]
    u.trajectory.add_transformations(*workflow)

    print(f"Writing: {output}")
    with mda.Writer(output, n_atoms=universe.atoms.n_atoms) as writer:
        for ts in universe.trajectory:
            writer.write(universe)

    print(f"PBC fix for {trajectory} complete!")


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply PBC correction in MD trajectories.")
    parser.add_argument("topology", help="Topology file (e.g. .gro, .tpr, .psf, .pdb & .prmtop)")
    parser.add_argument("trajectories", nargs="+",
                        help="Trajectory files (e.g. .xtc, .dcd & .nc)")
    parser.add_argument(
        "--selection",
        default="protein",
        help=('Atom selection (default: "protein"; '
              'e.g., "protein or resname LIG")')
    )

    args = parser.parse_args()

    for traj in args.trajectories:
        if traj.endswith(".xtc"):
            outname = traj.replace(".xtc", "_fixed.xtc")
        elif traj.endswith(".dcd"):
            outname = traj.replace(".dcd", "_fixed.dcd")
        elif traj.endswith(".nc"):
            outname = traj.replace(".nc", "_fixed.nc")
        else:
            outname = traj + "_fixed"

        fix_pbc(args.topology, traj, outname, args.selection)


if __name__ == "__main__":
    main()
