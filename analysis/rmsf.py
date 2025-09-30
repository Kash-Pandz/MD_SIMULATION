import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple, Optional


def calculate_rmsf(
    topology: str,
    trajectory: str,
    atom_sel: str = "protein and name CA",
    ref_frame: int = 0,
    align_in_memory: bool = True,
    #pdb_out: Optional[str] = None
) -> pd.DataFrame:
    """
    Calculate the Root Mean Square Fluctuation of CA atoms over a trajectory.

    Args:
        topology: 
            Topology file (e.g. .pdb, .gro, .tpr, .psf & .prmtop).
        trajectory:
            Trajectory file (e.g. .xtc, .trr, .dcd & .nc).
        atom_sel:
            MDAnalysis selection string for atoms to calculate RMSF with.
            Default is "protein and name CA". Selects the CA atoms of the protein.
        ref_frame:
            Frame index to which the average structure is initially aligned.
            Default is 0 (first frame)
        align_in_memory:
            Whether to keep alignment in memory. If False, this will write a new trajectory.
            Default is True.

    Returns:
        pd.DataFrame:
            Pandas dataframe with residue and rmsf value in angstroms.
    """
    universe = mda.Universe(topology, trajectory)
    sel = universe.select_atoms(atom_sel)

    if len(sel) == 0:
        raise ValueError(f"No atoms selected with selection '{atom_sel}'")

    # Get average structure 
    avg = align.AverageStructure(
        universe, universe, select=atom_sel, ref_frame=ref_frame
    ).run()
    ref = avg.results.universe

    # Align trajectory to average structure
    align.AlignTraj(
        universe, ref, select=atom_sel, in_memory=align_in_memory
    ).run()

    rmsf = rms.RMSF(sel).run()
    rmsf_values = rmsf.results.rmsf

    # Store in pandas dataframe
    df = pd.DataFrame({
        "resid": [atom.resid for atom in sel]
        "resn": [atom.resname for atom in sel]
        "rmsf": rmsf_values
    })

    # Write PDB with RMSF values in B-Factor column
    #if pdb_out is not None:
        #universe.add_TopologyAtrr("tempfactors")
        #for atom, val in zip(sel, rmsf_values):
            #atom.tempfactor = val
        #universe.atoms.write(pdb_out)

    return df
  

def plot_rmsf(
    df: pd.DataFrame,
    xlabel: str = "Residue",
    ylabel: str = "RMSF (Å)",
    figsize: Tuple[float, float] = (5, 5),
    dpi: int = 300
) -> None:
    """
    Plot RMSF data.

    Args:
        df: 
          Pandas Dataframe with RMSF results.
        xlabel: 
          X axis label. Default is "Residue".
        ylabel:
          Y axis label. Deafult is "RMSF (Å)".
        figsize:
          Figure size in inches (width, height). Default is (5,5).
        dpi: 
          Resolution of the figure. Default is 300.
    """
    plt.figure(figsize=figsize, dpi=dpi)
    plt.plot(df["resid"], df["rmsf"], "-o", markersize=3)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.show()
