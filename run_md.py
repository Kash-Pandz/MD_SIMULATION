#!/usr/bin/env python3
"""
MD Simulation with OpenMM — Standalone Script (Amber/tleap input)
=================================================================

Protocol (run independently for each temperature):
1. Load solvated system from tleap (Amber prmtop + inpcrd)
2. Two-stage energy minimisation (restrained → unrestrained)
3. NVT heating (500 ps, heavy-atom restraints)
4. NPT equilibration (5 ns, gradual backbone restraint release)
5. Production NPT (unrestrained, 3 replicas × 100 ns)

Prerequisites:
    Run tleap to generate system.prmtop and system.inpcrd
    (solvated, neutralised, with correct protonation states)

Usage:
    python run_md.py
    nohup python run_md.py > md.log 2>&1 &
"""

import openmm as mm
from openmm import app, unit
from openmm.app import (
    Simulation,
    AmberPrmtopFile,
    AmberInpcrdFile,
    PDBFile,
    StateDataReporter,
    DCDReporter,
)
import sys
import os
import random
import time

# ============================================================
# USER CONFIGURATION
# ============================================================

PRMTOP_FILE   = "system_APO.prmtop"            # Amber topology from tleap
INPCRD_FILE   = "system_APO.inpcrd"            # Amber coordinates from tleap
#LIGAND_RESNAME = "SUB"                     # Ligand residue name in topology
TEMPERATURES  = [308.15, 338.15]           # K — both temperatures run sequentially

# Note: restraint force constants in kcal/mol/Å²
# (converted to OpenMM units internally)

# Minimisation — two stage
MIN_RESTRAINED_K = 10.0                # kcal/mol/Å², heavy atom restraints
MIN_FREE_K       = 0.0                 # unrestrained

# NVT heating (500 ps, strong heavy-atom restraints)
NVT_STEPS        = 250_000             # 500 ps at 2 fs
NVT_RESTRAINT_K  = 10.0                # kcal/mol/Å²
NVT_RESTRAINT_MODE = "heavy"           # all non-H atoms

# NPT equilibration — gradual backbone restraint release
# Each stage: (steps, k in kcal/mol/Å², mode)
NPT_STAGES = [
    (500_000,  5.0, "backbone"),        # 1 ns
    (500_000,  2.0, "backbone"),        # 1 ns
    (500_000,  1.0, "backbone"),        # 1 ns
    (500_000,  0.1, "backbone"),        # 1 ns
    (500_000,  0.0, "backbone"),        # 1 ns unrestrained validation
]

# Production replicas
N_REPLICAS        = 3
PROD_STEPS        = 50_000_000         # 100 ns at 2 fs per replica
TIMESTEP          = 0.002              # ps (2 fs)

# Reporting
LOG_FREQ         = 5_000               # 10 ps
TRJ_FREQ         = 50_000              # 100 ps
STDOUT_FREQ      = 50_000              # 100 ps


# ============================================================
# PLATFORM DETECTION
# ============================================================

def get_platform():
    """Detect CUDA if available, otherwise fall back to CPU."""
    try:
        mm.Platform.getPlatformByName("CUDA")
        platform_name = "CUDA"
        platform_props = {
            "Precision": "mixed",
            "DeviceIndex": "0",
        }
        print("Using CUDA platform with mixed precision.")
    except Exception:
        platform_name = "CPU"
        platform_props = {}
        print("CUDA not found. Using CPU platform.")
    return platform_name, platform_props


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def load_amber_system(prmtop_path, inpcrd_path):
    """Load Amber prmtop/inpcrd files from tleap."""
    print(f"  Loading topology: {prmtop_path}")
    print(f"  Loading coordinates: {inpcrd_path}")

    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    print(f"  System: {prmtop.topology.getNumAtoms()} atoms, "
          f"{prmtop.topology.getNumResidues()} residues")

    if inpcrd.boxVectors is not None:
        print(f"  Box vectors from inpcrd: {inpcrd.boxVectors}")
    else:
        print("  WARNING: No box vectors found in inpcrd file.")

    return prmtop, inpcrd


def build_system(prmtop):
    """Create an OpenMM System from Amber topology with PME and HBond constraints."""
    kwargs = dict(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )
    return prmtop.createSystem(**kwargs)


def set_integrator(temperature, friction=1 / unit.picosecond, timestep=0.002):
    """Create a Langevin Middle integrator."""
    return mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin, friction, timestep * unit.picoseconds
    )


def apply_restraint(system, topology, positions, k=0.0, chain_indices=None, mode="backbone",
                    exclude_resnames=None):
    """Apply positional restraints to selected protein atoms.

    mode: 'all', 'heavy' (non-H), or 'backbone' (CA, C, N, O)
    k: force constant in kcal/mol/Å²
    exclude_resnames: set of residue names to skip (e.g. ligand, ions)
    """
    if k <= 0.0:
        return None
    if exclude_resnames is None:
        exclude_resnames = set()

    restraint = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
    # Convert from kcal/mol/Å² to OpenMM internal units
    # 1 kcal/mol/Å² = 100 kcal/mol/nm² = 418.4 kJ/mol/nm²
    k_openmm = k * unit.kilocalorie_per_mole / unit.angstrom ** 2
    restraint.addGlobalParameter("k", k_openmm)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    chains = list(topology.chains())
    if chain_indices is None:
        chain_indices = range(len(chains))

    n_restrained = 0
    for chain_index in chain_indices:
        if chain_index >= len(chains):
            print(f"  Warning: Chain index {chain_index} out of range. Skipping.")
            continue
        for atom in chains[chain_index].atoms():
            if atom.residue.name in exclude_resnames:
                continue
            if mode == "heavy" and atom.element.symbol == "H":
                continue
            if mode == "backbone" and atom.name not in ("CA", "C", "N", "O"):
                continue
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
            n_restrained += 1

    if n_restrained > 0:
        system.addForce(restraint)
        print(f"  Restrained {n_restrained} protein atoms (mode={mode}, k={k:.2f} kcal/mol/Å²)")
    return restraint


def apply_ligand_restraint(system, topology, positions, ligand_resname, k=0.0):
    """Apply heavy-atom positional restraints to ligand residues only.

    k: force constant in kcal/mol/Å²
    """
    if k <= 0.0:
        return None

    restraint = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
    k_openmm = k * unit.kilocalorie_per_mole / unit.angstrom ** 2
    restraint.addGlobalParameter("k", k_openmm)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    n_restrained = 0
    for atom in topology.atoms():
        if atom.residue.name != ligand_resname:
            continue
        if atom.element.symbol == "H":
            continue
        pos = positions[atom.index].value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, [pos[0], pos[1], pos[2]])
        n_restrained += 1

    if n_restrained > 0:
        system.addForce(restraint)
        print(f"  Restrained {n_restrained} ligand heavy atoms ({ligand_resname}, k={k:.2f} kcal/mol/Å²)")
    return restraint


def create_simulation(topology, system, integrator, platform_name, platform_props, box_vectors=None):
    """Build a Simulation and set box vectors if available."""
    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = Simulation(
        topology, system, integrator,
        platform=platform, platformProperties=platform_props,
    )
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)
    return simulation


def load_checkpoint(simulation, checkpoint_file):
    """Load a checkpoint if it exists. Returns True if loaded."""
    if checkpoint_file and os.path.exists(checkpoint_file):
        with open(checkpoint_file, "rb") as f:
            simulation.context.loadCheckpoint(f.read())
        print(f"  Resumed from checkpoint: {checkpoint_file}")
        return True
    return False


def save_checkpoint(simulation, checkpoint_file):
    """Save a checkpoint."""
    if checkpoint_file:
        with open(checkpoint_file, "wb") as f:
            f.write(simulation.context.createCheckpoint())
        print(f"  Checkpoint saved: {checkpoint_file}")


# ============================================================
# MAIN SIMULATION PIPELINE
# ============================================================

def run_simulation(temperature, prmtop_file, inpcrd_file, platform_name, platform_props):
    """Run the full simulation pipeline for a single temperature."""

    temp_dir = f"T_{temperature:.0f}K"
    os.makedirs(temp_dir, exist_ok=True)
    start_dir = os.getcwd()
    os.chdir(temp_dir)

    print(f"\n{'='*70}")
    print(f"  TEMPERATURE: {temperature} K")
    print(f"  Working directory: {os.getcwd()}")
    print(f"{'='*70}")

    # ----------------------------------------------------------
    # Load Amber topology & coordinates from tleap
    # ----------------------------------------------------------
    print("\n--- Loading Amber prmtop/inpcrd ---")
    prmtop_path = os.path.join(start_dir, prmtop_file)
    inpcrd_path = os.path.join(start_dir, inpcrd_file)
    prmtop, inpcrd = load_amber_system(prmtop_path, inpcrd_path)
    topology = prmtop.topology
    positions = inpcrd.positions
    velocities = None
    box_vectors = inpcrd.boxVectors

    # ----------------------------------------------------------
    # 1. Two-stage energy minimisation
    # ----------------------------------------------------------
    print("\n" + "="*70)
    print("  STAGE 1: ENERGY MINIMISATION")
    print("="*70)

    min_stages = [
        ("min_restrained", MIN_RESTRAINED_K, "heavy"),
        ("min_free",       MIN_FREE_K,       "heavy"),
    ]

    for prefix, k, mode in min_stages:
        print(f"\n--- Minimisation: {prefix} (k={k:.1f} kcal/mol/Å², mode={mode}) ---")
        system = build_system(prmtop)
        if box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*box_vectors)
        apply_restraint(system, topology, positions, k=k, mode=mode, exclude_resnames=None)
        #apply_restraint(system, topology, positions, k=k, mode=mode, exclude_resnames={LIGAND_RESNAME})
        #apply_ligand_restraint(system, topology, positions, LIGAND_RESNAME, k=k)
        integrator = set_integrator(temperature)

        sim = create_simulation(topology, system, integrator, platform_name, platform_props,
                                box_vectors=box_vectors)
        sim.context.setPositions(positions)
        sim.minimizeEnergy()

        state = sim.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions()
        print(f"  Potential energy: {state.getPotentialEnergy()}")

        with open(f"{prefix}.pdb", "w") as f:
            PDBFile.writeFile(topology, positions, f)
        print(f"  Wrote {prefix}.pdb")

    velocities = None
    print("\nMinimisation complete.")

    # ----------------------------------------------------------
    # 2. NVT heating
    # ----------------------------------------------------------
    print("\n" + "="*70)
    print(f"  STAGE 2: NVT HEATING ({NVT_STEPS * TIMESTEP / 1000:.1f} ns)")
    print("="*70)

    nvt_ns = NVT_STEPS * TIMESTEP / 1000
    print(f"  Restraints: k={NVT_RESTRAINT_K} kcal/mol/Å², mode={NVT_RESTRAINT_MODE}")

    system = build_system(prmtop)
    if box_vectors is not None:
        system.setDefaultPeriodicBoxVectors(*box_vectors)
    apply_restraint(system, topology, positions, k=k, mode=mode, exclude_resnames=None)
    #apply_restraint(system, topology, positions, k=NVT_RESTRAINT_K, mode=NVT_RESTRAINT_MODE, exclude_resnames={LIGAND_RESNAME})
    #apply_ligand_restraint(system, topology, positions, LIGAND_RESNAME, k=NVT_RESTRAINT_K)
    integrator = set_integrator(temperature, timestep=TIMESTEP)

    sim = create_simulation(topology, system, integrator, platform_name, platform_props,
                            box_vectors=box_vectors)
    sim.context.setPositions(positions)
    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    sim.reporters.append(
        StateDataReporter(sys.stdout, STDOUT_FREQ,
                          step=True, temperature=True,
                          potentialEnergy=True, speed=True)
    )
    sim.reporters.append(
        StateDataReporter(f"NVT_{temperature:.0f}K.log", LOG_FREQ,
                          step=True, temperature=True,
                          potentialEnergy=True, kineticEnergy=True)
    )
    sim.reporters.append(DCDReporter(f"NVT_{temperature:.0f}K.dcd", TRJ_FREQ))

    sim.step(NVT_STEPS)

    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()
    print(f"\nNVT complete. PE: {state.getPotentialEnergy()}")

    # ----------------------------------------------------------
    # 3. NPT equilibration — gradual restraint release
    # ----------------------------------------------------------
    total_npt_steps = sum(s for s, _, _ in NPT_STAGES)
    total_npt_ns = total_npt_steps * TIMESTEP / 1000

    print("\n" + "="*70)
    print(f"  STAGE 3: NPT EQUILIBRATION ({total_npt_ns:.1f} ns, {len(NPT_STAGES)} stages)")
    print("="*70)
    print(f"  Restraint schedule: {[(k, m) for _, k, m in NPT_STAGES]}")

    for stage_idx, (stage_steps, stage_k, stage_mode) in enumerate(NPT_STAGES, 1):
        stage_ns = stage_steps * TIMESTEP / 1000
        print(f"\n--- NPT Stage {stage_idx}/{len(NPT_STAGES)}: "
              f"{stage_ns:.1f} ns, k={stage_k:.1f} kcal/mol/Å², mode={stage_mode} ---")

        system = build_system(prmtop)
        if box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*box_vectors)
        if stage_k > 0:
            apply_restraint(system, topology, positions, k=k, mode=mode, exclude_resnames=None)
            apply_restraint(system, topology, positions, k=stage_k, mode=stage_mode, exclude_resnames={LIGAND_RESNAME})
            apply_ligand_restraint(system, topology, positions, LIGAND_RESNAME, k=stage_k)
        system.addForce(mm.MonteCarloBarostat(1 * unit.bar, temperature * unit.kelvin))
        integrator = set_integrator(temperature, timestep=TIMESTEP)

        sim = create_simulation(topology, system, integrator, platform_name, platform_props,
                                box_vectors=box_vectors)
        sim.context.setPositions(positions)
        if velocities is not None:
            sim.context.setVelocities(velocities)
        else:
            sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)

        sim.reporters.append(
            StateDataReporter(sys.stdout, STDOUT_FREQ,
                              step=True, temperature=True,
                              potentialEnergy=True, volume=True,
                              density=True, speed=True)
        )
        sim.reporters.append(
            StateDataReporter(f"NPT_{temperature:.0f}K_stage{stage_idx}.log", LOG_FREQ,
                              step=True, temperature=True, potentialEnergy=True,
                              kineticEnergy=True, volume=True, density=True)
        )
        sim.reporters.append(DCDReporter(f"NPT_{temperature:.0f}K_stage{stage_idx}.dcd", TRJ_FREQ))

        sim.step(stage_steps)

        state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()
        print(f"  Stage {stage_idx} complete. PE: {state.getPotentialEnergy()}")

    save_checkpoint(sim, f"npt_{temperature:.0f}K.chk")
    print(f"\nNPT equilibration complete.")
    print(f"Box vectors: {box_vectors}")

    # ----------------------------------------------------------
    # 4. Production NPT — 3 replicas
    # ----------------------------------------------------------
    prod_time_ns = PROD_STEPS * TIMESTEP / 1000

    print("\n" + "="*70)
    print(f"  STAGE 4: PRODUCTION ({N_REPLICAS} replicas × {prod_time_ns:.0f} ns)")
    print("="*70)

    for rep in range(1, N_REPLICAS + 1):
        rep_dir = f"PROD_{temperature:.0f}K/rep{rep}"
        os.makedirs(rep_dir, exist_ok=True)
        seed = random.randint(1, 2**31 - 1)

        print(f"\n--- Replica {rep}/{N_REPLICAS} (seed={seed}) ---")
        t0 = time.time()

        system = build_system(prmtop)
        if box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*box_vectors)
        system.addForce(mm.MonteCarloBarostat(1 * unit.bar, temperature * unit.kelvin))
        integrator = set_integrator(temperature, timestep=TIMESTEP)

        sim = create_simulation(topology, system, integrator, platform_name, platform_props,
                                box_vectors=box_vectors)
        sim.context.setPositions(positions)
        sim.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)

        prefix = f"{rep_dir}/PROD_{temperature:.0f}K_rep{rep}"
        sim.reporters.append(
            StateDataReporter(sys.stdout, STDOUT_FREQ,
                              step=True, temperature=True,
                              potentialEnergy=True, volume=True,
                              density=True, speed=True)
        )
        sim.reporters.append(
            StateDataReporter(f"{prefix}.log", LOG_FREQ,
                              step=True, temperature=True, potentialEnergy=True,
                              kineticEnergy=True, volume=True, density=True)
        )
        sim.reporters.append(DCDReporter(f"{prefix}.dcd", TRJ_FREQ))

        chk = f"{prefix}.chk"
        load_checkpoint(sim, chk)
        sim.step(PROD_STEPS)
        save_checkpoint(sim, chk)

        elapsed = time.time() - t0
        print(f"  Replica {rep} complete. Wall time: {elapsed/3600:.1f} h")

    print(f"\nAll {N_REPLICAS} replicas complete at {temperature} K.")

    os.chdir(start_dir)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    print(f"OpenMM version: {mm.__version__}")
    print(f"Available platforms: "
          f"{[mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]}")

    PLATFORM_NAME, PLATFORM_PROPS = get_platform()

    wall_start = time.time()

    for temp in TEMPERATURES:
        run_simulation(temp, PRMTOP_FILE, INPCRD_FILE, PLATFORM_NAME, PLATFORM_PROPS)

    total_h = (time.time() - wall_start) / 3600
    print(f"\n{'='*70}")
    print(f"  ALL SIMULATIONS COMPLETE")
    print(f"  Total wall time: {total_h:.1f} h")
    print(f"{'='*70}")
