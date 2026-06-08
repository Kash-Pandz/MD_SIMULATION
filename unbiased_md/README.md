# Unbiased Molecular Dynamics Simulation

A modular OpenMM pipeline for preparing and simulating antibody structures starting from a crystal, homology or predicted structure.

# Steps

[ ] 1.  Input PDB obtained using crystal, homology or predicted structure. All structures need to be IMGT, Chothia etc numbered.
[ ] 2.  System Preparation
      [ ]  Disulphide detection
      [ ]  ACE/NME terminal caps 
      [ ]  Protonation (pdb2pqr + PROPKA)
      [ ]  Solvation + neutralisation
[ ] 3.  Two-stage minimisation (heavy atom -> unrestrained)
[ ] 4.  Heating ramp using NVT ensemble (100 K -> Target Temp)
[ ] 5.  NPT equilibration (5 stages by default and gradual restraint release)
[ ] 6.  Production NPT (n_replicas x prod_ns)
