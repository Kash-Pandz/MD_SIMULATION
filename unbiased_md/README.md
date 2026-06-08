# Unbiased Molecular Dynamics Simulation

A modular OpenMM pipeline for preparing and simulating antibody structures starting from a crystal, homology or predicted structure.


## Steps

- [ ] 1. Input PDB obtained using crystal, homology or predicted structure. All structures need to be IMGT, Chothia etc numbered.
- [ ] 2. System Preparation (Disulphide detection, ACE/NME terminal caps, Protonation (pdb2pqr + PROPKA) and Solvation + neutralisation)
- [ ] 3. Two-stage minimisation (heavy atom -> unrestrained)
- [ ] 4. Heating ramp using NVT ensemble (100 K -> Target Temp)
- [ ] 5. NPT equilibration (5 stages by default and gradual restraint release)
- [ ] 6. Production NPT (n_replicas x prod_ns)


## Configuration

All numeric settings live in two dataclasses, `PrepConfig` and `MDConfig`. Defaults are tuned for an antibody Fv at 300K, 1 bar, 0.15M NaCl:

|     Parameter     |                Default                |              Notes              |
|:-----------------:|:-------------------------------------:|:-------------------------------:|
| pH                | 7.4                                   | pdb2pqr titration               |
| Force field       | amber14-all.xml                       | OpenMM-bundled                  |
| Water model       | TIP3P                                 | amber14/tip3p.xml               |
| Box padding       | 1.0 nm                                | distance to edge of box         |
| Ionic strength    | 0.15 M NaCl                           | with neutralisation             |
| Disulfide cutoff  | 2.5 Å                                 | SG-SG; bonded ~ 2.05 Å          |
| Timestep          | 4 fs                                  | with HMR + HBonds + rigid water |
| Hydrogen mass     | 4 amu                                 | HMR mass                        |
| NVT heating       | 100 K → T, 0.5 ns, 50 windows         | backbone restraints             |
| NPT equilibration | 5 × 1 ns, k = 5/2/1/0.1/0 kcal/mol/Å² | backbone restraints             |
| Production        | 100 ns × 3 replicas                   | unrestrained                    |


## Methodology

- [ ] Disulfide detection. Pairs of cysteine SG atoms within 2.5 Å are identified with MDAnalysis and the corresponding residues are renamed CYS → CYX so downstream tools treat them as bonded.
- [ ] Terminal capping. N- and C- termini capped with ACE (acetyl) and NME (N-methylamide) respectively. Cap atom positions are placed from standard Engh & Hugh peptide-bond parameters (C-N = 1.329 Å, C=O = 1.231 Å, C-CA = 1.522 Å, CA-C-N = 116.2°,
C-N-CA = 121.7°) using the Natural Extension Reference Frame (NeRF) algorithm (Parsons et al., 2005). Peptide bonds are placed
trans (ω = 180°); φ of the first residue is set to -120°, and ψ of
the last residue is derived from the existing carbonyl O so the new
amide lies in the correct peptide plane.
- [ ] Protonation. pdb2pqr3.x with PROPKA-driven titration-state prediction and the AMBER naming convention.
- [ ] Solvation. OpenMM `Modeller.addSolvent` with TIP3P water, 1.0 nm padding, and 0.15 M Na⁺/Cl⁻ (charge-neutralising)
