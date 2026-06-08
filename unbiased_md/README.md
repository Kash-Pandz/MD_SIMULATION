# Unbiased Molecular Dynamics Simulation

A modular OpenMM pipeline for preparing and simulating antibody structures starting from a crystal, homology or predicted structure.


## Install
```bash
git clone https.

## Steps

- [ ] 1. Input PDB obtained using crystal, homology or predicted structure. All structures need to be IMGT, Chothia etc numbered.
- [ ] 2. System Preparation (Disulphide detection, ACE/NME terminal caps, Protonation (pdb2pqr + PROPKA) and Solvation + neutralisation)
- [ ] 3. Two-stage minimisation (heavy atom -> unrestrained)
- [ ] 4. Heating ramp using NVT ensemble (100 K -> Target Temp)
- [ ] 5. NPT equilibration (5 stages by default and gradual restraint release)
- [ ] 6. Production NPT (n_replicas x prod_ns)


## Configuration

# All numeric settings live in two dataclasses, `PrepConfig` and `MDConfig`. Defaults are tuned for an antibody Fv at 300K, 1 bar, 0.15M NaCl:


|     Parameter     |                Default                |              Notes              |
|:-----------------:|:-------------------------------------:|:-------------------------------:|
| pH                | 7.0                                   | pdb2pqr titration               |
| Force field       | amber14-all.xml                       | OpenMM-bundled                  |
| Water model       | TIP3P                                 | amber14/tip3p.xml               |
| Box padding       | 1.0 nm                                | edge of largest extent          |
| Ionic strength    | 0.15 M NaCl                           | with neutralisation             |
| Disulfide cutoff  | 2.5 Å                                 | SG–SG; bonded ~ 2.05 Å          |
| Timestep          | 4 fs                                  | with HMR + HBonds + rigid water |
| Hydrogen mass     | 4 amu                                 | HMR mass                        |
| NVT heating       | 100 K → T, 0.5 ns, 50 windows         | np.linspace                     |
| NPT equilibration | 5 × 1 ns, k = 5/2/1/0.1/0 kcal/mol/Å² | backbone restraints             |
| Production        | 100 ns × 3 replicas                   |                                 |

