# MD_SIMULATION

## Overview

### Split PDB complex 



### Protonate Ligand Using OpenBabel
`obabel -ipdb ligand.pdb -opdb -O ligand_H.pdb -p 7`


`antechamber -i lig.pdb -fi pdb -o lig.prep -fo prepi -j 4 -at amber -c bcc -nc 0`

`parmchk2 -i lig.prep -f prepi -o lig.frcmod`

