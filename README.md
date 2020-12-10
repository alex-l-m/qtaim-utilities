# qtaim-utilities

Analysis based on Bader's quantum theory of atoms in molecules (QTAIM) can be done for free using a combination of GAMESS and AIMAll.
However, GAMESS input files can take some effort to prepare, and AIMAll outputs plaintext files designed for easy reading by humans, not for downstream statistical analysis.
In this repository I will be collecting utilities I've written to make my life easier during QTAIM analysis.

Currently includes a Python script for reading AIMAll mgp files and outputting critical point properties as comma separated values.
Also includes a Python script for parsing AIMAll sum files, but the output does not yet include all tables in the sum file.
Also includes a script for downloading named molecules from Chembl and creating an rdkit Mol object, as a first step towards automatically preparing GAMESS input from molecule names.
