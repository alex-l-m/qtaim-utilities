# qtaim-utilities

Analysis based on Bader's quantum theory of atoms in molecules (QTAIM) can be done for free using a combination of GAMESS and AIMAll.
However, GAMESS input files can take some effort to prepare, and AIMAll outputs plaintext files designed for easy reading by humans, not for downstream statistical analysis.
In this repository I will be collecting utilities I've written to make my life easier during QTAIM analysis.

Currently includes a Python script for reading AIMAll mgp files and outputting critical point properties as comma separated values.
I am working towards parsing the AIMAll sum files, and the repository includes functions for parsing the fixed-width tables they contain.
