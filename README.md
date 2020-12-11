# qtaim-utilities

Analysis based on Bader's quantum theory of atoms in molecules (QTAIM) can be done for free using a combination of GAMESS and AIMAll.
However, GAMESS input files can take some effort to prepare, and AIMAll outputs plaintext files designed for easy reading by humans, not for downstream statistical analysis.
In this repository I will be collecting utilities I've written to make my life easier during QTAIM analysis.

My workflow is as follows.

Create a GAMESS input file for a molecule, by downloading the structure from ChEMBL (which yields a graph representation), and generating atom coordinates with RDKit:

    python3 download_mol.py ammonia > ammonia.inp

Then, using my GAMESS installation, run GAMESS to get a wavefunction, stored in a ".dat" file:

    gms ammonia.inp

Extract the AIMPAC format wavefunction from the GAMESS output:

    python3 extract_wfn.py ammonia.dat > ammonia.wfn

Then, using my AIMAll installation, run AIMAll, generating a ".mgp" file containing critical point information, and a ".sum" file containing atomic properties:

    aimqb -nogui ammonia.wfn

Extract critical point properties into a comma-separated values (CSV) file, with columns "critical\_point", "property", and "value" (header included):

    python3 parse_mgp.py ammonia.mgp > ammonia_critpoints.csv

Extract atomic properties:

    python3 parse_sum.py ammonia.sum ammonia

This last step generates three output files, for properties of single atoms, pairs of atoms, and triplets of atoms, respectively.
The file names are constructed from a given prefix (the second argument), plus suffixes "\_oneatom.csv", "\_twoatom.csv", and "\_threeatom.csv".
The columns are a column containing the name of the table from the AIMAll ".sum" file, columns containing the atom names, a column containing the name of the property, and a column containing the value of the property.
No header is included.
Currently, not all tables from the ".sum" file are captured.
