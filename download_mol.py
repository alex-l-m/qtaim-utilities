import functools
import sys
import re
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import rdDistGeom

def name2mol(name):
    smiles = new_client.molecule.search(name)[0]["molecule_structures"]["canonical_smiles"]
    partial_molecule = Chem.MolFromSmiles(smiles)
    full_molecule = Chem.AddHs(partial_molecule)
    Chem.rdDistGeom.EmbedMolecule(full_molecule)
    full_molecule.SetProp("_Name", name)
    return full_molecule

def parse_inp(indir):
    with open(indir, 'r') as f:
        text = f.read()
    outdict = dict()
    text = text.replace(" =", "=").replace("= ", "=")
    section = None
    for token in text.split():
        section_match = re.match(r"\$(.*)", token)
        if token == "$END":
            section = None
            continue
        if section_match is not None and section_match.group(1).upper() != "DATA":
            section = section_match.group(1).upper()
            outdict[section] = dict()
            continue
        variable_match = re.match(r"(.*)=(.*)", token)
        if variable_match is not None and section is not None:
            variable = variable_match.group(1).upper()
            value = variable_match.group(2)
            outdict[section][variable] = value
    return outdict

def dict2text(indict):
    return "\n".join(functools.reduce(list.__add__, [ \
            [" ${}".format(section.upper())] + \
                ["{}={}".format(par.upper(), str(val).upper()) for par, val in params.items()] + \
                [" $END", ""] \
            for section, params in indict.items()]))

def mol2data(mol):
    name = mol.GetProp("_Name")
    atom_coords = mol.GetConformers()[0].GetPositions().tolist()
    atom_symbols = [i.GetSymbol() for i in mol.GetAtoms()]
    atom_atomic_numbers = [i.GetAtomicNum() for i in mol.GetAtoms()]
    structure_csv = [("{}{}".format(symbol, i), str(atomic_number), str(coords[0]), str(coords[1]), str(coords[2])) \
            for i, symbol, atomic_number, coords in zip(range(len(atom_coords)), atom_symbols, atom_atomic_numbers, atom_coords)]
    return "\n".join([" $DATA", name, "C1"] + [", ".join(line) for line in structure_csv] + [" $END", ""])

inmol = name2mol(sys.argv[1])
ingamess = parse_inp(sys.argv[2])
outfile = dict2text(ingamess) + "\n" + mol2data(inmol)
sys.stdout.write(outfile)
