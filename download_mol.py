import functools
import sys
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

example_file = { \
        "contrl": { \
            "runtyp": "optimize",
            "scftyp": "rhf",
            "aimpac": ".TRUE."
        },
        "basis": {\
            "gbasis": "n31",
            "ngauss": 6,
            "ndfunc": 1,
            "npfunc": 1
        }, \
        "statpt": {\
            "method": "nr",
            "opttol": "0.00001",
            "nstep": 1000
        }}

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
outfile = dict2text(example_file) + "\n" + mol2data(inmol)
sys.stdout.write(outfile)
