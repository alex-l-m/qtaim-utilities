import functools
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import rdDistGeom

def name2mol(name):
    smiles = new_client.molecule.search(name)[0]["molecule_structures"]["canonical_smiles"]
    partial_molecule = Chem.MolFromSmiles(smiles)
    full_molecule = Chem.AddHs(partial_molecule)
    Chem.rdDistGeom.EmbedMolecule(full_molecule)
    return full_molecule
