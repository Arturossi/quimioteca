#!/usr/bin/python3

from rdkit import Chem
from rdkit.Chem import PandasTools

my_sdf = './Conformer3D_CID_2519.sdf'

frame = PandasTools.LoadSDF(my_sdf, smilesName='SMILES', molColName='Molecule', includeFingerprints=False)

frame.to_dict()

# usar esse
suppl = Chem.SDMolSupplier(my_sdf)
for mol in suppl:
	print(mol.GetNumAtoms())


