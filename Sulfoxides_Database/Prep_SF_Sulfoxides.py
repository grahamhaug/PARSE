import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

"""
- Calculate properties for all of the remaining residues
- Output parsed data into a finalized Sulfoxide CSV for plotting
"""

incoming_data = pd.read_csv('parsed_scifinder_sulfoxides.csv')

#pull in the SMILEs data from the combined DF
smiles = incoming_data['SMILEs'].to_list()

#count how many elements are in current SMILEs list
numData = []
molcount = len(smiles)
for value in range(1,molcount+1):
	numData.append(value)

### Mol. Properties Calculations ###
#lists to hold data
molWeights = []
molFPS3s = []
numAtoms = []
mols = []
for smile in smiles:
	#convert smile to mol object
	mol = Chem.MolFromSmiles(smile)

	#calculate MW
	calcMW = Descriptors.MolWt(mol)
	molWeights.append(calcMW)

	#calculate FSP3
	calcFSP3 = rdkit.Chem.Lipinski.FractionCSP3(mol)
	molFPS3s.append(calcFSP3)

	#count numAtoms (no H's)
	atomCount = mol.GetNumAtoms()
	numAtoms.append(atomCount)

	#add mol object to mols list
	mols.append(mol)

### Bottcher Calculations ###
#write all the mols to sdf file for reading with BottchScore.py
#(slow but at least it works)
with Chem.SDWriter('Sulfmols.sdf') as w:
	for mol in mols:
		w.write(mol)
#calculate Bottcher Score
Bscore = BottchScore.runBS('Sulfmols.sdf')
os.remove("Sulfmols.sdf")

#write data to a CSV
outData = pd.DataFrame(columns=['Number', 'SMILEs', 'NumAtoms', 'MW', 'FSP3', 'Cm'])
outData['Number'] = numData
outData['SMILEs'] = smiles
outData['NumAtoms'] = numAtoms
outData['MW'] = molWeights
outData['FSP3'] = molFPS3s
outData['Cm'] = Bscore
outData.to_csv('properties_scifinder_sulfoxides.csv', index=False)





