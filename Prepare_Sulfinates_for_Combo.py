import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

"""
- Calculates molecular properties for SMILEs in incoming data
- Outputs a formatted .csv that serves as input for building
- sulfoxides with iterative_sulfoxide_maker.py
""" 

import_file_name = 'scifinder_sulfinates.csv'

#import a test case of smiles
incoming_data = pd.read_csv(f'{import_file_name}')

#get the smiles data and put it into a list
incoming_smiles = incoming_data['SMILEs'].to_list()

smiles_count = len(incoming_smiles)

#count how many SMILEs
numData = []
molcount = len(incoming_smiles)
for value in range(1,molcount+1):
	numData.append(value)

#Calculate properties for each SMILE
print(f"\nCalculating properties for {smiles_count} sulfinate SMILEs...")
#lists to store params
molWeights = []
molFPS3s = []
numAtoms = []
mols = []

for smile in incoming_smiles:
	#convert smile to mol object
	mol = Chem.MolFromSmiles(smile)
	#print(f"here is the first smile: {mol}")

	#calculate MW
	calcMW = Descriptors.MolWt(mol)
	#print(f'MW = {calcMW}')
	molWeights.append(calcMW)

	#calculate FSP3
	calcFSP3 = rdkit.Chem.Lipinski.FractionCSP3(mol)
	molFPS3s.append(calcFSP3)

	#count numAtoms (no H's)
	atomCount = mol.GetNumAtoms()
	numAtoms.append(atomCount)

	#add mol object to mols list
	mols.append(mol)

#write all the mols to sdf file for reading with BottchScore.py
#(slow but at least it works)
with Chem.SDWriter('testmols.sdf') as w:
	for mol in mols:
		w.write(mol)
#calculate Bottcher Score
Bscore = BottchScore.runBS('testmols.sdf')
os.remove("testmols.sdf")

#write data to a CSV
outData = pd.DataFrame(columns=['Number', 'SMILEs', 'NumAtoms', 'MW', 'FSP3', 'Cm'])
outData['Number'] = numData
outData['SMILEs'] = incoming_smiles
outData['NumAtoms'] = numAtoms
outData['MW'] = molWeights
outData['FSP3'] = molFPS3s
outData['Cm'] = Bscore
#outData['InChIs'] = alkylInchis

#uniqueDataInchis = outData.drop_duplicates(subset='InChIs')
outData.to_csv('scifinder_sulfinates_calcprops.csv', index=False)