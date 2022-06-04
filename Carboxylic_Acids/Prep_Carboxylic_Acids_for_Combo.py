import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

"""
- Combine the primary/secondary/tertiary sub-CSV's into one CSV
- Retain only single * residues (only 1 site of substitution)
- Calculate properties for all of the remaining residues
- Output parsed data into a finalized CA CSV ready for combination
"""

# #Import the primary, secondary, tertiary acids as DFs
# primary_csv = 'parsed_carboxylic_primary.csv'
# secondary_csv = 'parsed_carboxylic_secondary.csv'
# tertiary_csv = 'parsed_carboxylic_tertiary.csv'
# #convert CSV data to DFs
# primary_smiles = pd.read_csv(f'{primary_csv}') 
# secondary_smiles = pd.read_csv(f'{secondary_csv}') 
# tertiary_smiles = pd.read_csv(f'{tertiary_csv}')

# #combine primary, secondary, tertiary DFs into one working DF
# print("\nCombining input CSVs...")
# incoming_data = pd.concat(
# 	[primary_smiles, secondary_smiles, tertiary_smiles],
# 	ignore_index=True, axis=0)

# #get rid of the now-useless "Number" column
# incoming_data = incoming_data.drop(columns='Number')
# # #for print debug
# # incoming_smiles.to_csv(f'test.csv', index=False)

incoming_csv = 'combined_carboxylic_acids.csv'

incoming_data = pd.read_csv(f'{incoming_csv}') 

#pull in the SMILEs data from the combined DF
smiles = incoming_data['SMILEs'].to_list()

print("\nRemoving multi-replacement Carboxylic Acids...")
#check that '*' occurs only once in each SMILE
#only want to work with mono-sub acids
single_star = []
for smile in smiles:
	#count the number of '*'s in each smile
	star_count = 0 
	for char in smile:
		if char == '*':
			star_count = star_count +1

	#Save only molecules with a single B
	if star_count <= 1:
		single_star.append(smile)
# #count for debug
# print(len(single_star))

#count how many elements are in current SMILEs list
numData = []
molcount = len(single_star)
for value in range(1,molcount+1):
	numData.append(value)

### Mol. Properties Calculations ###
#lists to hold data
molWeights = []
molFPS3s = []
numAtoms = []
mols = []
for smile in single_star:
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
with Chem.SDWriter('CAmols.sdf') as w:
	for mol in mols:
		w.write(mol)
#calculate Bottcher Score
Bscore = BottchScore.runBS('CAmols.sdf')
os.remove("CAmols.sdf")

#write data to a CSV
outData = pd.DataFrame(columns=['Number', 'SMILEs', 'NumAtoms', 'MW', 'FSP3', 'Cm'])
outData['Number'] = numData
outData['SMILEs'] = single_star
outData['NumAtoms'] = numAtoms
outData['MW'] = molWeights
outData['FSP3'] = molFPS3s
outData['Cm'] = Bscore
outData.to_csv('scifinder_carboxylic_acids_calcprops.csv', index=False)





