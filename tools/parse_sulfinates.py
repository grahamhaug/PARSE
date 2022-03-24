import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

'''
Parses sdf file containing sulfinate salts
-removes metal from the salt
-replaces SOOH with SH
-Outputs a csv with smiles/inchi/MW/FSP3/BottchScore,etc.
'''

#specify SDF file containing sulfinate salts
fileName = 'sulfinate_salts_combined.sdf'

#parse the sulfinate salts sdf file
molecules = Chem.SDMolSupplier(fileName)

#make a list to capture raw smiles (with cations)
rawSmiles = []
#convert mols in sdf to smiles
for molecule in molecules:
		if molecule is None: continue
		#convert to smiles
		smi = Chem.MolToSmiles(molecule)
		rawSmiles.append(smi)

#make a list of smiles with metals chopped out
#these will be sulfinic acid smiles
sulfinicAcids = []
for smile in rawSmiles:
	#first calc. MW and remove anything over 1000 MU 
	calcMW = Descriptors.MolWt(Chem.MolFromSmiles(smile))
	if calcMW < 1000:
		#metals are at the end of SMILE in '.[X]' format
		#can remove using '.' delimiter
		choppedSmile = smile.rsplit('.',1)[0]
		sulfinicAcids.append(choppedSmile)

#trim off SOOH from SMILE and store the residue
#this list stores alkyl residues with "SH" to preseve conn.
alkylSmiles = []
#smiles expressions for sulfinic acid component
SA_Flags = ['S(=O)O', 'O=S(O)']
for acid in sulfinicAcids:
	for flag in SA_Flags:
		#replaces any/all instances of flags
		acid = acid.replace(flag, "B")
	alkylSmiles.append(acid)

#need to check if there are any repeat residues
alkylInchis = []
for candidate in alkylSmiles:
	#convert smile to inchikey - can use moltoinchi but it has a lot of nonsensical output
	inchi = Chem.MolToInchiKey(Chem.MolFromSmiles(candidate,), options='-SNon')#, logLevel=None, treatWarningAsError=False )
	alkylInchis.append(inchi)

#count how many elements are in alkylResidues (or the updated list with repeats removed)
numData = []
molcount = len(alkylSmiles)
for value in range(1,molcount+1):
	numData.append(value)

#lists to store data for csv
molWeights = []
molFPS3s = []
numAtoms = []
#holds mol objects
mols = []

#now should have a trimmed alkylSmiles list to work with
for molecule in alkylSmiles:
	#convert smile to mol object
	mol = Chem.MolFromSmiles(molecule)

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

#write all the mols to sdf file for reading with BottchScore.py
#(slow but at least it works)
with Chem.SDWriter('testmols.sdf') as w:
	for mol in mols:
		w.write(mol)
#calculate Bottcher Score
Bscore = BottchScore.runBS('testmols.sdf')
#calculate Bottcher per Atom
BPA = [float(bottch) / int(atoms) for bottch,atoms in zip(Bscore, numAtoms)]

#write data to a CSV
outData = pd.DataFrame(columns=['Num', 'SMILEs', 'InChIs', 'NumAtoms', 'MW', 'FSP3', 'Bottcher', 'BottchPA'])
outData['Num'] = numData
outData['SMILEs'] = alkylSmiles
outData['NumAtoms'] = numAtoms
outData['MW'] = molWeights
outData['FSP3'] = molFPS3s
outData['Bottcher'] = Bscore
outData['BottchPA'] = BPA
outData['InChIs'] = alkylInchis

#exclude nothing
#outData.to_csv('sulf-alykl.csv', index=False)

#exclude based on duplicate SMILEs
#uniqueDataSmiles = outData.drop_duplicates(subset='SMILEs')
#uniqueDataSmiles.to_csv('xsmile_sulf.csv', index=False)

#exclude based on duplicate InChIs
uniqueDataInchis = outData.drop_duplicates(subset='InChIs')
uniqueDataInchis.to_csv('xinchi_sulf.csv', index=False)
