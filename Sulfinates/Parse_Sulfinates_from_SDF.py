import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

'''
- Reads incoming SD file and converts each mol to SMILE
- Excludes any SMILE with MW > 1000 MU
- Replaces COOH's with dummy atom(s) * to mark connectivity
- Verifies that each processed SMILE has a * point of connectivity
- Removes any duplicate entries based on matching inchis
- outputs the processed SMILEs to a .csv 
'''

#specify SDF containing mol data
fileName = 'sulfinates.sdf'
chopName = fileName.split('.',1)[0]

#SMILEs flags to parse for
#change according to the SMILEs strings you're looking to replace
#Here, I am replacing sulfinate residues with * dummy atom
chemFlags = ['S(=O)O', 'O=S(O)', 'O=S([O-])', 'S(=O)[O-]']

#parse the input sdf
molecules = Chem.SDMolSupplier(fileName)
#raw count of mols pulled in
raw_count = len(molecules)

#convert all of the molecules from sdf into smiles
print(f"\nConverting {raw_count} mols to SMILEs format...")
inputSmiles = []
#convert mols in sdf to smiles
for molecule in molecules:
	#ensure the mol is readable
	if molecule is None: continue
	#convert to SMILEs format
	smi = Chem.MolToSmiles(molecule)
	inputSmiles.append(smi)

print("\nRemoving any molecules over 1000 MU...")
#Remove anything > 1000 MU
#in this case, also remove metal cations
removed_mw = []
truncedSmiles = []
for inpSmile in inputSmiles:
	#first calc. MW and remove anything over 1000 MU 
	calcMW = Descriptors.MolWt(Chem.MolFromSmiles(inpSmile))
	if calcMW < 1000:
		#Remove any of the metals present
		#metals are at the end of SMILE in '.[X]' format
		#can remove using '.' delimiter
		inpSmile = inpSmile.rsplit('.',1)[0]
		truncedSmiles.append(inpSmile)
	else:
		removed_mw.append(inpSmile)

#how many SMILEs were removed due to MW criteria
removed_mw_count = len(removed_mw)
print(f'Removed {removed_mw_count} SMILEs due to MW exclusion criteria.')

#need to check if there are any repeat residues
print("\nChecking for repeated residues...")
#check using InChIs
captured_inchis = []
for candidate in truncedSmiles:
	#convert smile to inchikey - can use moltoinchi but it has a lot of nonsensical output
	inchi = Chem.MolToInchiKey(Chem.MolFromSmiles(candidate,), options='-SNon')
	captured_inchis.append(inchi)

print("\nDetermining SMILEs junctions...")
#loop through truncated smiles; replace all flags with wildcard
#holds parsed residues for export
captured_residues = []
#loop looking for each flagged string
for candidate in truncedSmiles:
	for flag in chemFlags:
		#replaces any/all instances of flags
		candidate = candidate.replace(flag, "*")
	captured_residues.append(candidate)	

#count how many SMILEs are left
numData = []
molcount = len(captured_residues)
for value in range(1,molcount+1):
	numData.append(value)

#output a pandas DF csv
outData = pd.DataFrame(columns=['Number', 'SMILEs', 'InChIs'])
outData['Number'] = numData
outData['SMILEs'] = captured_residues
outData['InChIs'] = captured_inchis
pre_length = len(outData['SMILEs'])

#first remove any duplicate residues based on InChI match
outData.drop_duplicates(subset='InChIs', inplace=True)
post_length = len(outData['SMILEs'])
duplicates = pre_length - post_length
print(f'Number of duplicate residues: {duplicates}')

#check that each smile has at least one * (integrity check)
parsedDF = outData[outData["SMILEs"].str.contains("\*")==True]

#this might break
final_number = len(parsedDF['SMILEs'].to_list())

num_removed = raw_count - final_number
print(f'\nTotal number of SMILEs lost to parsing: {num_removed}')

#change this to write the correct output
print(f'\nWriting {final_number} SMILEs to scifinder_sulfinates.csv')
#write the reduced output to csv
parsedDF.to_csv(f'scifinder_sulfinates.csv', index=False)



