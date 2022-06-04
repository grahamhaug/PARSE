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
- Ensure no duplicate residues via InChi check
- Verifies that each processed SMILE has a * point of connectivity
- outputs the processed SMILEs to a .csv 
'''

sd_files = [file for file in os.listdir('.') if os.path.isfile(file) and file.endswith('.sdf')]

for sdfile in sd_files:

	#specify SDF containing mol data
	fileName = sdfile
	chopName = fileName.split('.',1)[0]

	#SMILEs flags to parse for
	#change according to the SMILEs strings you're looking to replace
	#Here, I am replacing carboxylic acid residues with * dummy atom
	chemFlags = ['C(O)=O', 'C(CC)=O', 'O=C(O)',
				'O=C(C)O', 'C(=O)O']

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
	removed_mw = []
	truncedSmiles = []
	for inpSmile in inputSmiles:
		#first calc. MW and remove anything over 1000 MU 
		calcMW = Descriptors.MolWt(Chem.MolFromSmiles(inpSmile))
		if calcMW < 1000:
			truncedSmiles.append(inpSmile)
		else:
			removed_mw.append(inpSmile)

	#how many SMILEs were removed due to MW criteria
	removed_mw_count = len(removed_mw)
	print(f'Removed {removed_mw_count} SMILEs due to MW exclusion criteria.')

	print("\nChecking for repeated residues...")
	#need to check if there are any repeat residues
	#check using InChIs
	captured_inchis = []
	for candidate in truncedSmiles:
		#convert smile to inchikey - can use moltoinchi but it has a lot of nonsensical output
		inchi = Chem.MolToInchiKey(Chem.MolFromSmiles(candidate,), options='-SNon')
		captured_inchis.append(inchi)

	print("\nDetermining SMILEs connectivity...")
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

	print(f'\nWriting {final_number} SMILEs to parsed_{fileName}.csv')
	#write the reduced output to csv
	parsedDF.to_csv(f'{chopName}.csv', index=False)



