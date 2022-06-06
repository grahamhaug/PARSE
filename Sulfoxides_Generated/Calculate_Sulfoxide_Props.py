import os
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

"""
- Use with Iterative_Sulfoxide_Generator.py output files
- Imports CSV files corresponding to constructed Sulfoxides
- Calculates #Atoms/MW/FSP3 Carbons/Bottcher Score/BottcherPerAtom Parameters
- Outputs a series of new CSV files called "{inputName}-CalcdProps.csv"
"""

#round up all the csv's in current directory
csv_files = [file for file in os.listdir('.') if os.path.isfile(file) and file.endswith('.csv')]

#report the findings
if csv_files:
	print(f"\nFound {len(csv_files)} .csv file(s) in local directory.")
else:
	print("No .csv file(s) found in local directory.")

#loop through each csv file rounded up
for csv_file in csv_files:

	#store the fileName for naming output file
	fileName = csv_file.rsplit('.',1)[0]

	#convert CSV to DF
	incoming_data = pd.read_csv(csv_file) 

	# #for debug
	# print("\nRaw input format (first several rows):")
	# sample = incoming_data.head()  
	# print(sample)

	#extract the SMILEs column
	smiles = incoming_data['SMILEs'].to_list()

	smiles_count = len(smiles)

	# #count how many SMILEs
	# numData = []
	# molcount = len(smiles)
	# for value in range(1,molcount+1):
	# 	numData.append(value)

	#Calculate properties for each SMILE
	print(f"\nCalculating properties for {smiles_count} sulfoxide SMILEs...")
	#lists to store params
	molWeights = []
	molFPS3s = []
	numAtoms = []
	mols = []

	frowns = []
	good_smiles = []
	for smile in smiles:
		#convert smile to mol object
		mol = Chem.MolFromSmiles(smile)
		#print(f"here is the first smile: {mol}")

		#gather the bad smiles
		if mol is None:
			frowns.append(smile)
			continue

		#debug record good smiles names in preferred format
		good_smiles.append(smile)

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

	num_frowns = len(frowns)
	print(f"\n{num_frowns} SMILEs could not be read.")

	#capture the bad smiles
	frownData = pd.DataFrame(columns=['SMILEs'])
	frownData['SMILEs'] = frowns
	frownData.to_csv(f'{fileName}_frowns.csv', index=False)


	# print("Converting processed mols into smiles for export...")
	# saved_smiles = []
	# for mol in mols:
	# 	converted_smile = Chem.MolToSmiles(mol)
	# 	saved_smiles.append(converted_smile)

		#count how many SMILEs
	numData = []
	molcount = len(good_smiles)
	for value in range(1,molcount+1):
		numData.append(value)

	#write all the mols to sdf file for reading with BottchScore.py
	#(slow but at least it works)
	counter = 0
	interval = 3000
	total = 0
	print("\nWriting mol data to temporary sd file...")
	with Chem.SDWriter('testmols.sdf') as w:
		for mol in mols:
			w.write(mol)
			if counter == interval:
				w.flush()
				counter = 0
				total = total + interval
				print(total)
			counter +=1
	w.close()
			
	print("\nCalculating Bottcher Score...")
	#calculate Bottcher Score
	Bscore = BottchScore.runBS('testmols.sdf')

	#write data to a CSV
	outData = pd.DataFrame(columns=['Number', 'SMILEs', 'NumAtoms', 'MW', 'FSP3', 'Cm'])
	outData['Number'] = numData
	outData['SMILEs'] = good_smiles
	outData['NumAtoms'] = numAtoms
	outData['MW'] = molWeights
	outData['FSP3'] = molFPS3s
	outData['Cm'] = Bscore
	outData.to_csv(f'{fileName}_CalcdProps.csv', index=False)



