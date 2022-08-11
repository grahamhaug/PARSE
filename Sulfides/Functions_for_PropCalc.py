import os
import re
import math
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
import BottchScore

def process_in_batches(big_list):
	"""
	Chops a list into sequential lists of ~10,000 entries 
	for incremental processing 
	Input: A big list
	Output: A list of smaller lists (~10,000)
	"""
	size_of_list = len(big_list)
	if size_of_list <= 10000:
		size_n = 1
	if size_of_list > 10000:
		chunk = 5000
		size_n = math.ceil(size_of_list / chunk)

	chunked_list = np.array_split(big_list, size_n)
	return chunked_list


def calc_mw(smiles_list):
	"""
	Calculate molecular weight for each SMILE in a list of SMILEs
	Input: A list of SMILEs
	Returns: A list of molecular weights
	"""

	print(f"\nCalculating Molecular Weight (MW)...")

	molWeights = []

	#work in batches of ~10,000 for memory
	chunked_list = process_in_batches(smiles_list)	
	num_chunks = len(chunked_list)
	print(f'Processing MW data in {num_chunks} batches...')
	batch_num = 1
	for chunk in chunked_list:
		print(f'\tBatch {batch_num}')
		for smile in chunk:

			#print(smile)
			mol = Chem.MolFromSmiles(smile)

			#check for none types (unreadable smiles)
			if mol is None:
				print(f'Frown: {smile}')
				continue

			#calculate the molecular weight using RDkit	
			calcMW = Descriptors.MolWt(mol)

			#append the MW to a list
			molWeights.append(calcMW)

		batch_num += 1	

	return molWeights


def calc_fsp3(smiles_list):
	"""
	Calculate fraction of SP3-hydridized carbons for each SMILE in a list of SMILEs
	Input: A list of SMILEs
	Returns: A list of FSP3s
	"""

	print(f"\nCalculating Fraction of SP3 Carbons (FSP3)...")

	molFSP3 = []

	#work in batches of ~10,000 for memory
	chunked_list = process_in_batches(smiles_list)	
	num_chunks = len(chunked_list)
	print(f'Processing MW data in {num_chunks} batches...')
	batch_num = 1
	for chunk in chunked_list:
		print(f'\tBatch {batch_num}')
		for smile in chunk:

			#print(smile)
			mol = Chem.MolFromSmiles(smile)

			#check for none types (unreadable smiles)
			if mol is None:
				print(f'Frown: {smile}')
				continue

			#calculate the FSP3 using Lipinski/RDkit
			calcFSP3 = rdkit.Chem.Lipinski.FractionCSP3(mol)

			#append the MW to a list
			molFSP3.append(calcFSP3)

		batch_num += 1	

	return molFSP3

def calc_Bottcher(smiles_list):
	"""
	Calculate Bottcher complexity score for each SMILE in a list of SMILEs
	Input: A list of SMILEs
	Returns: A list of Cm's 
	"""

	print(f"\nCalculating Bottcher Complexity (Cm)...")

	#this stores a list of lists
	#each run calculates Bscore for a fraction of the big list, returned as list
	#this is a list of those lists, it's unpacked/concatenated, below
	batches_of_Cm = []

	#work in batches of ~10,000 for memory
	chunked_list = process_in_batches(smiles_list)	
	num_chunks = len(chunked_list)
	print(f'Processing Cm data in {num_chunks} batches...')
	batch_num = 1
	print(f'\tBatch {batch_num}')

	#working in chunks of ~ 10k SMILEs
	for chunk in chunked_list:

		chunk_mols = []
		#for every smile of the <10k smiles in chunk, convert smile to a mol
		#add it to a list of mols
		for smile in chunk:

			#print(smile)
			mol = Chem.MolFromSmiles(smile)

			#check for none types (unreadable smiles)
			if mol is None:
				print(f'Frown: {smile}')
				continue

			chunk_mols.append(mol)

		print(f'\t\tNumber of mols in chunk: {len(chunk_mols)}')

		#Current BottchScore needs an sdf to calculate Cm
		#write this temp SDF in batches of 5k as it's quite resource intensive
		counter = 0
		interval = 3000
		total = 0
		print("\t\tWriting mol data to temporary sd file...")
		with Chem.SDWriter(f'temp_mols.sdf') as w:
			for mol in chunk_mols:
				w.write(mol)
				if counter == interval:
					w.flush()
					counter = 0
					total = total + interval
					#print(total)
				counter +=1
		w.close()

		#calculate the Cm for all the mols in the temp_mols.sdf
		chunk_Cm = BottchScore.runBS(f'temp_mols.sdf')
		print(f'Length of chunk_Cm: {len(chunk_Cm)}')

		#append this to the big list of Cms (for all the chunks)
		batches_of_Cm.append(chunk_Cm)

		#delete the temp file
		os.remove('temp_mols.sdf')

	#this unpacks the list of lists into a single list
	molCm = [ x for y in batches_of_Cm for x in y]	
	return molCm









