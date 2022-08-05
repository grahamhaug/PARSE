import os
import re
import math
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


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
		chunk = 10000
		size_n = math.ceil(size_of_list / chunk)

	chunked_list = np.array_split(big_list, size_n)
	return chunked_list


def pan_for_residues(smiles_list, chem_flags):

	""" 
	Parses raw PubChem SMILEs, returns only sulfinate residues
	- Splits raw SMILEs into constituent fragments based on '.'
	- Parses each fragment for hits on predefined chemical flags
	- Returns a list of flag-containing residues
	Input: list of SMILEs featuring metals/fragments/ligands/etc
		   list of Chemical Flags (SMILEs fragments) to parse for (ex: sulfinate characters)
	Output: list of SMILEs containing only sulfinate residues
	"""

	#print(f'Identifying fragments containing specified chemical flags...')
	sulfinates = []
	#loop over all the raw SMILEs
	for smile in smiles_list:
		#split SMILE into fragments
		fragment_list = smile.split('.')
		#remove duplicate fragments
		fragment_list = list(set(fragment_list))
		
		#from the list of frags, determine sulfinate-containing frag
		for fragment in fragment_list:
			#look for each sulfinate flag; count matches
			#counting prevents issues where multiple flags might exist in one frag
			search_hits = 0
			for flag in chem_flags:
				#have to escape flags for regex
				flag_escaped = re.escape(flag)

				#if a fragment contains any # of a given flag, iterate hit counts
				if bool(re.search(flag_escaped,fragment)):
					search_hits += 1

			#retain fragment if it has sulfinate flag(s)
			if search_hits >= 1:
				#print(fragment)
				sulfinates.append(fragment)
	return sulfinates


def drop_canon_smiles(smile_list):
    """
    - Test to see how RDKit's CanonSmiles() function works after reviewer suggestion
    - This function is ~2x faster than the previous 'drop_duplicate_smiles'
    - likely because no need to convert to RDKit mol object first
    """

    print(f'\nChecking for duplicate SMILEs using RDKit CanonSmiles...')

    pre = len(smile_list)
    
    frowns = []
    canons = [] #boom, ignore the spelling
    #work in batches of ~10,000 for memory
    chunked_list = process_in_batches(smile_list)   
    num_chunks = len(chunked_list)
    print(f'Processing CanonSmiles data in {num_chunks} batches...')
    batch_num = 1
    for chunk in chunked_list:
        print(f'\tBatch {batch_num}')
        for smile in chunk:
            #convert the mol to canonical SMILE
            canon = Chem.CanonSmiles(smile)

            #check that it's actually converted to something readable
            if canon is None:
                frowns.append(canon)
                canons.append(np.nan)
            else:
                #add the canonical smile to list of canons
                canons.append(canon)

        batch_num += 1

    #output a list of the canonical smiles
    canon_smiles = pd.DataFrame(columns=['SMILES', 'Canon_SMILES'])
    canon_smiles['SMILES'] = smile_list
    canon_smiles['Canon_SMILES'] = canons

    #check length of df first
    pre_length = len(canon_smiles['SMILES'])

    #remove any of the inchis that were non-convertible smiles (NaNs)
    parsed_df = canon_smiles.dropna().reset_index(drop=True)
    
    #remove any duplicate residues based on CanonSmiles match
    parsed_df = parsed_df.drop_duplicates(subset='Canon_SMILES')

    #parsed_df = parsed_df.reset_index(drop=True)
    post_length = len(parsed_df['SMILES'])
    duplicates = pre_length - post_length
    print(f'Duplicate Canon_SMILES removed: {duplicates}')

    if len(frowns) >= 1:
        print(f'{len(frowns)} Unreadable SMILES (FROWNs):')
        for frown in frowns:
            print(f'\t{frown}')

    #parsed_smiles = parsed_df[convert_df['SMILEs'].notnull()]
    smiles_list = parsed_df['SMILES'].to_list()

    #output a csv
    parsed_df.to_csv(f'chopped_smiles.csv', index=False)
    return smiles_list

def drop_duplicate_smiles(smile_list):
	"""
	Drops duplicate SMILEs from list based on Sets/InchiKeys
	- drops any exact duplicate SMILEs strings with set()
	- converts smiles to RDkit mol objects
	- prints out any SMILEs that couldn't be converted: "frowns"
	- converts good mol objects to InchiKeys
	- drops duplicate InchiKeys (more robust than SMILEs delet.)
	Input: list of sulfinates from pan_for_sulfinates
	Output: 2-column .csv file with SMILEs/InChI keys 
	Returns: List of SMILEs with most duplicates deleted
	"""

	print(f'\nChecking for duplicate SMILEs/InChIs...')

	#first, drop any duplicates from incoming list
	#first pass removal of duplicate SMILEs			
	pre = len(smile_list)
	smile_list = list(set(smile_list))
	post = len(smile_list)
	duplicates = pre - post
	print(f'Duplicate SMILEs removed: {duplicates}')

	#Now check for duplicate InChIs:
	#convert list of smiles to list of inchis
	inchis = []
	frowns = []
	#work in batches of ~10,000 for memory
	chunked_list = process_in_batches(smile_list)	
	num_chunks = len(chunked_list)
	print(f'Processing InChI data in {num_chunks} batches...')
	batch_num = 1
	for chunk in chunked_list:
		print(f'\tBatch {batch_num}')
		for smile in chunk:
			mol = Chem.MolFromSmiles(smile)

			#check for none types (unreadable smiles)
			if mol is None:
				frowns.append(smile)
				#print(f'Unfreadable Frown: {smile}')
				#ensures DF length is OK
				inchis.append(np.nan)
				continue

			#convert mol object to inchikey
			inchi = Chem.MolToInchiKey(mol)
			
			#add the inchi to a list of inchis
			inchis.append(inchi)

		batch_num += 1

	#construct a SMILE:InchiKey DF
	#retaining SMILEs without back-converting from inchis is likely safer	
	convert_df = pd.DataFrame(columns=['SMILEs', 'InChIs'])
	convert_df['SMILEs'] = smile_list
	convert_df['InChIs'] = inchis

	#check length of df first
	pre_length = len(convert_df['SMILEs'])

	#remove any of the inchis that were non-convertible smiles (NaNs)
	parsed_df = convert_df.dropna().reset_index(drop=True)
	
	#remove any duplicate residues based on InChI match
	parsed_df.drop_duplicates(subset='InChIs', inplace=True)

	parsed_df = parsed_df.reset_index(drop=True)
	post_length = len(parsed_df['SMILEs'])
	duplicates = pre_length - post_length
	print(f'Duplicate InChIs removed: {duplicates}')

	if len(frowns) >= 1:
		print(f'{len(frowns)} Unreadable SMILEs (frowns):')
		for frown in frowns:
			print(f'\t{frown}')

	#parsed_smiles = parsed_df[convert_df['SMILEs'].notnull()]
	smiles_list = parsed_df['SMILEs'].to_list()


	#output a csv
	parsed_df.to_csv(f'chopped_smiles.csv', index=False)
	return smiles_list


def remove_high_mw(smiles_list, mw_cutoff):
	"""
	Applies a MW cutoff to remove fragments with MW > 1500 MU
	Input: List of SMILEs, a variable MW value
	Returns: New list with all SMILES < cutoff MW
	"""

	#can set the cutoff MW, here:
	cutoff = mw_cutoff

	#Some lists for counting
	good_smiles = []
	removed_mw = []

	#remove anything over the cutoff MW value
	print(f"\nRemoving SMILES over {cutoff} MW...")

	#work in batches of ~10,000 for memory
	chunked_list = process_in_batches(smiles_list)	
	num_chunks = len(chunked_list)
	print(f'Processing InChI data in {num_chunks} batches...')
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

			calcMW = Descriptors.MolWt(mol)
			#print(calcMW)

			if calcMW < cutoff: #if good, save
				good_smiles.append(smile)

			else: # if too big, append elsewhere
				removed_mw.append(smile)

		batch_num +=1		

	removed_mw_count = len(removed_mw)
	print(f'SMILEs removed due to MW exclusion criteria: {removed_mw_count}')

	return good_smiles


def subs_chem_flags(smiles_list, chem_flags):
	"""
	Replaces chemical flags with wildcards for joining residues
	- Parses SMILEs for specific fragments
	- Replaces those fragments with "*" for later combination
	- Checks that no duplicate fragments exist after substitution
	- Returns a list of the substituted 
	- appends the SMILEs/InChI .csv with this new column (maybe)
	Input: a semi-parsed list of SMILEs (output from drop_duplicate_smiles)
	       list of chemical flags to replace
	Output: Returns a truncated list with all fragments wildcarded
	"""

	print(f"\nReplacing Substructure Flags with '*' Joints...")
	
	#count incoming SMILEs before substitution
	pre_count = len(smiles_list)
	print(f'Incoming SMILEs: {pre_count}')

	#replace flags with *'s
	substituted_smiles = []
	#loop looking for each flagged string
	for smile in smiles_list:
		for flag in chem_flags:
		#replaces any/all instances of flags
			smile = smile.replace(flag, "*")
		#append the substituted smile to a list
		substituted_smiles.append(smile)

	#print(f'Substituted SMILES: {len(substituted_smiles)}')

	#drop any duplicate SMILEs after subs.
	substituted_smiles = list(set(substituted_smiles))
	post_count = len(substituted_smiles)
	duplicates = pre_count - post_count
	print(f'Duplicate SMILEs after substitution: {duplicates}')
	print(f'Unique SMILEs: {post_count}')

	#return the truncated list
	return substituted_smiles


def remove_charged_smiles(smiles_list):
	"""
	Determines the net charge of each SMILE in list
	- Parses for # of +/- in each string, compares the value
	- if the sum = -1, the molecule is net -1 charged and is retained
	- Outputs the # of non -1 molecules removed
	Input: a list of SMILEs
	Returns: A truncated list of net -1 sulfinates
	"""
	print(f'\nRemoving SMILEs with net charges...')
	
	#Keep only net neutral fragments
	neutral_smiles = []
	#keep these, too, for diff
	charged_smiles = []

	for smile in smiles_list:
		#count charge types found
		pos_count = smile.count("-")
		neg_count = smile.count("+")

		#sulfinates should be net -1
		if pos_count == neg_count:
			neutral_smiles.append(smile)
		else:
			charged_smiles.append(smile)

	#give some diag on # charged
	num_charged = len(charged_smiles)
	print(f'{num_charged} non-neutral SMILEs removed.')		
	
	return neutral_smiles


def remove_wonky_smiles(smiles_list):
	"""
	Removes SMILEs containing strange fragments in DB
	Input: List of SMILEs
	Output: truncated list of SMILEs with no wonky substrings
	- Can also modify to remove whatever substring(s) of interest
	"""
	print(f'\nChecking for wonky substrings/isotopes in SMILEs...')
	
	pre_length = len(smiles_list)

	#I don't want SMILES containing any of these substrings
	#Isotopes, weird valences, metals
	wonky_substrings = ['[2H]', '[3H]', '[11C]', '[11CH3]', '[12CH3]', '[13C]', '[13CH]',
						'[13CH2]', '[13CH3]', '[14C]', '[14CH3]', '[15N]', '[15NH]', 
						'[17O]', '[18O]', '[18F]', '[19F]', '[C]', '[C+]', '[CH]',
						'[CH-]', '[CH+]', '[CH2]', '[CH2-]', '[O]', '[N]', '[NH+]',
						'[NH-]', '[S]', '[S+]', '[76Br]', '[B]', '[B-]', '[SiH]',
						'[Co]', '[P+]', '[Pt]', '[Fe]', '[Pd]', '[Po]', '[Ni]', '[Al]',
						'[As]', '[Cr]', '[Ga]', '[Ge]', '[Hg]', '[Pb]', '[Sb]', '[Se]',
						'[Si-]', '[Te]', '[Tl]', '[V]', '[W]'
						]

	new_list = []
	for smile in smiles_list: 
		if not any(wonky in smile for wonky in wonky_substrings):
			new_list.append(smile)

	post_length = len(new_list)
	removed = pre_length - post_length
	print(f'{removed} SMILEs containing isotopes/weird valences removed.')

	return new_list
				
def remove_wonky_FGs(smiles_list):
	"""
	Removes SMILEs containing strange functional groups in DB
	Input: List of SMILEs
	Output: truncated list of SMILEs with no wonky substrings
	- Can also modify to remove whatever substring(s) of interest
	"""
	print(f'\nChecking for wonky functional groups in SMILEs...')
	
	pre_length = len(smiles_list)

	#I don't want SMILES containing any of these substrings
	#Isotopes, weird valences, metals
	wonky_FGs = ['S(=O)(=O)O', 'C(O)(O)O', 'S(O)(O)O', 'OO', 'N(O)O', 'C(O)O',
				'N=C=O', 'N=C=S', 'N1Cl', 'NCl', 'N1Br', 'NBr', 'NI', 'N1I',
				'N=C=N', 'N=S', 'C(=S)', 'C(S)S', '[N+]#N', 'C(N)O', 'N=O',
				'C=C=O', 'C=C(=O)', 'S(=O)(=O)[O-]',
				'OF', 'OCl', 'OBr', 'OI',
				'P(=O)(O)O', 'N\\Cl', 'CPC', 'SO', 'SS',
				'N(F)', 'N(Cl)', 'N(Br)', 'N(I)',
				'[Si](F)', '[Si](Cl)', '[Si](Br)', '[Si](I)', 
				'SiH', 'OPN', 'C(I)', 'C(Br)', 'O[O' 'PO', ]
	
	new_list = []
	for smile in smiles_list: 
		if not any(wonky in smile for wonky in wonky_FGs):
			new_list.append(smile)

	post_length = len(new_list)
	removed = pre_length - post_length
	print(f'{removed} SMILEs containing strange functional groups removed.')

	return new_list				
