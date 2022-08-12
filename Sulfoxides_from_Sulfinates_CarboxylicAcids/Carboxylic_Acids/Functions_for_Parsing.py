import os
import re
import math
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import time


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
    Parses raw PubChem SMILEs, returns only fragments containing 'S'
    - Splits raw SMILEs into constituent fragments based on '.'
    - Parses each fragment for chem_flags
    - Returns a list of flag-containing residues
    Input: list of SMILEs featuring metals/fragments/ligands/etc
           list of Chemical Flags (SMILEs fragments) to parse for (ex: sulfinate characters)
    Output: list of SMILEs containing only chem-flagged residues
    """

    #print(f'Identifying fragments containing specified chemical flags...')
    panned_structures = []
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
                panned_structures.append(fragment)
                
    return panned_structures


def remove_wonky_smiles(smiles_list):
    """
    Removes SMILEs containing strange functional groups in DB
    Input: List of SMILEs
    Output: truncated list of SMILEs with no wonky substrings
    - Can also modify to remove whatever substring(s) of interest
    """
    print(f'\nChecking for wonky functional groups in SMILES...')
    
    pre_length = len(smiles_list)

    #I don't want SMILES containing any of these substrings
    #Isotopes, weird valences, metals, crazy things
    wonky_FGs = ['S(=O)(=O)O', 'C(O)(O)O', 'S(O)(O)O', 'OO', 'N(O)O', 'C(O)O',
                'N=C=O', 'N=C=S', 'N1Cl', 'NCl', 'N1Br', 'NBr', 'NI', 'N1I',
                'N=C=N', 'N=S', 'C(=S)', 'C(S)S', '[N+]#N', 'C(N)O', 'N=O',
                'C=C=O', 'C=C(=O)', 'S(=O)(=O)[O-]', 'Cl=', 'Br=', 'F=', 'I=',
                'OF', 'OCl', 'OBr', 'OI', '[Cl-]', '[Br-]', '[F-]', '[I-]',
                'P(=O)(O)O', 'N\\Cl', 'CPC', 'SO', 'SS',
                'N(F)', 'N(Cl)', 'N(Br)', 'N(I)',
                '[Si](F)', '[Si](Cl)', '[Si](Br)', '[Si](I)', 
                'SiH', 'OPN', 'C(I)', 'C(Br)', 'O[O' 'PO', '=S=O',
                'O=S=', '=S=S', 'S=S=', '=N=N', 'N=N=', '=[N]=[N]', '[N]=[N]=', '[NH3+]',
                '[2H]', '[3H]', '[11C]', '[11CH3]', '[12CH3]', '[13C]', '[13CH]',
                '[13CH2]', '[13CH3]', '[14C]', '[14CH]', '[14CH3]', '[15N]', '[15NH]', 
                '[17O]', '[18O]', '[18F]', '[19F]', '[C]', '[C+]', '[CH]',
                '[CH-]', '[CH+]', '[CH2]', '[CH2-]', '[O]', '[N]', '[NH+]',
                '[NH-]', '[S]', '[S+]', '[76Br]', '[B]', '[B-]', '[SiH]',
                '[Co]', '[P+]', '[Pt]', '[Fe]', '[Pd]', '[Po]', '[Ni]', '[Al]',
                '[As]', '[Cr]', '[Ga]', '[Ge]', '[Hg]', '[Pb]', '[Sb]', '[Se]',
                '[Si-]', '[Te]', '[Tl]', '[V]', '[W]', '=Cl', '=Br', '=I', '=F',
                'Cl(=O)', 'Br(=O)', 'I(=O)', 'F(=O)'
                ]
    
    new_list = []
    for smile in smiles_list: 
        if not any(wonky in smile for wonky in wonky_FGs):
            new_list.append(smile)

    post_length = len(new_list)
    removed = pre_length - post_length
    print(f'SMILES containing strange functional groups/valences removed: {removed}')

    return new_list             


def remove_charged_smiles(smiles_list):
    """
    Determines the net charge of each SMILE in list
    - Parses for # of +/- in each string, compares the value
    - if the sum = 0, the molecule is net 0 charged and is retained
    Input: a list of SMILEs
    Returns: A truncated list of net neutral structures
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

        #split based on net charge
        if pos_count == neg_count:
            neutral_smiles.append(smile)
        else:
            charged_smiles.append(smile)

    #give some diag on # charged
    num_charged = len(charged_smiles)
    print(f'{num_charged} non-neutral SMILEs removed.')     
    
    return neutral_smiles


def remove_high_mw(smiles_list, mw_cutoff):
    """
    Applies a MW cutoff to remove fragments with MW > 1500 MU
    Input: List of SMILEs, a variable MW value
    Returns: New list with all SMILES < cutoff MW
    """

    #timeit
    start=time.time()

    #set the cutoff MW, here:
    cutoff = mw_cutoff

    #Some lists for counting
    good_smiles = []
    removed_mw = []

    #remove anything over the cutoff MW value
    print(f"\nRemoving SMILES over {cutoff} MW...")

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
                print(f'\t\tFROWN: {smile}')
                continue

            calcMW = Descriptors.MolWt(mol)
            #print(calcMW)

            if calcMW < cutoff: #if good, save
                good_smiles.append(smile)

            else: # if too big, append elsewhere
                removed_mw.append(smile)

        batch_num +=1

    removed_mw_count = len(removed_mw)
    print(f'SMILES removed due to MW exclusion criteria: {removed_mw_count}')

    #timeit
    end = time.time()
    elapsed = end-start
    format_elapsed = "{:.2f}".format(elapsed)
    print(f'Sec Elapsed: {format_elapsed}')

    return good_smiles


def look_for_smarts(smiles_list):
    """
    Input: 
        A list of SMILES
    Output:
        A list of SMILES containing confirmed sulfides
    - Searches through the list of SMILES and retains only
      those structures that contain the SMARTS string of interst
    """

    #SMARTS definition
    #Hits only carboxylic acids; omits CAs bound to aromatic carbons
    carboxacid_smarts = Chem.MolFromSmarts('[CX3&!$(Cc)](=O)[OX2H1]')
    #sulfide_smarts = Chem.MolFromSmarts('[#16X2H0&!$([$([sr5]:[nr5,or5,sr5]),$([sr5]:[cr5]:[nr5,or5,sr5])])&!$(s1cccc1)&!$(SS)&!$(S-[CX3]=[OX1])]')

    #timeit
    start=time.time()

    frowns = []
    structures = [] #boom, ignore the spelling
    #work in batches of ~10,000 for memory
    chunked_list = process_in_batches(smiles_list)   
    num_chunks = len(chunked_list)
    print(f'\nSubstructure searching candidate SMILES in {num_chunks} batches...')
    batch_num = 1
    for chunk in chunked_list:
        print(f'\tBatch {batch_num}')
        for smile in chunk:
            #convert to mol
            mol = Chem.MolFromSmiles(smile)

            #Substructure search with SMARTS definition
            matches = mol.GetSubstructMatches(carboxacid_smarts)

            #count number of matches
            num_subs = len(matches)

            # #Store the SMILE if it contains at least one SMARTS match
            # if num_subs >= 1:
            #     structures.append(smile)

            #for carboxylic acids, want exactly 1 SMARTS match
            if num_subs ==1:
                structures.append(smile)

        batch_num += 1

    print(f'Structures containing SMARTS substring: {len(structures)}')
    
    #timeit
    end = time.time()
    elapsed = end-start
    format_elapsed = "{:.2f}".format(elapsed)
    print(f'Sec Elapsed: {format_elapsed}')

    return structures


def canonicalize_smiles(smiles_list,file_name):

    """
    Input:
        List of processed SMILES
    Output: 
        List of unique canonicalized SMILES (confirmed RDkit mols)
        A .csv containing the processed SMILES
    - Canonicalizes SMILES for consistency
    - Removes duplicate SMILES
    - Outputs a csv containing the truncated list of canonical SMILES
    """

    #timeit
    start=time.time()

    frowns = []
    canon_smiles = [] #boom, ignore the spelling
    #work in batches of ~10,000 for memory
    chunked_list = process_in_batches(smiles_list)   
    num_chunks = len(chunked_list)
    print(f'\nCanonicalizing SMILES in {num_chunks} batches...')
    batch_num = 1
    for chunk in chunked_list:
        print(f'\tBatch {batch_num}')
        for smile in chunk:
            canon_smile = Chem.CanonSmiles(smile)

            if canon_smile is None:
                frowns.append(canon_smile)
                #canons.append(np.nan)
            else:
                #add the canonical smile to list of canons
                canon_smiles.append(canon_smile)

        batch_num += 1

    #print a summary of SMILES processed and FROWNS if applicable
    print(f'SMILES canonicalized: {len(canon_smiles)}')
    if len(frowns) >= 1:
        print(f'{len(frowns)} Unreadable SMILES (FROWNs):')
        for frown in frowns:
            print(f'\t{frown}')

    ### Output CSV ###
    print(f'\nChecking for/removing duplicate SMILES...')
    #output a csv of the canonical smiles
    canons = pd.DataFrame(columns=['SMILES'])
    canons['SMILES'] = canon_smiles

    #check length of df before dropping duplicates
    pre_length = len(canons['SMILES'])
    
    #remove any duplicate residues based on CanonSmiles match
    parsed_df = canons.drop_duplicates(subset='SMILES')

    #parsed_df = parsed_df.reset_index(drop=True)
    post_length = len(parsed_df['SMILES'])
    duplicates = pre_length - post_length
    print(f'Duplicate SMILES removed: {duplicates}')

    if len(frowns) >= 1:
        print(f'{len(frowns)} Unreadable SMILES (FROWNs):')
        for frown in frowns:
            print(f'\t{frown}')

    #parsed_smiles = parsed_df[convert_df['SMILEs'].notnull()]
    smiles_list = parsed_df['SMILES'].to_list()

    #output a csv
    parsed_df.to_csv(f'{file_name}.csv', index=False)
    print(f'Processed SMILES sent to {file_name}.csv')

    #timeit
    end = time.time()
    elapsed = end-start
    format_elapsed = "{:.2f}".format(elapsed)
    print(f'Sec Elapsed: {format_elapsed}')

    return smiles_list