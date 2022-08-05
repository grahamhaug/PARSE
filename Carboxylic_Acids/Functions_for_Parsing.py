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
	Output: list of SMILEs containing only CA residues
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
	convert_df = pd.DataFrame(columns=['SMILES', 'InChIs'])
	convert_df['SMILES'] = smile_list
	convert_df['InChIs'] = inchis

	#check length of df first
	pre_length = len(convert_df['SMILES'])

	#remove any of the inchis that were non-convertible smiles (NaNs)
	parsed_df = convert_df.dropna().reset_index(drop=True)
	
	#remove any duplicate residues based on InChI match
	parsed_df.drop_duplicates(subset='InChIs', inplace=True)

	parsed_df = parsed_df.reset_index(drop=True)
	post_length = len(parsed_df['SMILES'])
	duplicates = pre_length - post_length
	print(f'Duplicate InChIs removed: {duplicates}')

	if len(frowns) >= 1:
		print(f'{len(frowns)} Unreadable SMILES (frowns):')
		for frown in frowns:
			print(f'\t{frown}')

	#parsed_smiles = parsed_df[convert_df['SMILEs'].notnull()]
	smiles_list = parsed_df['SMILES'].to_list()


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


def get_bound_info_alpha(last_char, preceding_frag, smile):
    """ 
    Use when last_char (char before *) is type alpha
    - Checks if last char is C or not
        - need to add a condition if/if not C 
    - Returns Tuple "bound_info": (bound_atom, atom_index)
        - bound_atom is 'C', 'N', etc. 
        - atom_index is 0-based index of the bound_atom in original SMILE
    """

    #print(f"\talpha")

    # #double check incoming char is alpha
    # if last_char.isalpha():
    #     #print("Candidate immediately preceded by letter")
    # else: 
    #     print(f"Error: last_char {last_char} is not alpha")

    #check if bound atom is C
    if last_char == 'C':
        #print("Carbon bound")
        #print everything before the *
        #print(preceding_frag)
        
        # #add an arrow underneath to locate bound atom
        # space = ' '
        # carrot_pos = (len(preceding_frag) - 1)
        # print(f'{space*carrot_pos}^')

        #get the index of the bound atom
        bound_index = len(preceding_frag) - 1

    else:
        #print("NOT Carbon bound")
        #print(preceding_frag)
        
        # #add an arrow underneath to locate bound atom
        # space = ' '
        # carrot_pos = (len(preceding_frag) - 1)
        # print(f'{space*carrot_pos}^')

        #get the index of the bound atom
        bound_index = len(preceding_frag) - 1

        bound_info = (smile[bound_index], bound_index)

    #create a tuple with pertinent info
    bound_info = (smile[bound_index], bound_index) 

    #return bound atom and its index
    return bound_info


def get_bound_info_numeric(last_char, preceding_frag, smile):
    """ 
    Use when last_char (char before *) is type numeric
    - Checks if last char is C or not
        - need to add a condition if/if not C 
    - Returns Tuple "bound_info": (bound_atom, atom_index)
        - bound_atom is 'C', 'N', etc. 
        - atom_index is 0-based index of the bound_atom in original SMILE
    """

    #print(f"\tnumeric")

    #last character is a number; what is to the left of the number?

    #get the incoming fragment minus the preceding
    final_char_index = len(preceding_frag) - 1
    # print(final_char_index)

    one_left_frag = preceding_frag[:-1]
    #reasign the last_char one to the left of the number
    last_index = len(one_left_frag) - 1
    one_left_frag = preceding_frag[:last_index+1]

    last_char = one_left_frag[last_index]
    #print(f'char to left of #: {last_char}')

    #look 1 left of number if last char is number. 
    #there may be very rare cases where two numbers precede a number
    if last_char.isnumeric():
        last_index = last_index - 1
        #print("Numeric last char!!!!")
        #two_left_frag = one_left_frag[:-1]
        last_char = one_left_frag[last_index]


    # print(last_char)
    if last_char == '%':
        last_index = last_index - 1
        #print("Numeric last char!!!!")
        #two_left_frag = one_left_frag[:-1]
        last_char = one_left_frag[last_index]

    one_left_frag = preceding_frag[:last_index+1]

    #print(last_char)

    #check that the new last_char is alpha
    if last_char.isalpha():
        bound_info = get_bound_info_alpha(last_char, one_left_frag, smile)

    #number can also immediately follow a bracket
    elif last_char == ']':
        #print(last_char)
        #Count left until the bound character is a non-H alpha
        num_index = len(one_left_frag)-1
        while not(last_char.isalpha()) or last_char == 'H':
            num_index = num_index - 1

            #assign last_char as the first non-H alpha to the left of #
            last_char = smile[num_index]

        #chop a new frag up to the bound index
        num_index_frag = smile[:num_index+1]
        
        #now call func on alpha char
        bound_info = get_bound_info_alpha(last_char, num_index_frag, smile)
            
    else: #this is a bad error but I don't think it exists in practice...
        print("Error: Char to left of # is unanticipated...")
        print(smile)

    #return the bound atom and its index
    return bound_info


def get_atom_before_parenth(smile):
    """If a candidate follows a ')', have to count back (left) to see 
    what atom the candidate is connected to"""
    space = " "

    text = smile
    istart = []  # stack of indices of opening parentheses
    d = {}

    for i, c in enumerate(text):
        if c == '(':
             istart.append(i)
        if c == ')':
            try:
                d[i] = istart.pop()
            except IndexError:
                print('Too many closing parentheses')
    # if istart:  # check if stack is empty afterwards
    #     print('Too many opening parentheses')

    #print(d)

    #get the maximum-valued key
    #the max key has value = minimum index in SMILE
    #the min value is the earliest index prior to parentheses
    max_value = max(d)
    #print(max_value)

    #get the index of the character to the left of the opening (
    #(this is a key/value pair; the max key has the min value
    #the min value is the earliest index prior to parentheses
    #ex: C1=CC(=C(C(=C1)C(=O)O)C(=O)[O-]) returns index 4 for "C"
    bound_index = d.get(max_value) - 1
    bound_atom = smile[bound_index]

    #Make sure the ID'ed location is a letter
    #If the bound atom is not a a-z, look to the left of it
    #repeat until it is a-z
    while not(bound_atom.isalpha()):
        #print("looking to the left")
        #print(f'bound atom: {bound_atom}')
        
        #if the "atom" is a number, look to the char to its left
        if bound_atom.isnumeric():
            #print("hit number; one to the left")
            bound_index = bound_index - 1
            bound_atom = smile[bound_index]
            while not(bound_atom.isalpha()):
                #print(f'hit {bound_atom} - not alpha, look left')
                bound_index = bound_index - 1
                bound_atom = smile[bound_index]
            break

        #if bound_atom is a ), find index before preceding parentheses "capsule"
        if bound_atom == ')':
            bound_index = d[bound_index] - 1
            bound_atom = smile[bound_index]

        #if it's a bracket, need to be careful for H's and @s.
        if bound_atom == ']':
            #print("closing bracket hit")
            #continue looping until the bound character is a non-H letter
            while not(bound_atom.isalpha()) or bound_atom == 'H':
                bound_index = bound_index - 1
                bound_atom = smile[bound_index]

        else:
            #print(bound_atom)
            break
        #print(f'end loop bound atom: {bound_atom}')


    # if bound_atom == 'C':
    #     print(f"Carbon bound")
    #     print(smile)
    #     print(f"{space*bound_index}^")

    # else:
    #     #print("NOT Carbon bound")
    #     print(smile)
    #     print(f"{space*bound_index}^")

    #print(bound_atom)
    bound_info = (bound_atom, bound_index)
    #print(bound_info)
    return bound_info   


def get_bound_atom(smile,start,end):
    #pass start and end instead of chem flags
    #

    """determine if a given flag is bound to a carbon or not"""
    # preceded by a couple of options:
    # Atom: Is it C or not?
    # Number: What is the char to the left of # (should be a letter, I think)
    # (: Need to look to char to the left of this
    # ): Count back to the atom before the paired opening parenth (outside function)

    #Returns: a tuple with (bound_atom, atom_index)

    space = ' '

    #print("\nIncoming string:")
    #print(f'{smile}')

    preceding_frag = smile[:start]
    #print("String before located frag:")
    #print(preceding_frag)

    trailing_frag = smile[end:]
    #print("String after located frag:")
    #print(trailing_frag)

    #what is the last character of the fragment
    last_char = preceding_frag[-1]
    #print(f'Last char is: {last_char}')

    #last char is a letter: is it C or not?
    if last_char.isalpha():
        bound_info = get_bound_info_alpha(last_char, preceding_frag, smile)
        #print(bound_info)


    #last character is a number; what is to the left of the number?
    #it should always be an atom...I think
    elif last_char.isnumeric():
        #print("numeric")
        #print(preceding_frag)
        bound_info = get_bound_info_numeric(last_char, preceding_frag, smile)
        #print(bound_info)

    elif last_char == ')':
        #print("Candidate preceded by ')' parentheses")
        bound_info = get_atom_before_parenth(preceding_frag)
        #print(bound_info)

    elif last_char == '(':
        #print(f'\tNew branch pass')
        #look to the left of the opening parentheses
        one_left_frag = preceding_frag[:-1]
        last_char = one_left_frag[-1]
        #print(f"Char to left of '(': {last_char}")
        
        if last_char.isalpha():
            bound_info = get_bound_info_alpha(last_char, one_left_frag, smile)
            #print(bound_info)

        elif last_char.isnumeric():
            #print("foo")
            bound_info = get_bound_info_numeric(last_char, one_left_frag, smile)
            #print(bound_info)

        elif last_char == ')':
            #print(f'\tOut to get_atom_before_parenth: {preceding_frag}')
            bound_info = get_atom_before_parenth(preceding_frag)
            #print(bound_info)

        elif last_char == ']':
            #print("closing bracket hit")
            #get the index of the bracket
            bound_index = len(preceding_frag) - 1
            #print(bound_index)
            #continue looping until the bound character is a non-H letter
            while not(last_char.isalpha()) or last_char == 'H':
                bound_index = bound_index - 1
                last_char = smile[bound_index]

            #last_char should now be alpha; pass to alpha function
            #print(last_char)
            one_left_frag = preceding_frag[:bound_index+1]
            #print(one_left_frag)
            if last_char.isalpha():
                bound_info = get_bound_info_alpha(last_char, one_left_frag, smile)
                #print(bound_info)

    #number can also immediately follow a bracket
    elif last_char == ']':
        #print('new bracket')
        one_left_frag = preceding_frag[:-1]
        last_char = one_left_frag[-1]
        #print(last_char)
        #Count left until the bound character is a non-H alpha
        num_index = len(one_left_frag)-1
        while not(last_char.isalpha()) or last_char == 'H':
            num_index = num_index - 1

            #assign last_char as the first non-H alpha to the left of #
            last_char = smile[num_index]

        #chop a new frag up to the bound index
        num_index_frag = smile[:num_index+1]
        
        #now call func on alpha char
        bound_info = get_bound_info_alpha(last_char, num_index_frag, smile)

    #last char can be a slash if it's attached to an unsaturated C
    elif last_char == '\\' or last_char == '/':
        #print(f"\tSlash case")
        #need to count left to find the first processable char
        one_left_frag = preceding_frag[:-1]
        #print(f'one_left_frag: {one_left_frag}')
        last_char = one_left_frag[-1]
        #print(f'one char to left: {last_char}')

        # count to the left iteratively until a non-/ or \ char is hit
        n = -1
        i = 1 
        #print(f'\tloop{i}: {one_left_frag[n]}')
        while last_char == '\\' or last_char == '/' or last_char == '=':
            i = i + 1
            last_char = one_left_frag[n]
            n = n - 1 
            #print(f'\tloop{i}: {one_left_frag[n]}')

        #print(f'n = {n}')
        if n == -1:
            #print("foo")
            n_left_frag = preceding_frag[:n]
        elif n < -1:
            #print("faa")
            #have to count back n back chars, want n+1 to parse, though
            n_left_frag = preceding_frag[:n+1]

        #print(f'n_left_frag: {n_left_frag}')

        #determine the outcome
        if last_char.isalpha():
            #print("alpha)")
            bound_info = get_bound_info_alpha(last_char, n_left_frag, smile)
            #print(bound_info)

        elif last_char.isnumeric():
            #print("new case numeric")
            bound_info = get_bound_info_numeric(last_char, n_left_frag, smile)
            #print(bound_info)

        elif last_char == ')':
            #print("Candidate preceded by ')' parentheses")
            bound_info = get_atom_before_parenth(preceding_frag)
            #print(bound_info)

        elif last_char == '(':
            #print("Candidate in encapsulated branch")
            #look to the left of the opening parentheses
            one_left_frag = n_left_frag[:-1]
            last_char = one_left_frag[-1]
            #print(f'one char to left: {last_char}')
        
            if last_char.isalpha():
                bound_info = get_bound_info_alpha(last_char, one_left_frag, smile)
                #print(bound_info)

            elif last_char.isnumeric():
                bound_info = get_bound_info_numeric(last_char, one_left_frag, smile)
                #print(bound_info)

            elif last_char == ')':
                #print("Candidate preceded by ')' parentheses")
                bound_info = get_atom_before_parenth(preceding_frag)
                #print(bound_info)

    if 'bound_info' not in locals():
        print(f'Error: {smile}')

    return bound_info


def is_flag_bound_to_sp3_C(smile,bound_info):
    """
    - Checks if chem_flag is bound to a Carbon
    - Checks if the bound carbon is SP3 hybridized
    - Returns "TRUE" if both of these are true
    - Returns "FALSE" if one/both of these is/are false
    - INPUT:    smile = a complete SMILE
                bound_info = tuple containing (atom_type,atom_index)
                where:
                    bound_info contains info about atom that a Flag is bound to
                    atom_type is the ID of the bound atom
                    atom_index is the 0-indexed pos of bound atom in SMILE 
    """

    #Regex for determining saturation/unsaturation
    #re1 = r'=[^C]'
    re2 = r'.='
    re3 = r'[^\\]\\'
    re4 = r'#[^C]'
    re5 = r'.#'
    rightregex = r'\(='
    regexes = [re2, re3, re4, re5]
    combinedRegex = re.compile('|'.join('(?:{0})'.format(x) for x in regexes))

    #unpack incoming tuple
    (atom_type, atom_index) = bound_info

    #we want to replace all COOH bound to SP3 C's
    #first check if the bound atom is C 
    if atom_type == 'C':
        end_range = atom_index + 3
        start_range = atom_index - 2
        #region_around_C = smile[start_range:end_range]
        region_before_C = smile[start_range:atom_index]
        region_after_C = smile[atom_index+1:end_range]
 
        #print(f'L/R of C: {region_around_C}')
        #print(f'Left of C: {region_before_C}')

        #counts regex matches
        match = 0
        
        #check 2 char to the left of the bound C 
        if re.search(combinedRegex,region_before_C):
            #print(f"\tUnsaturated!")
            match += 1

        #print(f'Right of C: {region_after_C}')
        if re.search(rightregex,region_after_C):
            #print(f"\tUnsaturated!")
            match += 1

        if match >= 1:
            return False
        elif match == 0:   
            return True
        else:
            print("debacle detected")

    #if atom is not C, we don't want to replace it with *
    elif atom_type != 'C':
        return False



def subs_carboxylic_acids(smiles_list, chem_flags):
    """
    Replaces chemical flags with wildcards for joining residues
    - Parses SMILEs for specific fragments
    - Replaces those fragments with "*" for later combination
    - Checks that no duplicate fragments exist after substitution
    - Returns a list of the substituted SMILEs
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
    removed_frowns = []
    #loop looking for each flagged string
    for smile in smiles_list:
        #print(f'\nProcessing SMILE:')
        #print(f'{smile}')

        for flag in chem_flags:


            if flag in smile:

                #the incoming smile could just be the flag
                if flag == smile:
                    smile = smile.replace(flag,"*")

                #print(f'Flag found: {flag}')
                flag_length = len(flag)
                #print(flag_length)

                #count number of potential replacements
                max_replace = smile.count(flag)
                #print(f'Potential replacements: {max_replace}')
                replace_count = 0 

                #setup flags for regex
                flag_escaped = re.escape(flag)


                ### Possible Replacement Cases ### 
                #regexA: 'flag' is at the end of the line
                regexA = rf'{flag_escaped}$'

                start_pos = 0
                while len(re.findall(regexA,smile[start_pos:])) > 0:

                    #print(f"\nCase 1 - replace at the end")

                    m = re.search(regexA, smile[start_pos:])
                    # print("match1")
                    # print(m.group())
                    start = m.start() + start_pos
                    end = m.end() + start_pos

                    #determine what the flag is bound to 
                    bound_info = get_bound_atom(smile,start,end)

                    #check that the flag is bound to an SP3 C
                    test = is_flag_bound_to_sp3_C(smile,bound_info)

                    #Replace it with * if it is, don't if it's not
                    if test == True:
                        #print("CSP3 - Replace Flag")
                        smile = smile[:start] + "*" + smile[end:]
                        #print(smile)
                        replace_count += 1 
                    #if test == False:
                        #print("Unsaturated - Do not replace Flag")
                        #break

                    start_pos = end
                    remainder = smile[start_pos:]
                    #print(f'Leftover string: \n{remainder}')


                #regexB: 'flag' is between "(" and ")"
                regexB = rf'\({flag_escaped}\)'

                start_pos = 0
                while len(re.findall(regexB,smile[start_pos:])) > 0:
                    
                    #print(f"\nCase 2 - replace encapsulated")

                    #search returns the first (left-most) match
                    m = re.search(regexB, smile[start_pos:])

                    #calculate the indices of the match
                    start = m.start() + start_pos
                    end = m.end() + start_pos 

                    #determine what the flag is bound to 
                    bound_info = get_bound_atom(smile,start,end)

                    #check that the flag is bound to an SP3 C
                    test = is_flag_bound_to_sp3_C(smile,bound_info)

                    #Replace it with * if it is, don't if it's not
                    if test == True:
                        #print("CSP3 - Replace Flag")
                        smile = smile[:start] + "(*)" + smile[end:]
                        #print(smile)
                        start_pos = end - (flag_length - 2)
                        replace_count += 1
                    if test == False:
                        #print("Unsaturated - Do not replace Flag")
                        #break
                        start_pos = end

                    remainder = smile[start_pos:]
                    #print(f'Leftover string: \n{remainder}')
     

                #regexC: 'flag' is immediately followed by closing parentheses
                regexC = rf'{flag_escaped}\)'
                #return all the matches
                #match3 = re.findall(regexC,smile)
                #while matches are still being found, keep substituting
                #calculate new indices at each iteration
                start_pos = 0
                while len(re.findall(regexC,smile[start_pos:])) > 0:

                    #print(f"\nCase 3 - replace branched")

                    #search returns the first (left-most) match
                    m = re.search(regexC, smile[start_pos:])
                    #print("match3")
                    #print(m.group())
                    #print(smile)
                    #calculate the indices of the match
                    #print('Indices of located flag:')
                    start = m.start() + start_pos
                    #print(start)
                    end = m.end() + start_pos
                    #print(end)
                    #print("fragment check by indices:")
                    #print(smile[start:end])
                    #print("trailing:")
                    #print(smile[end:])


                    #print(f'start: {start}')
                    #print(f'end: {end}')
                    #determine what the flag is bound to 
                    bound_info = get_bound_atom(smile,start,end)

                    #check that the flag is bound to an SP3 C
                    test = is_flag_bound_to_sp3_C(smile,bound_info)

                    #Replace it with * if it is, don't if it's not
                    if test == True:
                        #print("CSP3 - Replace Flag")
                        smile = smile[:start] + "*)" + smile[end:]
                        #print("\nReplaced SMILE:")
                        #print(smile)
                        #print(f'end: {end}')
                        #need to offset by the length of the flag and parenth.
                        # start_pos = end - 5
                        start_pos = end - (flag_length - 1) 
                        #print(f'adjusted start_pos: {start_pos}')
                        replace_count += 1 
                    if test == False:
                        #print("Unsaturated - Do not replace Flag") 
                        start_pos = end
                        #break    

                    remainder = smile[start_pos:]
                    #print(f'Leftover string: \n{remainder}')
 
                #print("Processing complete.")
                #print("Final SMILE:")
                #print(smile)
                #print(f"Replacements made: {replace_count}")
            
            replacement = '*'
            if replacement in smile:
                #append the substituted smile to a list
                substituted_smiles.append(smile)
            # else:
            #     #append them to a different list for checking
            #     removed_frowns.append(smile)

    # #output the mono-substituted CAs as csv
    # frown_df = pd.DataFrame(columns=['SMILEs'])
    # frown_df['SMILEs'] = removed_frowns
    # frown_df.to_csv(f'removed_frowns.csv', index=False)            


    #drop any duplicate SMILEs after subs.
    dup_pre = len(substituted_smiles)
    substituted_smiles = list(set(substituted_smiles))
    post_count = len(substituted_smiles)
    duplicates = dup_pre - post_count
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
						'[13CH2]', '[13CH3]', '[14C]', '[14CH]', '[14CH3]', '[15N]', '[15NH]', 
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
				'SiH', 'OPN', 'C(I)', 'C(Br)', 'O[O' 'PO', '=S=O',
				'O=S=', '=S=S', 'S=S=', '=N=N', 'N=N=', '=[N]=[N]', '[N]=[N]=']
	
	new_list = []
	for smile in smiles_list: 
		if not any(wonky in smile for wonky in wonky_FGs):
			new_list.append(smile)

	post_length = len(new_list)
	removed = pre_length - post_length
	print(f'{removed} SMILEs containing strange functional groups removed.')

	return new_list				

	

def split_by_substitution_num(substituted_smiles):
	"""
	Splits the list of substituted SMILEs based on # of subs
	- One list with one substitution
	- One list with n > 1 substitutions
	Input: A list of substituted SMILEs
	Returns: A tuple of three lists, one mono-, one multi-subs, one error list
	"""

	main_list_length = len(substituted_smiles)

	mono_subs = []
	multi_subs = []
	#just in case
	bad_subs = []

	for smile in substituted_smiles:
		star_count =  0
		for char in smile:
			if char == '*':
				star_count = star_count + 1

		if star_count == 1:
			mono_subs.append(smile)
		if star_count >= 2:
			multi_subs.append(smile)
		if star_count == 0:	
			bad_subs.append(smile)

	print(f'\n{len(mono_subs)} mono-substituted SMILEs.')
	print(f'{len(multi_subs)} multi-substituted SMILEs.')
	print(f'{len(bad_subs)} non-substituted SMILEs (FROWNs).')

	return mono_subs, multi_subs, bad_subs



