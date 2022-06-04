import os
import re
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from helper_functions import *


def biggest_num_in_smile(smile):
	'''
	- Returns [int]: highest_number: highest # in a smile as int
	- If no number, returns 0
	- Updated to handle weird %double_digit format of smiles..should be good to %99
	'''
	test_list = re.findall(r'%(\d{2})\d*', smile)

	leftover_string = re.sub('%\d{2}', '', smile)

	#To do: edit this to consider #'s > 9 of format %10, %19, %29, etc.
	if True in [char.isdigit() for char in leftover_string]:
		#determine the highest# in string
		num_smile = ""
		for m in leftover_string:
			if m.isdigit():
				num_smile = num_smile + m

		highest_single = max([int(c) for c in str(num_smile)])
		
	else:
		highest_single = 0
		
	test_list.append(highest_single)
	test_list = [int(i) for i in test_list]
	biggest_number = max(test_list)
	return biggest_number




def biggest_of_two_smiles(smile_a, smile_b):
	'''
	- RETURNS [int]: highest_a/b: highest # between two smiles as int
	'''
	highest_a = biggest_num_in_smile(smile_a)
	highest_b = biggest_num_in_smile(smile_b)

	if highest_a >= highest_b:
		return highest_a
	else:
		return highest_b




def get_numbers_to_sub(smile_a, smile_b):
	'''
	- alpha is a smile that will have beta residue tacked on
	- alpha can have multiple substitution sites => counts how many
	- substitution sites must be numbered, but those numbers must start
	  at an index higher than any number already present in alpha or beta
	- RETURNS [tuple]: sub_counter, mod_values
		- sub_counter [int]: count of how many subs. points are in alpha
	 	- mod_values [list]: a list of ints (substitution markers) that don't 
	 		conflict with any of the numbers already present in alpha/beta smiles
	'''
	#determine highest number in two smiles
	highest_number = biggest_of_two_smiles(smile_a, smile_b)
	#print(f"highest# in alpha and beta: {highest_number}")
	
	#determine number of subs needed in alpha smile
	sub_counter = smile_a.count('*')

	#if highest number n = 0, start numbering at 1
	if highest_number != 0: 
		mod_values = []
		for value in range(sub_counter):
			new_value = value + highest_number + 1
			if new_value > 9:
				new_value = f"%{new_value}"
				mod_values.append(new_value)
			else:
				mod_values.append(new_value)
			#mod_values.append(new_value)
	#if highest number n = x | x != 0, start numbering at n+1
	elif highest_number == 0:
		mod_values = []
		for value in range(sub_counter):
			new_value = value + 1
			if new_value > 9:
				new_value = f"%{new_value}"
				mod_values.append(new_value)
			else:
				mod_values.append(new_value)
			#mod_values.append(new_value)

	#return the substitution values 
	return sub_counter, mod_values 






def alpha_subs(smile, mod_values):
	'''
	INPUT [str, list]: alpha smile and list of # where
						fragments will be appended (starting
						at highest # between alpha/beta smiles)
	RETURNS [str]: smile: original smile with * replaced by
	structure of interest (here, 'S(=O)')
	'''
	#substitue in the sulfoxide structure & #s
	for x in mod_values:
	 	rstring = f'S{x}(=O)'
	 	smile = smile.replace("*", rstring, 1)
	return(smile)





def parse_beta(smile, alpha_joints):
	'''
	INPUT [str, list]: beta smile and a list of numbers where
						fragments will be connected to alpha
	RETURNS [list]: unique_frags: a list of length alpha_joints of
					prepared beta fragments w/ numbered connectivity
					corresponding to a given alpha
	'''
	#case1: SMILE starts with *
	if smile.startswith("*"):
		if smile[1].isalpha(): #is a letter
		#print("case 1")
			smile = smile[1] + f"Q" + smile[2:]
		else: 
			position = smile.find(']')
			if position == -1:
				print(f"catastrophic error with {smile}")
			smile = smile[1:position+1] + f"Q" + smile[position+1:]


	#case2: SMILE ends with *
	#this hopefully works - have a backup of the old one
	elif smile.endswith("*"):
		smile = f'{smile[0:-1]}Q'

	#case3: SMILE contains * but not at start/end
	elif "*" in smile:
		#print("case 2")
		if smile.find("(*)") > 0:
			#remove parentheses, replace with #
			smile = smile.replace("(*)", f"Q")
		else:
			smile = smile.replace("*", f"Q")	

	else: #no * in the SMILE
		return

	unique_frags = []
	for num in alpha_joints:
		unique_frag = smile.replace("Q", str(num))
	
		#return the formatted fragment for addition to end of alpha
		unique_frag = f".{unique_frag}"
		unique_frags.append(unique_frag)

	return unique_frags




# smile_a = 'CC(N(C(C)*)C(C)*)*'
# smile_b = '*C12c3c4c5c6c7c8c(c9c%10c1c1c3c3c%11c4c4c5c5c7c7c%12c8c8c9c9c%10c%10c1c1c3c3c%11c%11c4c4c5c7c5c7c%12c8c8c9c9c%10c1c1c3c3c%11c4c5c4c7c8c9c1c34)C62'

# num_of_frags, alpha_joints = get_numbers_to_sub(smile_a, smile_b)
# print(num_of_frags)
# print(alpha_joints)

# out_alpha = alpha_subs(smile_a, alpha_joints)
# print(out_alpha)

# beta_frags = parse_beta(smile_b, alpha_joints)
# print(beta_frags)

# sulfoxide = f'{out_alpha}'
# for beta in beta_frags:
# 	sulfoxide += str(beta)
# print()
# print(f"\n{sulfoxide}")