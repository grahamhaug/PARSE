import os
import re
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
#I know, I know (the star)
from helper_functions import *
from datetime import datetime

"""
- Builds a set of sulfoxides SMILEs from sets of processed CAs and Sulfinates
- Processes the data in chunks as the lists can become quite memory intensive
- Outputs a series of numbered .csv files (50k each) that can be combined into 
- one giant one with csv_stitcher.py
"""

#pull in SA csv (alphas)
SA_csv = pd.read_csv('selected_sulfinates_19.csv')
#get the smiles data and put it into a list
smiles_SA = SA_csv['SMILEs'].to_list()
num_of_alphas = len(smiles_SA)

#pull in a subset of CA data (beta)
CA_csv = pd.read_csv('scifinder_carboxylic_acids_calcprops.csv')
#get the smiles data and put it into a list
smiles_CA = CA_csv['SMILEs'].to_list()
num_of_betas = len(smiles_CA)

#determine how many are to be made and give start time
#starting time
start = datetime.now()
now = datetime.now().strftime("%H:%M:%S:%f")
purported_sulfoxides = num_of_alphas * num_of_betas
print(f'\n{now}: Expected total sulfoxides: {purported_sulfoxides}\n')

alphas = smiles_SA
betas = smiles_CA

output_alphas = []
output_betas = []
output_sulfoxides = []

#output cvs header:
output_name = "sulfoxide_smiles_"

#set the chunk interval
output_chunk = 500000
chunk_count = 0
append_count = 1

for alpha in alphas:
	for beta in betas:
		#determine #beta fragments needed and how
		#to number alpha fragment for substitution
		num_of_frags, alpha_joints = get_numbers_to_sub(alpha, beta)

		#alpha fragment is prepped for substitution, here:
		output_alpha = alpha_subs(alpha, alpha_joints)

		### BETA prep ###
		beta_frags = parse_beta(beta, alpha_joints)
		#output_betas.append(beta_frags)

		#print(sub_values)
		#output_smile = alpha_subs(alpha, alpha_joints)
		output_alphas.append(output_alpha)

		#build the sulfoxide output
		sulfoxide = f'{output_alpha}'
		for beta in beta_frags:
			sulfoxide += str(beta)

		#add built to a list up to 1000ish
		output_sulfoxides.append(sulfoxide)

		#once 500,000 sulfoxides in list, write to a numbered csv
		if len(output_sulfoxides) == output_chunk:
			#iteratively named csv_file
			csv_name = f"{output_name}{append_count}"
			#output sulfoxides to a csv
			outData = pd.DataFrame()
			outData['SMILEs'] = output_sulfoxides
			outData.to_csv(f'{csv_name}.csv', index=False, header=True)
			
			#clear list after writing a chunk
			output_sulfoxides = []	

			chunk_count = output_chunk + chunk_count 
			append_count += 1 
			#diagnostics
			now = datetime.now().strftime("%H:%M:%S:%f")
			print(f'{now}: {output_chunk} sulfoxides written to {csv_name}.csv. Total: {chunk_count}')

#output remaining sulfoxides to a csv
outDatac = pd.DataFrame(columns=['SMILEs'])
outDatac['SMILEs'] = output_sulfoxides
outDatac.to_csv(f'{output_name}last.csv', index=False, header=True)

#output details/diagnostics
output_count = chunk_count + len(output_sulfoxides)
end = datetime.now()
now = datetime.now().strftime("%H:%M:%S:%f")
time_elapsed = end - start
print(f"\n{now}: Completed parsing")
print(f"Total time: {time_elapsed}")
print(f"\n{output_count} sulfoxides written.")

