import pandas as pd
from Functions_for_Parsing import pan_for_residues, drop_duplicate_smiles	
from Functions_for_Parsing import remove_high_mw, remove_charged_smiles
from Functions_for_Parsing import remove_wonky_smiles, remove_wonky_FGs
from Functions_for_Parsing import subs_carboxylic_acids, split_by_substitution_num

"""
- Reads in raw SMILEs from Pubchem .csv file
- Retains sulfinate fragments only; removes Metals/Ligands/etc
- Outputs a formatted .csv that serves as input for building
- sulfoxides with iterative_sulfoxide_maker.py
""" 

### Chem Flags to Process ###
chem_flags = ['C(O)=O', 'C(CC)=O', 'O=C(O)', 'O=C(C)O', 'C(=O)O', 'C(=O)[O-]']

### Data Import ###
import_file_name = 'Pubchem_primary_acids.csv'
#import_file_name = 'small.csv'
chop_name = import_file_name.split('.',1)[0]
#using fields to not pull in the full mega frame
fields = ['isosmiles']
incoming_data = pd.read_csv(f'{import_file_name}', usecols=fields)

### Print a header ###
print(f'\nProcessing: {import_file_name}')

#convert 'isosmiles' column to list, drop null values
raw_smiles = incoming_data[incoming_data['isosmiles'].notnull()]
incoming_smiles = raw_smiles['isosmiles'].to_list()

#count the number of incoming smiles
incoming_length = len(incoming_smiles)
print(f'Incoming raw SMILEs: {incoming_length}')

### Pull in raw SMILEs, split multicomponents, drop metals,
#initial 
#Parse raw SMILEs from PubChem for sulfinate residues
structures = pan_for_residues(incoming_smiles, chem_flags)
pre = len(structures)


### Remove anything weird or over 1000 MW
#remove anything over 1000 mw
processed_structures = remove_high_mw(structures, 1000)

#remove any wonky functional Groups
processed_structures = remove_wonky_FGs(processed_structures)

#remove any wonky SMILEs
processed_structures = remove_wonky_smiles(processed_structures)


### Drop any duplicates by InchiKeys
#drop duplicate SMILEs using InchiKeys
candidate_structures = drop_duplicate_smiles(processed_structures)

#count the remaining structures
num_parsed = len(candidate_structures)
print(f'\nUnique structures before * substitution: {num_parsed}')

#replace chemical flags with * for combination
#also drops some duplicate strings
#Avoids replacing any esters
processed_structures = subs_carboxylic_acids(candidate_structures, chem_flags)

#split the list into mono/multi/error subs lists
mono_subs, multi_subs, error_subs = split_by_substitution_num(processed_structures)

#output the mono-substituted CAs as csv
mono_df = pd.DataFrame(columns=['SMILEs'])
mono_df['SMILEs'] = mono_subs
mono_df.to_csv(f'MonoSubs_{chop_name}.csv', index=False)

#output the multi-substituted CAs as csv
multi_df = pd.DataFrame(columns=['SMILEs'])
multi_df['SMILEs'] = multi_subs
multi_df.to_csv(f'MultiSubs_{chop_name}.csv', index=False)

#if there were any ester-only structures, output those, too. 
if len(error_subs) != 0:
	error_df = pd.DataFrame(columns=['SMILEs'])
	error_df['SMILEs'] = error_subs
	error_df.to_csv(f'EstersOnly_{chop_name}.csv', index=False)

#count the remaining structures
post_parsed = len(mono_subs) + len(multi_subs)
print(f'\nTotal unique structures: {post_parsed}')	
print(f'\nTotal Mono-substituted CAs: {len(mono_subs)}')

# #output the structures
# print(f'\n{len(mono_subs)} processed SMILES written to: Processed_{chop_name}.csv')
# output_df = pd.DataFrame(columns=['SMILEs'])
# output_df['SMILEs'] = mono_subs
# output_df.to_csv(f'Processed_{chop_name}.csv', index=False)	