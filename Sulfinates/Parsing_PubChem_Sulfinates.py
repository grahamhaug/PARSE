import pandas as pd
from Functions_for_Parsing import pan_for_residues, drop_canon_smiles	
from Functions_for_Parsing import remove_high_mw, remove_charged_smiles
from Functions_for_Parsing import remove_wonky_smiles, remove_wonky_FGs 
from Functions_for_Parsing import subs_chem_flags

"""
- Reads in raw SMILEs from Pubchem .csv file
- Retains sulfinate fragments only; removes Metals/Ligands/etc
- Outputs a formatted .csv that serves as input for building
- sulfoxides with iterative_sulfoxide_maker.py
""" 

### Chem Flags to Process ###
sulf_flags = ['S(=O)O', 'O=S(O)', 'O=S([O-])', 'S(=O)[O-]', 'S(=O)[O]']

### Data Import ###
import_file_name = 'PubChem_Sulfinates_InputData.csv'
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

#parse raw SMILEs from PubChem for sulfinate residues
sulfinates = pan_for_residues(incoming_smiles, sulf_flags)
pre = len(sulfinates)

### Processing for weird/very heavy structures ###
#remove anything over 1500 mw
sulfinates = remove_high_mw(sulfinates, 1000)

#remove any wonky functional Groups
processed_structures = remove_wonky_FGs(sulfinates)

#remove any wonky SMILEs
processed_structures = remove_wonky_smiles(processed_structures)


### Do a final check before substitution ### 
#drop duplicate SMILEs using InchiKeys
processed_structures = drop_canon_smiles(processed_structures)


### Preparing to substitute sulfoxides ###
#count the remaining sulfinates
num_parsed = len(processed_structures)
print(f'\nUnique sulfinates before * substitution: {num_parsed}')

#replace chemical flags with * for combination
#also drops some duplicate strings
substituted_structures = subs_chem_flags(processed_structures, sulf_flags)

#remove any non-neutral SMILEs
substituted_structures = remove_charged_smiles(substituted_structures)

#count the remaining sulfinates
post_parsed = len(substituted_structures)
print(f'\nTotal unique sulfinates: {post_parsed}')


#output the sulfinates
print(f'\n{len(substituted_structures)} processed SMILES written to: Processed_{chop_name}.csv')
output_df = pd.DataFrame(columns=['SMILEs'])
output_df['SMILEs'] = substituted_structures
output_df.to_csv(f'Processed_{chop_name}.csv', index=False)
