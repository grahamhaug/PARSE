import pandas as pd
from Functions_for_Parsing import pan_for_residues, drop_canon_smiles	
from Functions_for_Parsing import remove_high_mw, remove_charged_smiles
from Functions_for_Parsing import remove_wonky_smiles#, subs_chem_flags

"""
- Reads in raw SMILEs from Pubchem .csv file
- Retains sulfinate fragments only; removes Metals/Ligands/etc
- Outputs a formatted .csv that serves as input for building
- sulfoxides with iterative_sulfoxide_maker.py
""" 

### Chem Flags to Process ###
sulf_flags = ['S(=O)']

### Data Import ###
import_file_name = 'PubChem_Known_Sulfoxides.csv'
chop_name = import_file_name.split('.',1)[0]
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

#parse raw SMILEs from PubChem for sulfoxides
sulfoxides = pan_for_residues(incoming_smiles, sulf_flags)
pre = len(sulfoxides)

#drop duplicate SMILEs using InchiKeys
sulfoxides = drop_canon_smiles(sulfoxides)

#remove anything over 1938 mw
processed_structures = remove_high_mw(sulfoxides, 1938)

#count the remaining sulfoxides
num_parsed = len(sulfoxides)
# print(f'\nUnique sulfoxides before * substitution: {num_parsed}')

#count the remaining sulfoxides
# post_parsed = len(processed_structures)
print(f'\nTotal unique sulfoxides: {num_parsed}')

#remove any non-neutral SMILEs
processed_structures = remove_charged_smiles(processed_structures)

#remove any wonky SMILEs
processed_structures = remove_wonky_smiles(processed_structures)

#output the processed_structures
print(f'{len(processed_structures)} SMILEs written to: Processed_{chop_name}.csv')
output_df = pd.DataFrame(columns=['SMILEs'])
output_df['SMILEs'] = processed_structures
output_df.to_csv(f'Processed_{chop_name}.csv', index=False)
