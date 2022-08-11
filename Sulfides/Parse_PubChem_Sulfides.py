import time
import pandas as pd
from Functions_for_Parsing import pan_for_residues, remove_wonky_smiles
from Functions_for_Parsing import remove_charged_smiles, remove_high_mw
from Functions_for_Parsing import look_for_smarts, canonicalize_smiles	

"""
- Reads in raw SMILEs from Pubchem .csv file
- Retains sulfide fragments only; removes Metals/Ligands/junk/etc
- Outputs a .csv file of canonical SMILES for confirmed sulfides
""" 

### Chem Flags to Process ###
#Make sure a candidate fragment actually contains an S
chem_flags = ['S']

### Data Import ###
import_file_name = 'PubChem_sulfides_alpha_sec.csv'
chop_name = import_file_name.split('.',1)[0]
#using fields to not pull in the full mega frame
fields = ['isosmiles']
incoming_data = pd.read_csv(f'{import_file_name}', usecols=fields)

#timeit
start=time.time()

### Print a header ###
print(f'\nProcessing: {import_file_name}')
 
#convert 'isosmiles' column to list, drop null values
raw_smiles = incoming_data[incoming_data['isosmiles'].notnull()]
incoming_smiles = raw_smiles['isosmiles'].to_list()

#count the number of incoming alleged SMILES from Pubchem
#not an accurate num since they contain fragments delimited by '.'
incoming_length = len(incoming_smiles)
print(f'Incoming raw SMILEs: {incoming_length}')

### Pull in raw SMILEs, split multicomponents delimited by '.'
#Parse raw SMILEs from PubChem for fragments containing 'S'
structures = pan_for_residues(incoming_smiles, chem_flags)
pre = len(structures)

#remove SMILEs containing wonky FGs
processed_structures = remove_wonky_smiles(structures)

#remove SMILES containing a net charge !=0
processed_structures = remove_charged_smiles(processed_structures)

#remove any molecules over 1000 mw
#also checks that SMILES can be converted to RDkit Mol objects
#(sometimes not the case)
processed_structures = remove_high_mw(processed_structures, 1000)

#Substructure search to confirm structures contain FG of interest
#need to edit the function with SMARTS definition of interest
smarts_structures = look_for_smarts(processed_structures)

#canonicalize the remaining SMILES
#drop duplicates; outputs a .csv in local dir
candidate_structures = canonicalize_smiles(smarts_structures)

#count the surviving structures
num_parsed = len(candidate_structures)
print(f'\nUnique structures after parsing: {num_parsed}')

#timeit
end = time.time()
elapsed = end-start
format_elapsed = "{:.2f}".format(elapsed)
print(f'Total sec elapsed: {format_elapsed}')