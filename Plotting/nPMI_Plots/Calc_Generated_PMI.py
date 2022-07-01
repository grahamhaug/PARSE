import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
import matplotlib as mpl
from matplotlib import pyplot as plt
from datetime import datetime

"""
- Calculates the first two NPRs for a series of data
- outputs a series of csv files that can be combined with stitcher
- output gets read in to plotting script
"""

### Data Import ###
### Generated Sulfoxide Data ###
gen_name = 'PubChem_generated_sulfoxides_Properties.csv'
raw_gen = pd.read_csv(f'{gen_name}')
trimmed_gen = raw_gen.where(raw_gen['Cm'] <= 1400).dropna().reset_index(drop=True)
gen_df = trimmed_gen.sample(frac=0.1, random_state=1, ignore_index=True)
gen_smiles = gen_df['SMILEs'].to_list()

### Convert to mol ### 
#Generated data
gen_molecules = []
for smile in gen_smiles:
    mol = Chem.MolFromSmiles(smile)
    gen_molecules.append(mol)
#check that length matches (no lost SMILEs on conversion)
print(len(gen_molecules))

#do the same for the generated sulfoxide set
output_name = "generated_sulfoxide_nprs_"

#proceed by chunks of 1000 mols for memory
output_chunk = 1000
chunk_count = 0
append_count = 1

#keep track of any error molecules
failed_embeds = []
count = 1
failCount = 0

#lists to store NPRs
gen_x1s = []
gen_y1s = []

for molecule in gen_molecules:
    #might remove H's after embedding; have to check
    molecule_wH = Chem.AddHs(molecule)
    #give all possible chances for embedding success
    test = AllChem.EmbedMolecule(molecule_wH, 
        enforceChirality=True, 
        useRandomCoords=True,
        ignoreSmoothingFailures=True)

    #RDkit returns '-1' for failed embeddings; catch them
    if test < 0:
        print(f"fail: {count}")
        #ensure that spacing/continuity is preserved
        gen_x1s.append('NaN')
        gen_y1s.append('NaN')
        failCount +=1
        continue

    #calculate NPRs
    x1 = Chem.Descriptors3D.NPR1(molecule_wH)
    y1 = Chem.Descriptors3D.NPR2(molecule_wH)
    #add them into a running list
    gen_x1s.append(x1)
    gen_y1s.append(y1)
    count += 1

    #this will output them into a separate csv..it's not ideal
    #but at least it wil be easy to paste into the full csv
    if len(gen_x1s) == output_chunk:
        csv_name = f"{output_name}{append_count}"
        outData = pd.DataFrame()
        outData['X1s'] = gen_x1s
        outData['Y1s'] = gen_y1s
        outData.to_csv(f'{csv_name}.csv', index=False, header=True)

        #clear lists after writing a chunk
        gen_x1s = [] 
        gen_y1s = [] 

        chunk_count = output_chunk + chunk_count 
        append_count += 1 
        #diagnostics
        now = datetime.now().strftime("%H:%M:%S:%f")
        print(f'{now}: {output_chunk} data written to {csv_name}.csv. Total: {chunk_count}')

#output remaining sulfoxides to a csv
outDatac = pd.DataFrame(columns=['X1s', 'Y1s'])
outDatac['X1s'] = gen_x1s
outDatac['Y1s'] = gen_y1s
outDatac.to_csv(f'{output_name}last.csv', index=False, header=True)
        
print('Gen PMIs done')
print(f'Failed sulfoxides: {failCount}')

