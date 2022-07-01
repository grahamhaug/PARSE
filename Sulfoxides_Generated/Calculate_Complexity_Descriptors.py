import pandas as pd
# import rdkit
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Descriptors
# from rdkit.Chem import Lipinski
# import BottchScore
from Functions_for_PropCalc import calc_mw, calc_fsp3, calc_Bottcher

"""
Calculate chemical complexity descriptors for a list of SMILEs
- MW = molecular weight
- FSP3 = fraction of SP3-hybridized carbon atoms
- Cm = Bottcher complexity score
Input: A csv file containing a list of SMILEs
Output: A csv file containing 4 columns: SMILE, MW, FSP3, Cm
"""

### Data Import ###
import_file_name = 'PubChem_generated_sulfoxides.csv'
chop_name = import_file_name.split('.',1)[0]
#using fields to not pull in the full mega frame
fields = ['SMILEs']
incoming_data = pd.read_csv(f'{import_file_name}', usecols=fields)
incoming_smiles = incoming_data['SMILEs'].to_list()

#count incoming SMILEs
smiles_count = len(incoming_smiles)
print(f'\n{smiles_count} SMILEs read from: {import_file_name}')

#get a column to index the csv
numData = []
for value in range(1,smiles_count+1):
	numData.append(value)


### Calculate Properties ###
#calculate molecular weight
molWeights = calc_mw(incoming_smiles)

#calculate FSP3
molFSP3 = calc_fsp3(incoming_smiles)

#calculate Cm
molCm = calc_Bottcher(incoming_smiles)


### Data Output ###
#Output data to a csv 
outData = pd.DataFrame(columns=['Number', 'SMILEs', 'MW', 'FSP3', 'Cm'])
outData['Number'] = numData
outData['SMILEs'] = incoming_smiles
outData['MW'] = molWeights
outData['FSP3'] = molFSP3
outData['Cm'] = molCm
outData.to_csv(f'{chop_name}_Properties.csv', index=False)
