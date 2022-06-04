import pandas as pd
from sklearn.preprocessing import StandardScaler
'''
-Reads in a CSV file containing calculated properties: 
  [# of atoms, molecular weight, fraction C's that are SP3,
  Bottcher score (from Forli Labs), Bottcher Complexity per atom]
- Standardizes the data (z-score normalization) in each column
- Return a CSV file containing standardized data
'''

#import sulfinate csv containing calculated properties
#incoming file name
import_file_name = 'scifinder_sulfinates_calcprops.csv'
incoming_data = pd.read_csv(f'{import_file_name}')

#retain number and smiles columns
header_df = incoming_data.drop(columns=['NumAtoms', 'MW', 'FSP3', 'Cm'])

#transform only numerical data
num_df = incoming_data.drop(columns=['Number', 'SMILEs'])

#apply standardization to numerical column data
std_scaler = StandardScaler()
std_scaler

#create and output formatted csv
df_std = pd.DataFrame(std_scaler.fit_transform(num_df), columns=num_df.columns)
output_df = pd.concat([header_df, df_std], axis=1)
output_df.to_csv(f'standardized_sulfinates.csv', index=False, header=True)

