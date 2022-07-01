import pandas as pd

"""
Extract a representative sample of carboxylic acids
- Applies SRS with a 0.25% sampling rate to each of 
	primary, secondary, tertiary acids
- Adds the sampled structures to a new DF/csv
"""

#name incoming CSVs
primary_acid_csv = 'MonoSubs_Pubchem_primary_acids.csv' 
secondary_acid_csv = 'MonoSubs_Pubchem_secondary_acids.csv' 
tertiary_acid_csv = 'MonoSubs_Pubchem_tertiary_acids.csv' 

#read in the data
primary_raw = pd.read_csv(f'{primary_acid_csv}')
secondary_raw = pd.read_csv(f'{secondary_acid_csv}')
tertiary_raw = pd.read_csv(f'{tertiary_acid_csv}')

# #convert to DFs
# primary_acids = primary_raw.where(primary_raw['Cm'] <= 1400).dropna().reset_index(drop=True)
# secondary_acids = secondary_raw.where(secondary_raw['Cm'] <= 1400).dropna().reset_index(drop=True)
# tertiary_acids = tertiary_raw.where(tertiary_raw['Cm'] <= 1400).dropna().reset_index(drop=True)

#sample DFs at 0.25% using SRS 
srs_primary = primary_raw.sample(frac=0.0025, random_state=1)
srs_secondary = secondary_raw.sample(frac=0.0025, random_state=1)
srs_tertiary = tertiary_raw.sample(frac=0.0025, random_state=1)

#dataframe roundup
dataframes_to_append = [srs_primary, srs_secondary, srs_tertiary]

#make a new DF to store the sampled data
#sampled_acids = pd.DataFrame(columns=['Number', 'SMILEs', 'MW', 'FSP3', 'Cm'])
sampled_acids = pd.concat(dataframes_to_append)

#output a CSV with the sampled data
sampled_acids.to_csv(f'sampled_acids.csv', index=False)





