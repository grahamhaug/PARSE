import os
import pandas as pd
import glob

# read each csv
# drop number column
# concatenate them vertically
# check for duplicate inchikeys
# recalculate the #s

# get data file names
path =r'C:\Computational Work\Gaussian\Sulfoxide 2022\datasci\PubChem_Sulfides\pieces'
filenames = glob.glob(path + "/*.csv")

dfs = []
for filename in filenames:
    dfs.append(pd.read_csv(filename))

# Concatenate all data into one DataFrame
big_frame = pd.concat(dfs, ignore_index=True)

#dropped_nums = big_frame.drop(columns=['Number'])

#remove duplicate inchikeys
#parsedDF = dropped_nums.drop_duplicates(subset='InChIs')

big_frame.to_csv('oxidized_sulfides.csv', index=False)



