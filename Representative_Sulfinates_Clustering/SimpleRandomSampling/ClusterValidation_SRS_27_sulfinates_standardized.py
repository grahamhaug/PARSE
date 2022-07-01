from os import path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from scipy.spatial.distance import pdist
import itertools as it

#import set of sulfinate residues containing unscaled calc. props
import_file_name = 'PubChem_Sulfinates.csv'


### Standardize Data ### 
#Check if a csv containing standardized data already exists in local dir
std_file = path.exists(f'standardized_{import_file_name}')

if std_file: #if yes, don't remake
    print(f"Using standardized data file in local directory")
    #pull in scaled numerical data
    incoming_data = pd.read_csv(f'standardized_{import_file_name}')
    #scaled_data = incoming_data.iloc[:,0:4].values
 
else: #if no, make one from import_file_name. 
    #capture the numerical data only
    unscaled = pd.read_csv(f'{import_file_name}', usecols=['MW', 'FSP3', 'Cm'])
    #save this for writing to new DF 
    other_data = pd.read_csv(f'{import_file_name}', usecols=['Number', 'SMILEs'])

    #apply std. scaler to standardize the data
    std_scaler = StandardScaler()
    #scale numerical data
    scaled_data = pd.DataFrame(std_scaler.fit_transform(unscaled), columns=unscaled.columns)

    #output a new spreadsheet with the three standardized columns
    output_df = pd.concat([other_data, scaled_data], axis=1)
    output_df.to_csv(f'standardized_{import_file_name}', index=False, header=True)
    incoming_data = output_df
    print(f'Standardized Data written to standardized_{import_file_name}')

#this DF will store the results of each run
outData = pd.DataFrame(columns=['RunNum', 'PD_Min', 'PD_Max', 'PD_Avg'])

pairs = list(it.combinations(range(27),2))
i=1
mins = []
maxs = []
avgs = []
iterate = []
#do 10k random samples
while i < 10001:
    #pull 6 without replacement
    random_selection = incoming_data.sample(n=27, replace=False, weights=None)
    
    # #can write selections to csv's to look at individual selections but it's annoying
    # #if you do this make sure to reduce the number of runs so you don't generate
    # #a billion files
    # random_selection.to_csv(f'random_selection_{i}.csv', index=False, header=True)

    #random_selection.to_csv(f'random_selection_{i}.csv', index=False, header=True)
    num_data = random_selection[['MW', 'FSP3', 'Cm']].copy()
    data_np = num_data.to_numpy()

    #calculate dist. matrix based on pairs and pairwise metrics
    d = pdist(data_np)
    mins.append(d.min())
    maxs.append(d.max())
    avgs.append(d.mean()) 
    iterate.append(i)

    # #print metrics per run; nice if smaller set
    # print("The smallest distance is {:}, and it occurs between sulfinate {:} and sulfinate {:}".format(d.min(), *pairs[d.argmin(axis=0)]))
    # print("The largest distance is {:}, and it occurs between sulfinate {:} and sulfinate {:}".format(d.max(), *pairs[d.argmax(axis=0)]))
    # print("The average distance between sulfinates is {:}".format(d.mean()))
    # print("")

    i+=1

#this DF stores the results of each run
outData = pd.DataFrame(columns=['RunNum', 'PD_Min', 'PD_Max', 'PD_Avg'])
outData['RunNum'] = iterate
outData['PD_Min'] = mins
outData['PD_Max'] = maxs
outData['PD_Avg'] = avgs
outData.to_csv(f'Random_Sampling_Performance_standardized.csv', index=False, header=True)
