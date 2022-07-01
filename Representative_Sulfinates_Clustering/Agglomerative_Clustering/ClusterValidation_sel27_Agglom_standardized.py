from os import path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import pdist
import itertools as it
import matplotlib.cm as cm
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw

"""
- Perform 50 K-Means clustering runs and perfrom pairwise distance calcs.
- Reads in a 'properties' csv file containing SMILEs, MW, FSP3, Cm data 
- Standardizes data if it isn't already (optionally outputs .csv)
- Calculates Silhouette score and Min/Avg/Max pairwise distance for each run
- Outputs the results of each run to a single summary csv
- This script also will plot the data and gen a csv for each run with selected sulfs.
- This latter bit can be annoying over many runs, so comment out as needed
"""

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

prepped_data = incoming_data.iloc[:,2:5].values


### Prepare the Model ### 
#assign a number of clusters
nclusters = 27
model = AgglomerativeClustering(n_clusters=nclusters, affinity='euclidean', linkage='ward')

#fit the model to the selected data
model = model.fit(prepped_data)
clusters = model.labels_


### Silhouette Score ### 
s_scores = []
reporting = metrics.silhouette_score(prepped_data, clusters, metric='euclidean')
print(f'Silhouette Score: {reporting}')
s_scores.append(reporting)


### Separate data by clusters ###
#Add a column with each SMILEs' cluster #
incoming_data['Cluster'] = clusters
#group by the cluster # 
gb = incoming_data.groupby('Cluster')

centroids = []
for i in range(nclusters):
    #separate incoming data into bins by cluster number
    data = incoming_data[incoming_data.Cluster==i]
    data_array = data.iloc[:,2:5].values

    ### Also calculate Centroids ### 
    #calculate the arithmetic mean of the cluster (centroid)
    data_centroid = np.mean(data_array, axis=0)
    centroids.append(data_centroid)
    #print(data_centroid)

#convert to array
centroids = np.array(centroids)

#Loop over all clusters and find index of closest point to the cluster center and append to closest_pt_idx list.
closest_pt_idx = []
for iclust in range(model.n_clusters):
#for iclust in range(clusters):
    # get all points assigned to each cluster:
    cluster_pts = prepped_data[model.labels_ == iclust]
    # get all indices of points assigned to this cluster:
    cluster_pts_indices = np.where(model.labels_ == iclust)[0]

    cluster_cen = centroids[iclust]
    min_idx = np.argmin([euclidean(prepped_data[idx], cluster_cen) for idx in cluster_pts_indices])
    
    # Testing:    
    # print('\nclosest point to cluster center: ', cluster_pts[min_idx])
    # print('closest index of point to cluster center: ', cluster_pts_indices[min_idx])
    # print('  ', prepped_data[cluster_pts_indices[min_idx]])
    closest_pt_idx.append(cluster_pts_indices[min_idx])

#get the rows of incoming data corresponding to closest_pt_index numbers
representative_sulfinates = pd.DataFrame()
for index in closest_pt_idx:
    captured = incoming_data.iloc[[index]]
    representative_sulfinates = pd.concat([captured, representative_sulfinates], axis=0)

representative_sulfinates.to_csv(f'Selected_Sulfinates.csv', index=False, header=True)


### Output Image and SD File for Selected Sulfinates ### 
smiles_data = pd.read_csv(f'Selected_Sulfinates.csv', usecols=['SMILEs'])
selected_smiles = smiles_data['SMILEs'].to_list()

### Output an image of the selected SMILEs ### 
#convert to selected SMILEs into mol objects
mols = [Chem.MolFromSmiles(smile) for smile in selected_smiles]

#write the mols to a SD file
with Chem.SDWriter(f'Selected_Sulfinates.sdf') as w:
    for mol in mols:
        w.write(mol)
w.close()

num_list = []
for x in range(1,len(mols)+1):
    num_list.append(x)

image = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(550, 350), legends=[f'{x}' for x in num_list])
image.save(f'Selected_Sulfinates.png')


### Plotting the data ### 
# #initialize a 3D plot
fig = plt.figure(figsize=(16,16))
ax = fig.add_subplot(projection='3d')

#plot centroids
#kmeans.cluster_centers_ has dimensions [6,3]
#6 rows of 3 columns (6 instances of x,y,z data)
#kmeans.cluster_centers_[:,0] => get all rows, first column (0 column => x elements)
ax.scatter(
    centroids[:,0],
    centroids[:,1],
    centroids[:,2],
    s=250, c="y", edgecolor='k',
    label = "Centroids", marker='*',
    zorder=100)

#data colouuurs (pinkies out):
colors = cm.rainbow(np.linspace(0, 1, nclusters))

# add clustered data to plot
for i in range(nclusters):
    ax.scatter(gb.get_group(i).MW, gb.get_group(i).FSP3, gb.get_group(i).Cm,
        c=[colors[i]], label=f'Cluster {i+1}', linewidths=0.1, alpha=0.9, s=10)

# #tilt the plot
ax.view_init(10, 245)

# #set labels
ax.xaxis.set_rotate_label(False)
ax.yaxis.set_rotate_label(False)
ax.zaxis.set_rotate_label(False)
ax.set_xlabel(r'$\mathbf{MW}$', fontsize=20, fontweight='bold', labelpad=25)
ax.set_ylabel(r'$\mathbf{Fsp^3}$', fontsize=20, fontweight='bold', labelpad=20)
ax.set_zlabel(r'$\mathbf{C_m}$', fontsize=20, fontweight='bold', labelpad=25)

# #save the plot to a file
plt.legend(loc='best')
plt.savefig(f'Agglom_n{nclusters}_standardized_sulfinates.png', bbox_inches="tight")


### Calculate Pairwise Distance Data ### 
### Pairwise Distances ###
mins = []
maxs = []
avgs = []

file_name = f'Selected_Sulfinates.csv'

#within-run pairwise distances between selected sulfinates
num_data = pd.read_csv(f'{file_name}', usecols=['MW', 'FSP3', 'Cm'])

#convert numdata to numpy array
data_np = num_data.to_numpy()

#calculate the pairwise distance data
pairs = list(it.combinations(range(6),2))
d = pdist(data_np)
mins.append(d.min())
maxs.append(d.max())
avgs.append(d.mean()) 

#just 1 run in agglom
iterate = 1


### Output Summary CSV ###
#this DF summarizes results
outData = pd.DataFrame(columns=['RunNum', 'Silhouette', 'PD_Min', 'PD_Max', 'PD_Avg'])
outData['RunNum'] = iterate
outData['Silhouette'] = s_scores
outData['PD_Min'] = mins
outData['PD_Max'] = maxs
outData['PD_Avg'] = avgs
outData.to_csv(f'Agglom_Sampling_Performance-standardized.csv', index=False, header=True)