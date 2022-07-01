from os import path
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from matplotlib import pyplot as plt
from sklearn_extra.cluster import KMedoids
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
- Reads in a 'properties' csv file contai ning SMILEs, MW, FSP3, Cm data 
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
std_file = path.exists(f'minmaxed_{import_file_name}')

if std_file: #if yes, don't remake
    print(f"Using minmaxed data file in local directory")
    #pull in scaled numerical data
    incoming_data = pd.read_csv(f'minmaxed_{import_file_name}')
    #scaled_data = incoming_data.iloc[:,0:4].values
 
else: #if no, make one from import_file_name. 
    #capture the numerical data only
    unscaled = pd.read_csv(f'{import_file_name}', usecols=['MW', 'FSP3', 'Cm'])
    #save this for writing to new DF 
    other_data = pd.read_csv(f'{import_file_name}', usecols=['Number', 'SMILEs'])

    #apply std. scaler to standardize the data
    min_max_scaler = MinMaxScaler()
    #scale numerical data
    scaled_data = pd.DataFrame(min_max_scaler.fit_transform(unscaled), columns=unscaled.columns)

    #output a new spreadsheet with the three standardized columns
    output_df = pd.concat([other_data, scaled_data], axis=1)
    output_df.to_csv(f'minmaxed_{import_file_name}', index=False, header=True)
    incoming_data = output_df
    print(f'Minmaxed Data written to minmaxed_{import_file_name}')

prepped_data = incoming_data.iloc[:,2:5].values

#perform i iterations of K-means clustering
#outputs i files containing the molecules closest to cluster centroids
s_scores = []
num_runs = 20
for k in range(1,num_runs+1):
    #assign a number of clusters
    clusters = 27
    kmedoids = KMedoids(n_clusters = clusters, init='k-medoids++', max_iter=500)

    #fit the model to the selected data
    fit_model = kmedoids.fit(prepped_data)

    #keep clustering labels
    incoming_data['Cluster'] = kmedoids.labels_

    labels = kmedoids.labels_

    #retain the cluster centers (menoids)
    medoids = kmedoids.cluster_centers_
    #and the indices for output
    medoid_indices = kmedoids.medoid_indices_

    #get the rows of incoming data corresponding to closest_pt_index numbers
    representative_sulfinates = pd.DataFrame()
    for index in medoid_indices:
        captured = incoming_data.iloc[[index]]
        representative_sulfinates = pd.concat([captured, representative_sulfinates], axis=0)

    representative_sulfinates.to_csv(f'selected_sulfinates_{k}.csv', index=False, header=True)

    ### Output Image and SD File for Selected Sulfinates ### 
    smiles_data = pd.read_csv(f'Selected_Sulfinates_{k}.csv', usecols=['SMILEs'])
    selected_smiles = smiles_data['SMILEs'].to_list()

    ### Output an image of the selected SMILEs ### 
    #convert to selected SMILEs into mol objects
    mols = [Chem.MolFromSmiles(smile) for smile in selected_smiles]

    #write the mols to a SD file
    with Chem.SDWriter(f'Selected_Sulfinates_{k}.sdf') as w:
        for mol in mols:
            w.write(mol)
    w.close()

    num_list = []
    for x in range(1,len(mols)+1):
        num_list.append(x)

    image = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(550, 350), legends=[f'{x}' for x in num_list])
    image.save(f'Selected_Sulfinates_{k}.png')


    #print some clustering metrics
    reporting = metrics.silhouette_score(prepped_data, labels, metric='euclidean')
    print(f'Silhouette Score: {reporting}')
    s_scores.append(reporting)


    #add the cluster data to a new column in data
    cluster_data = {}
    for i in range(clusters):
        cluster_data[i] = incoming_data[incoming_data.Cluster==i]

    # #initialize a 3D plot
    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(projection='3d')


    #plot centroids
    #kmedoids.cluster_centers_ has dimensions [6,3]
    #6 rows of 3 columns (6 instances of x,y,z data)
    #kmedoids.cluster_centers_[:,0] => get all rows, first column (0 column => x elements)
    ax.scatter(
        kmedoids.cluster_centers_[:,0],
        kmedoids.cluster_centers_[:,1],
        kmedoids.cluster_centers_[:,2],
        s=250, c="y", edgecolor='k',
        label = "Medoids", marker='*',
        zorder=100)

    #data colouuurs (pinkies out):
    colors = cm.rainbow(np.linspace(0, 1, clusters))

    # add clustered data to plot
    for i in range(clusters):
        ax.scatter(cluster_data[i].MW, cluster_data[i].FSP3, cluster_data[i].Cm, 
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
    plt.savefig(f'kmedoids_k{clusters}_minmaxed_sulfinates_{k}.png', bbox_inches="tight")

    #k+=1

#j = 1
#evaluate pairwise distances between selected sulfinates
mins = []
maxs = []
avgs = []
iterate = []
for j in range(1,num_runs+1):
    file_name = f'selected_sulfinates_{j}.csv'
    #print(txt)

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
    iterate.append(j)

    # print(f'Dist.Min: {d.min()}')
    # print(f'Dist.Max: {d.max()}')
    # print(f'Dist.Mean: {d.mean()}')
    # print("The smallest distance is {:}, and it occurs between sulfinate {:} and sulfinate {:}".format(d.min(), *pairs[d.argmin(axis=0)]))
    # print("The largest distance is {:}, and it occurs between sulfinate {:} and sulfinate {:}".format(d.max(), *pairs[d.argmax(axis=0)]))
    # print("The average distance between sulfinates is {:}".format(d.mean()))
    # print("")

    #j+=1

#this DF stores the results of each run
outData = pd.DataFrame(columns=['RunNum', 'Silhouette', 'PD_Min', 'PD_Max', 'PD_Avg'])
outData['RunNum'] = iterate
outData['Silhouette'] = s_scores
outData['PD_Min'] = mins
outData['PD_Max'] = maxs
outData['PD_Avg'] = avgs
outData.to_csv(f'kmedoids_Sampling_Performance-minmaxed.csv', index=False, header=True)