import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import euclidean
from matplotlib import pyplot as plt

#for determining "optimal" number of clusters to fit
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)



### Import Data ###
#csv containing standardized mol prop data
import_file_name = 'standardized_sulfinates.csv'
#import prop data from csv
incoming_data = pd.read_csv(f'{import_file_name}')     
#store data: MW/FSP3/Bottcher
bottcher_data = incoming_data.iloc[:,3:6].values


### Configure Model ###
model = AgglomerativeClustering(n_clusters=6, affinity='euclidean', linkage='ward')
model = model.fit(bottcher_data)
clusters = model.labels_
reporting = metrics.silhouette_score(bottcher_data, clusters, metric='euclidean')
print(f'Silhouette Score: {reporting}')

### Dendrogram related ###
#dendrogram = sch.dendrogram(sch.linkage(bottcher_data, method='ward'))
#use a dendrogram to determine the optimal # of clusters
# plt.title("Hierarchical Clustering Dendrogram")
# plot_dendrogram(model, truncate_mode="level", p=5)
# plt.xlabel("# pts in each node")
# plt.show()

### Organize Output Data ###
#add the cluster data to a new column in data
incoming_data['Cluster'] = clusters
#separate the data by cluster
data1 = incoming_data[incoming_data.Cluster==0]
data2 = incoming_data[incoming_data.Cluster==1]
data3 = incoming_data[incoming_data.Cluster==2]
data4 = incoming_data[incoming_data.Cluster==3]
data5 = incoming_data[incoming_data.Cluster==4]
data6 = incoming_data[incoming_data.Cluster==5]

#calculate the centroid (arithmetic mean) of each cluster
data1_array = data1.iloc[:,3:6].values
data1_centroid = np.mean(data1_array, axis=0)

data2_array = data2.iloc[:,3:6].values
data2_centroid = np.mean(data2_array, axis=0)

data3_array = data3.iloc[:,3:6].values
data3_centroid = np.mean(data3_array, axis=0)

data4_array = data4.iloc[:,3:6].values
data4_centroid = np.mean(data4_array, axis=0)

data5_array = data5.iloc[:,3:6].values
data5_centroid = np.mean(data5_array, axis=0)

data6_array = data6.iloc[:,3:6].values
data6_centroid = np.mean(data6_array, axis=0)
print(data6_centroid)

#store cluster centroids in a list
centroids = [
    data1_centroid,
    data2_centroid,
    data3_centroid,
    data4_centroid,
    data5_centroid,
    data6_centroid,
    ]

#convert list to 2d np array (6x3)
centroids = np.array(centroids)

# Loop over all clusters and find index of closest point to the cluster center and append to closest_pt_idx list.
closest_pt_idx = []
for iclust in range(model.n_clusters):
    # get all points assigned to each cluster:
    cluster_pts = bottcher_data[model.labels_ == iclust]
    # get all indices of points assigned to this cluster:
    cluster_pts_indices = np.where(model.labels_ == iclust)[0]

    cluster_cen = centroids[iclust]
    min_idx = np.argmin([euclidean(bottcher_data[idx], cluster_cen) for idx in cluster_pts_indices])
    
    # Testing:    
    print(f'\ncluster center: {cluster_cen}')
    print('closest point to cluster center: ', cluster_pts[min_idx])
    print('closest index of point to cluster center: ', cluster_pts_indices[min_idx])
    print('  ', bottcher_data[cluster_pts_indices[min_idx]])
    closest_pt_idx.append(cluster_pts_indices[min_idx])

#get the rows of incoming data corresponding to closest_pt_index numbers
representative_sulfinates = pd.DataFrame()
for index in closest_pt_idx:
    captured = incoming_data.iloc[[index]]
    representative_sulfinates = pd.concat([captured, representative_sulfinates], axis=0)

representative_sulfinates.to_csv(f'selected_sulfinates.csv', index=False, header=True)

### Plotting Data ###
# #initialize a 3D plot
fig = plt.figure(figsize=(16,16))
ax = fig.add_subplot(projection='3d')

#plot centroids
#6 rows of 3 columns (6 instances of x,y,z data)
ax.scatter(
    centroids[:,0],
    centroids[:,1],
    centroids[:,2],
    s=250, c="y", edgecolor='k',
    label = "Centroids", marker='*',
    zorder=100)

# add clustered data to plot
ax.scatter(data1.MW, data1.FSP3, data1.Cm, c='blue', label='Cluster 1', linewidths=0.1, alpha=0.7)
ax.scatter(data2.MW, data2.FSP3, data2.Cm, c='magenta', label='Cluster 2', linewidths=0.1, alpha=0.7)
ax.scatter(data3.MW, data3.FSP3, data3.Cm, c='red', label='Cluster 3', linewidths=0.1, alpha=0.7)
ax.scatter(data4.MW, data4.FSP3, data4.Cm, c='cyan', label='Cluster 4', linewidths=0.1, alpha=0.7)
ax.scatter(data5.MW, data5.FSP3, data5.Cm, c='limegreen', label='Cluster 5', linewidths=0.1, alpha=0.7)
ax.scatter(data6.MW, data6.FSP3, data6.Cm, c='darkorange', label='Cluster 6', linewidths=0.1, alpha=0.7)

# #tilt the plot
ax.view_init(20, 125)

# #set labels
ax.set_xlabel('MW', fontsize=10, fontweight='bold')
ax.set_ylabel('FSP3', fontsize=10, fontweight='bold')
ax.set_zlabel('Cm (Bottcher)', fontsize=10, fontweight='bold')

# #save the plot to a file
plt.legend(loc='best')
plt.savefig('agglom6_bottch_euc_standardized.png', bbox_inches="tight")
