import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import euclidean
import itertools

#csv containing standardized mol prop data
#this dataset is z-score normalized/standardized
import_file_name = 'standardized_sulfinates.csv'

#import prop data from csv
incoming_data = pd.read_csv(f'{import_file_name}')     

#consider only the numerical data: MW/FSP3/Bottcher
bottcher_data = incoming_data.iloc[:,3:6].values
s_scores = []
iterate = []
i=1
while i < 21:
    #set up kmeans model (runs 10 times and returns lowest SSE model)
    kmeans = KMeans(n_clusters = 6, init="k-means++", n_init=10, max_iter=500)

    #fit the model to the selected data
    fit_model = kmeans.fit(bottcher_data)

    #keep clustering labels
    incoming_data['Cluster'] = kmeans.labels_

    labels = kmeans.labels_

    #retain the cluster centers
    centroids = kmeans.cluster_centers_
    # for centroid in centroids:
    #     print(centroid)
    # print("")
    # print(centroids)

    #print some clustering metrics
    reporting = metrics.silhouette_score(bottcher_data, labels, metric='euclidean')
    print(f'Silhouette Score: {reporting}')
    s_scores.append(reporting)

    # Loop over all clusters and find index of closest point to the cluster center and append to closest_pt_idx list.
    closest_pt_idx = []
    for iclust in range(kmeans.n_clusters):
        # get all points assigned to each cluster:
        cluster_pts = bottcher_data[kmeans.labels_ == iclust]
        # get all indices of points assigned to this cluster:
        cluster_pts_indices = np.where(kmeans.labels_ == iclust)[0]

        cluster_cen = kmeans.cluster_centers_[iclust]
        min_idx = np.argmin([euclidean(bottcher_data[idx], cluster_cen) for idx in cluster_pts_indices])
        
        # Testing:    
        # print('\nclosest point to cluster center: ', cluster_pts[min_idx])
        # print('closest index of point to cluster center: ', cluster_pts_indices[min_idx])
        # print('  ', bottcher_data[cluster_pts_indices[min_idx]])
        closest_pt_idx.append(cluster_pts_indices[min_idx])

    #get the rows of incoming data corresponding to closest_pt_index numbers
    representative_sulfinates = pd.DataFrame()
    for index in closest_pt_idx:
        captured = incoming_data.iloc[[index]]
        representative_sulfinates = pd.concat([captured, representative_sulfinates], axis=0)

    representative_sulfinates.to_csv(f'selected_sulfinates_{i}.csv', index=False, header=True)

    #add the cluster data to a new column in data
    # #separate the data by cluster
    data1 = incoming_data[incoming_data.Cluster==0]
    data2 = incoming_data[incoming_data.Cluster==1]
    data3 = incoming_data[incoming_data.Cluster==2]
    data4 = incoming_data[incoming_data.Cluster==3]
    data5 = incoming_data[incoming_data.Cluster==4]
    data6 = incoming_data[incoming_data.Cluster==5]

    # #initialize a 3D plot
    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(projection='3d')

    #call shape to determine structure of n-dimensional array: 
    #print(kmeans.cluster_centers_.shape)

    #plot centroids
    #kmeans.cluster_centers_ has dimensions [6,3]
    #6 rows of 3 columns (6 instances of x,y,z data)
    #kmeans.cluster_centers_[:,0] => get all rows, first column (0 column => x elements)
    ax.scatter(
        kmeans.cluster_centers_[:,0],
        kmeans.cluster_centers_[:,1],
        kmeans.cluster_centers_[:,2],
        s=250, c="y", edgecolor='k',
        label = "Centroids", marker='*',
        zorder=100)

    #data colouuurs (pinkies out):
    colors = itertools.cycle(['navy', 'magenta', 'red', 'cyan', 'limegreen', 'darkorange'])

    # add clustered data to plot
    ax.scatter(data1.MW, data1.FSP3, data1.Cm, c=next(colors), label='Cluster 1', linewidths=0.1, alpha=0.7)
    ax.scatter(data2.MW, data2.FSP3, data2.Cm, c=next(colors), label='Cluster 2', linewidths=0.1, alpha=0.7)
    ax.scatter(data3.MW, data3.FSP3, data3.Cm, c=next(colors), label='Cluster 3', linewidths=0.1, alpha=0.7)
    ax.scatter(data4.MW, data4.FSP3, data4.Cm, c=next(colors), label='Cluster 4', linewidths=0.1, alpha=0.7)
    ax.scatter(data5.MW, data5.FSP3, data5.Cm, c=next(colors), label='Cluster 5', linewidths=0.1, alpha=0.7)
    ax.scatter(data6.MW, data6.FSP3, data6.Cm, c=next(colors), label='Cluster 6', linewidths=0.1, alpha=0.7)

    # #tilt the plot
    ax.view_init(25, 125)

    # #set labels
    ax.set_xlabel('MW', fontsize=10, fontweight='bold')
    ax.set_ylabel('FSP3', fontsize=10, fontweight='bold')
    ax.set_zlabel('Bottcher Cm', fontsize=10, fontweight='bold')

    # #save the plot to a file
    plt.legend(loc='best')
    plt.savefig(f'kmeans_k6_sulfinates_{i}.png', bbox_inches="tight")

    iterate.append(i)
    i+=1

#this DF stores the results of each run
outData = pd.DataFrame(columns=['RunNum', 'Silhouette'])
outData['RunNum'] = iterate
outData['Silhouette'] = s_scores
outData.to_csv(f'KMeans_Sampling_Performance.csv', index=False, header=True)