import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn_extra.cluster import KMedoids
from sklearn import metrics
from scipy.spatial.distance import euclidean
import itertools

#csv containing standardized mol prop data
#this dataset is z-score normalized/standardized
import_file_name = 'MinMaxScaled_sulfinates.csv'

#import prop data from csv
incoming_data = pd.read_csv(f'{import_file_name}')     

#consider only the numerical data: MW/FSP3/Bottcher
bottcher_data = incoming_data.iloc[:,3:6].values
s_scores = []
iterate = []
i=1
while i < 6:
    #set up kmeans model (runs 10 times and returns lowest SSE model)
    kmedoids = KMedoids(n_clusters = 6, init='k-medoids++', method='pam', max_iter=1000)

    #fit the model to the selected data
    fit_model = kmedoids.fit(bottcher_data)

    #keep clustering labels
    incoming_data['Cluster'] = kmedoids.labels_

    labels = kmedoids.labels_

    #retain the cluster centers (menoids)
    medoids = kmedoids.cluster_centers_
    #and the indices for output
    medoid_indices = kmedoids.medoid_indices_

    #print some clustering metrics
    reporting = metrics.silhouette_score(bottcher_data, labels, metric='euclidean')
    print(f'Silhouette Score: {reporting}')
    s_scores.append(reporting)

    #get the rows of incoming data corresponding to closest_pt_index numbers
    representative_sulfinates = pd.DataFrame()
    for index in medoid_indices:
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
    #print(kmedoids.cluster_centers_.shape)

    #plot medoids
    #kmedoids.cluster_centers_ has dimensions [6,3]
    #6 rows of 3 columns (6 instances of x,y,z data)
    #kmeans.cluster_centers_[:,0] => get all rows, first column (0 column => x elements)
    ax.scatter(
        kmedoids.cluster_centers_[:,0],
        kmedoids.cluster_centers_[:,1],
        kmedoids.cluster_centers_[:,2],
        s=250, c="y", edgecolor='k',
        label = "Medoids", marker='*',
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
    plt.savefig(f'kmedoids_k6_sulfinates_{i}b.png', bbox_inches="tight")

    iterate.append(i)
    i+=1

#this DF stores the results of each run
outData = pd.DataFrame(columns=['RunNum', 'Silhouette'])
outData['RunNum'] = iterate
outData['Silhouette'] = s_scores
outData.to_csv(f'KMedoids_Standardizded_Sampling_Performance.csv', index=False, header=True)