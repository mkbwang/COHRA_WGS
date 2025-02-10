import pandas as pd
import numpy as np
import leidenalg as la
import igraph as ig
import networkx as nx
import matplotlib.pyplot as plt
import random

def graph_cluster_summary(source="saliva", type="taxa"):

    print(f"Processing {source} {type}")
    # load lasso results and SPIEC-EASI results
    lasso_selection = pd.read_csv(f'{source}_{type}_yr1.csv', index_col=0)
    adj = pd.read_csv(f"{source}_{type}_adjmat.tsv", sep="\t")

    adj_mat = adj.to_numpy()
    degrees = np.sum(adj_mat, axis=1)

    print(f"There are {np.sum(degrees == 0)} features out of {len(degrees)} that are not correlated with any other features")

    # check the features which are connected to at least one other feature
    subset_index = np.where(degrees > 0)[0]
    adj_subset = adj.iloc[subset_index, subset_index]
    feature_names = adj_subset.columns.values
    lasso_selection_subset = lasso_selection.loc[feature_names, :]
    adj_mat_subset = adj_subset.to_numpy()

    # leiden clustering
    graph = ig.Graph.Adjacency((adj_mat_subset > 0).tolist())
    partition = la.find_partition(graph, la.ModularityVertexPartition)
    cluster_labels = np.array(partition.membership)
    nG = nx.from_numpy_array(adj_mat_subset)

    # Map cluster labels to colors
    unique_labels = list(set(cluster_labels))
    print(f"The graph is clustered into {len(unique_labels)} clusters")
    colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_labels)))
    label_color_map = {label: color for label, color in zip(unique_labels, colors)}
    node_colors = [label_color_map[label] for label in cluster_labels]

    # plot the whole graph
    fig0, ax0 = plt.subplots(figsize=(4, 4))
    seed = 1
    random.seed(seed)
    pos = nx.spring_layout(nG)  # You can use other layouts as well
    nx.draw(nG, pos, ax=ax0,  node_color=node_colors, with_labels=False, node_size=50, font_color='white')

    fig0.savefig(f'plots_{type}/{source}/network.png', dpi=500, bbox_inches='tight')
    plt.close(fig0)


    # plot the graph in each cluster
    for j in range(np.max(cluster_labels) + 1):
        # select subset of the adjacency matrix as the subnetwork
        selected_index = np.where(cluster_labels == j)[0]
        selected_network = adj_mat_subset[np.ix_(selected_index, selected_index)]
        node_names = feature_names[selected_index]
        node_to_labels = dict(zip(np.arange(len(node_names)), node_names))
        selected_G = nx.from_numpy_array(selected_network)
        fig1, ax1 = plt.subplots(figsize=(4, 3))
        pos = nx.spring_layout(selected_G)
        nx.draw(selected_G, pos, ax=ax1, node_color="#1f96ad", with_labels=False, node_size=50, font_color='black')
        nx.draw_networkx_labels(selected_G, pos, labels=node_to_labels, font_color='black', font_size=5)
        fig1.savefig(f"plots_{type}/{source}/cluster_{j}.png", dpi=500, bbox_inches='tight')
        plt.close(fig1)

    # combine the lasso result and the network clustering result
    num_clusters = len(unique_labels)
    clusters_details = [[] for j in range(num_clusters)]
    for j in range(num_clusters):
        clusters_details[j] = feature_names[np.where(cluster_labels == j)[0]]

    for j in range(num_clusters):
        selected_features = clusters_details[j]
        coefs_mat = lasso_selection_subset.loc[selected_features, :].to_numpy()
        selection_times = np.sum(coefs_mat != 0, axis=1)
        directions = np.sum(coefs_mat, axis=1) > 0
        selected_index = np.where(cluster_labels == j)[0]
        selected_network = adj_mat_subset[np.ix_(selected_index, selected_index)]
        degrees = np.sum(selected_network, axis=1)
        selected_feature_information = pd.DataFrame({'Name': selected_features,
                                                     'Selected Times': selection_times,
                                                     'Degree': degrees,
                                                     "Positive Caries": directions})

        selected_feature_information.to_csv(f'plots_{type}/{source}/Cluster_{j}_information.csv', index=False)



if __name__ == '__main__':

    graph_cluster_summary(source="saliva", type="taxa")
    graph_cluster_summary(source="saliva", type="ko")
    graph_cluster_summary(source="plaque", type="taxa")
    graph_cluster_summary(source="plaque", type="ko")

