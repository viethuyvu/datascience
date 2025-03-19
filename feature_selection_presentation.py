import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
import plotly.graph_objects as go

def assess_individual_cluster(cluster, nomralized_jaccard_index_matrix,cluster_center):
    if len(cluster) == 1:
        return 0
    total_squared_distance_from_center = 0
    for npi in cluster:
        if npi != cluster_center:
            distance= 1-nomralized_jaccard_index_matrix.loc[npi, cluster_center]
            total_squared_distance_from_center += distance**2
    
    return total_squared_distance_from_center/(len(cluster)-1)

def assign_cluster_modify(nomralized_jaccard_index_matrix, nps, cluster1, cluster2, cluster3):
    for npi in nps:
        if npi != cluster1[0] and npi != cluster2[0] and npi != cluster3[0]:
            cluster1_score = nomralized_jaccard_index_matrix.loc[npi, cluster1[0]]
            cluster2_score = nomralized_jaccard_index_matrix.loc[npi, cluster2[0]]
            cluster3_score = nomralized_jaccard_index_matrix.loc[npi, cluster3[0]]
            max_score = max(cluster1_score, cluster2_score, cluster3_score)

            candidates = []
            if cluster1_score == max_score:
                candidates.append(1)
            if cluster2_score == max_score:
                candidates.append(2)
            if cluster3_score == max_score:
                candidates.append(3)

            chosen_cluster = random.choice(candidates)
            if chosen_cluster == 1:
                cluster1.append(npi)
            elif chosen_cluster == 2:
                cluster2.append(npi)
            else:
                cluster3.append(npi)
    return cluster1, cluster2, cluster3

def generate_normalized_jaccard_index_matrix():
    hist1df = pd.read_excel('hist1region.xlsx')
    hist1df.set_index('Window', inplace=True)
    #filter out np with all 0
    hist1df = hist1df.loc[:,~(hist1df==0).all(axis=0)]
    nps = hist1df.columns
    nomralized_jaccard_index_matrix = pd.DataFrame(index=nps, columns=nps, dtype=float)

    for npi in nps:
        for npj in nps:
            if npi == npj:
                nomralized_jaccard_index_matrix.loc[npi, npj] = 1
            else:
                m11 = ((hist1df[npi]==1)&(hist1df[npj]==1)).sum()
                a = (hist1df[npi]==1).sum()
                b = (hist1df[npj]==1).sum()

                if a == 0 or b == 0:
                    nomralized_jaccard_index_matrix.loc[npi, npj] = np.nan
                else:
                    nomralized_jaccard_index_matrix.loc[npi, npj] = m11/min(a,b)

    nomralized_jaccard_index_matrix.to_excel('nomralized_jaccard_index_matrix.xlsx')
    return nomralized_jaccard_index_matrix, nps, hist1df

def kclustering(num_iterations = 100):
    nomralized_jaccard_index_matrix, nps, hist1df = generate_normalized_jaccard_index_matrix()
    best_clusters = None
    lowest_total_variation = float('inf')

    for i in range(num_iterations):
        print(f"Iteration {i+1}/{num_iterations}")
        #choose 3 np as random, random in python like this is true random
        random_nps = random.sample(list(nps), 3)

        cluster1_center = random_nps[0]
        cluster2_center = random_nps[1]
        cluster3_center = random_nps[2]

        cluster1 = [cluster1_center]
        cluster2 = [cluster2_center]
        cluster3 = [cluster3_center]

        new_cluster1_center = cluster1_center
        new_cluster2_center = cluster2_center
        new_cluster3_center = cluster3_center

        while True:
            cluster1 = [cluster1_center]
            cluster2 = [cluster2_center]
            cluster3 = [cluster3_center]
            cluster1, cluster2, cluster3 = assign_cluster_modify(nomralized_jaccard_index_matrix, nps, cluster1, cluster2, cluster3)
            if len(cluster1) > 1:
                max_avg = 0
                new_cluster1_center = cluster1_center
                for npi in cluster1:
                    sum = 0
                    for npj in cluster1:
                        if npi != npj:
                            sum += nomralized_jaccard_index_matrix.loc[npi, npj]
                    avg = sum/(len(cluster1)-1)
                    if avg > max_avg:
                        max_avg = avg
                        new_cluster1_center = npi

            if len(cluster2) > 1:
                max_avg = 0
                new_cluster2_center = cluster2_center
                for npi in cluster2:
                    sum = 0
                    for npj in cluster2:
                        if npi != npj:
                            sum += nomralized_jaccard_index_matrix.loc[npi, npj]
                    avg = sum/(len(cluster2)-1)
                    if avg > max_avg:
                        max_avg = avg
                        new_cluster2_center = npi

            if len(cluster3) > 1:
                max_avg = 0
                new_cluster3_center = cluster3_center
                for npi in cluster3:
                    sum = 0
                    for npj in cluster3:
                        if npi != npj:
                            sum += nomralized_jaccard_index_matrix.loc[npi, npj]
                    avg = sum/(len(cluster3)-1)
                    if avg > max_avg:
                        max_avg = avg
                        new_cluster3_center = npi

            if cluster1_center == new_cluster1_center and cluster2_center == new_cluster2_center and cluster3_center == new_cluster3_center:
                break
            else:
                cluster1_center = new_cluster1_center
                cluster2_center = new_cluster2_center
                cluster3_center = new_cluster3_center
                cluster1.clear()
                cluster2.clear()
                cluster3.clear()

        cluster1_var = assess_individual_cluster(cluster1, nomralized_jaccard_index_matrix, cluster1_center)
        cluster2_var = assess_individual_cluster(cluster2, nomralized_jaccard_index_matrix, cluster2_center)
        cluster3_var = assess_individual_cluster(cluster3, nomralized_jaccard_index_matrix, cluster3_center)
        total_variation = cluster1_var + cluster2_var + cluster3_var
        print(f"Iteration {i+1} variation: {total_variation}")
        if total_variation < lowest_total_variation:
            lowest_total_variation = total_variation
            best_clusters = (cluster1, cluster2, cluster3)

    print(f"Lowest total variation: {lowest_total_variation}")

    return best_clusters, hist1df

def get_percentage_for_1_feature_of_one_cluster(featuredf, feature_name, cluster, hist1df):
    feature_series = featuredf[feature_name]
    percentage_of_window_in_np = []
    for nps in cluster:
        np_presence_in_window = hist1df[nps]
        count11 = ((np_presence_in_window==1)&(feature_series!=0)).sum()
        countnp1 = (np_presence_in_window==1).sum()
        hist1_window_percentage = count11/countnp1
        percentage_of_window_in_np.append(hist1_window_percentage)

    return np.mean(percentage_of_window_in_np)

def get_percentage_for_all_feature_of_one_cluster(featuredf, cluster, hist1df):
    feature_percentage = {}
    for feature_name in featuredf.columns:
        feature_percentage[feature_name] = get_percentage_for_1_feature_of_one_cluster(featuredf, feature_name, cluster, hist1df)

    return feature_percentage

def feature_selection_data():
    featuredf = pd.read_csv('Hist1_region_features.csv')
    featuredf.set_index('name', inplace=True)
    best_clusters, hist1df = kclustering(100)
    cluster1, cluster2, cluster3 = best_clusters
    cluster1_feature_percentage = get_percentage_for_all_feature_of_one_cluster(featuredf, cluster1, hist1df)
    cluster2_feature_percentage = get_percentage_for_all_feature_of_one_cluster(featuredf, cluster2, hist1df)
    cluster3_feature_percentage = get_percentage_for_all_feature_of_one_cluster(featuredf, cluster3, hist1df)

    return cluster1_feature_percentage, cluster2_feature_percentage, cluster3_feature_percentage

def plot_a_cluster_radar_chart(cluster_feature_percentage, cluster_name, color):
    labels = list(cluster_feature_percentage.keys())
    values = list(cluster_feature_percentage.values())

    num_vars = len(labels)
    #Generates num_vars evenly spaced values between 0 and 2pi
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # The plot is a circle, so we need to "complete the loop" and append the start value to the end
    values += values[:1]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    plt.tight_layout(pad=10.0)

    ax.fill(angles, values, color=color, alpha=0.4)
    ax.plot(angles, values, color=color, linewidth=2)

    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])

    ax.set_xticklabels(labels, fontsize=11)
    # Move labels outward
    ax.tick_params(axis='x', pad=20)  # Increase padding between axis and labels

    plt.title(f'Radar chart for {cluster_name}')
    plt.savefig(f'{cluster_name}_radar_chart.png')
    plt.clf()

def plot_overlay_data_chart(cluster1, cluster2, cluster3, labels, cluster_names, colors ):
    num_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    plt.tight_layout(pad=10.0)

    for values, name, color in zip([cluster1, cluster2, cluster3], cluster_names, colors):
        values += values[:1]
        ax.fill(angles, values, color=color, alpha=0.4, label=name)
        ax.plot(angles, values, color=color, linewidth=2)

    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=11)
    ax.tick_params(axis='x', pad=20)
    
    plt.title('Radar Chart for all Clusters Overlay')
    plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.1))
    plt.savefig('overlay_radar_chart.png')
    plt.clf()

def plot_radar_charts():
    cluster1_feature_percentage, cluster2_feature_percentage, cluster3_feature_percentage = feature_selection_data()
    plot_a_cluster_radar_chart(cluster1_feature_percentage, 'Cluster 1', 'r')
    plot_a_cluster_radar_chart(cluster2_feature_percentage, 'Cluster 2', 'b')
    plot_a_cluster_radar_chart(cluster3_feature_percentage, 'Cluster 3', 'g')

    labels = list(cluster1_feature_percentage.keys())
    cluster1 = list(cluster1_feature_percentage.values())
    cluster2 = list(cluster2_feature_percentage.values())
    cluster3 = list(cluster3_feature_percentage.values())
    cluster_names = ['Cluster 1', 'Cluster 2', 'Cluster 3']
    colors = ['r', 'b', 'g']

    plot_overlay_data_chart(cluster1, cluster2, cluster3, labels, cluster_names, colors)



if __name__ == "__main__":
    plot_radar_charts()
    