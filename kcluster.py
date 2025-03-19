import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random

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

def makingheatmap(hist1df, cluster,cluster_name):
    clusterdf = hist1df.loc[:,cluster]
    clusterdf = clusterdf.T
    plt.figure(figsize=(11, 11))
    sns.heatmap(clusterdf, cmap=['white','red'], cbar=False)
    plt.xticks([])
    plt.yticks(rotation=0)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel('Genomic Windows', fontsize=10)
    plt.ylabel('NPs', fontsize=10)
    plt.title(f'Cluster {cluster_name} Heatmap', fontsize=12)
    plt.savefig(f'cluster_{cluster_name}.png')
    plt.clf()

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

def kclustering(num_iterations = 10):
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

    makingheatmap(hist1df, best_clusters[0], '1')
    makingheatmap(hist1df, best_clusters[1], '2')
    makingheatmap(hist1df, best_clusters[2], '3')

    return best_clusters, hist1df

def feature_selection():
    featuredf = pd.read_csv('Hist1_region_features.csv')
    featuredf.set_index('name', inplace=True)
    hist1_feature = featuredf['Hist1']
    best_clusters, hist1df = kclustering(100)
    cluster_percentage = []
    for cluster in best_clusters:
        percentagelist = []
        for nps in cluster:
            np_presence_in_window = hist1df[nps]
            count11 = ((np_presence_in_window==1)&(hist1_feature!=0)).sum()
            countnp1 = (np_presence_in_window==1).sum()
            hist1_window_percentage = count11/countnp1
            percentagelist.append(hist1_window_percentage)

        cluster_percentage.append(percentagelist)

    sns.boxplot(data=cluster_percentage, boxprops={'facecolor':'none', 'edgecolor':'black'}) 

    colors = ['red', 'blue', 'green']
    for i, percentages in enumerate(cluster_percentage):
        x = np.random.normal(i, 0.05, size=len(percentages))
        plt.scatter(x, percentages, alpha=0.7,c=colors[i%len(colors)], edgecolors='w')

    plt.title('Hist1 Percentage in Window')
    plt.savefig('boxplot_hist1.png')
    plt.clf()

    lad_features = featuredf['LAD']
    cluster_percentage.clear()
    for cluster in best_clusters:
        percentagelist = []
        for nps in cluster:
            np_presence_in_window = hist1df[nps]
            count11 = ((np_presence_in_window==1)&(lad_features!=0)).sum()
            countnp1 = (np_presence_in_window==1).sum()
            hist1_window_percentage = count11/countnp1
            percentagelist.append(hist1_window_percentage)

        cluster_percentage.append(percentagelist)

    sns.boxplot(data=cluster_percentage, boxprops={'facecolor':'none', 'edgecolor':'black'})
    colors = ['red', 'blue', 'green']
    for i, percentages in enumerate(cluster_percentage):
        x = np.random.normal(i, 0.05, size=len(percentages))
        plt.scatter(x, percentages, alpha=0.7,c=colors[i%len(colors)], edgecolors='w')


    plt.title('LAD Percentage in Window')
    plt.savefig('boxplot_lad.png')
    plt.clf()

    return best_clusters

def characterize_cluster_by_radial_pos(cluster,cluster_name):
    radial_pos_df = pd.read_csv('hist1radialpos.csv', index_col=0)
    categories = ['Strongly Apical', 'Somewhat Apical', 'Neither Apical nor Equatorial', 'Somewhat Equatorial', 'Strongly Equatorial']
    count_for_each_pos = [0] * 5
    for np in cluster:
        radial_pos = radial_pos_df.loc[np, 'Radial Position']
        if radial_pos in categories:
            index = categories.index(radial_pos)
            count_for_each_pos[index] += 1

    total = len(cluster)
    if total == 0:
        percentages = [0] * 5
    else:
        percentages = [count/total*100 for count in count_for_each_pos]

    data = pd.DataFrame({
        'Radial Position': categories,
        'Percentage': percentages
    })

    plt.figure(figsize=(12, 6))
    sns.barplot(x='Radial Position', y='Percentage', data=data)

    plt.xlabel('Radial Position')
    plt.ylabel('Percentage of NPs')
    plt.title('Radial Position of NPs in Cluster Distribution')

    plt.savefig(f'radial_pos_distribution_{cluster_name}.png')
    plt.clf()

def characterize_clusters_stacked_bar(best_clusters):
    radial_pos_df = pd.read_csv('hist1radialpos.csv', index_col=0)
    categories = ['Strongly Apical', 'Somewhat Apical', 'Neither Apical nor Equatorial', 'Somewhat Equatorial', 'Strongly Equatorial']

    cluster_percentages = []

    for cluster in best_clusters:
        count_for_each_pos = [0] * len(categories)
        
        for np in cluster:
            radial_pos = radial_pos_df.loc[np, 'Radial Position']
            if radial_pos in categories:
                index = categories.index(radial_pos)
                count_for_each_pos[index] += 1

        total = len(cluster)
        percentages = [(count / total) * 100 if total > 0 else 0 for count in count_for_each_pos]
        cluster_percentages.append(percentages)

    data = pd.DataFrame(cluster_percentages, columns=categories, index=[f"Cluster {i+1}" for i in range(len(best_clusters))])

    plt.figure(figsize=(20, 10))
    data.plot(kind="bar", stacked=True, colormap="Set2", edgecolor="black", alpha=0.8)

    plt.xlabel("Cluster")
    plt.ylabel("Percentage of NPs")
    plt.title("Stacked Bar Chart of Radial Position Distribution Across Clusters")
    plt.legend(title="Radial Position", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.xticks(rotation=0)
    plt.savefig("stacked_radial_position.png", bbox_inches="tight")
    plt.clf()


def feature_activity_3():
    best_clusters = feature_selection()
    characterize_cluster_by_radial_pos(best_clusters[0], '1')
    characterize_cluster_by_radial_pos(best_clusters[1], '2')
    characterize_cluster_by_radial_pos(best_clusters[2], '3')
    characterize_clusters_stacked_bar(best_clusters)


if __name__ == "__main__":
    feature_activity_3()
    