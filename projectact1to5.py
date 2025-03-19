import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import random

df = pd.read_csv('data.txt', sep='\t')
df['Window'] = df.iloc[:,0:3].apply(lambda x: f"{x.iloc[0]}:{x.iloc[1]}-{x.iloc[2]}", axis=1)
df.set_index('Window', inplace=True)
df.drop(df.columns[0:3], axis=1, inplace=True) # drop with inplace=True to avoid copy


windows_counts = (df.iloc[:,:] == 1).sum(axis=0)
# store number of np for each window
np_counts = (df.iloc[:,:] == 1).sum(axis=1)
np_counts_not_zero = np_counts[np_counts != 0]



def activity1():
    #print number of lines representing number of windows, read.csv() takes the first line as column name automatically
    print(f"Number of genomic windows: {df.shape[0]}")

    #print number of np = number of rows-3
    print(f"Number of NPs: {df.shape[1]}")

    

    avg_windows = windows_counts.mean()
    print(f"Average number of windows per NP: {avg_windows}")
    max_windows = windows_counts.max()
    print(f"Maximum number of windows for a NP: {max_windows}")
    min_windows = windows_counts.min()
    print(f"Minimum number of windows for a NP: {min_windows}")  

    

    avg_nps = np_counts_not_zero.mean()
    print(f"Average number of NPs per window: {avg_nps}")
    max_nps = np_counts_not_zero.max()
    print(f"Maximum number of NPs for a window: {max_nps}")
    min_nps = np_counts_not_zero.min()
    print(f"Minimum number of NPs for a window: {min_nps}")

def activity2():
    detection_freq = windows_counts*100/df.shape[0]
    detection_freq = detection_freq.rename('Detection Frequency')
    plt.figure(figsize=(15, 6))
    plt.scatter(detection_freq.index, detection_freq)
    plt.xlabel('NP')
    plt.ylabel('Detection Frequency')
    plt.title('Detection Frequency of NPs Scatter Plot')
    plt.xticks([])
    plt.savefig('detection_frequency.png')

    ascending_sorted_detection_freq = detection_freq.sort_values(ascending=True)
    print(f"Top 10 NPs with the highest detection frequency:\n{ascending_sorted_detection_freq[-10:]}")

    np_percentile_bins = [0,0.2,0.4,0.6,0.8,1.0]
    np_percentile_labels = ['Strongly Apical','Somewhat Apical','Neither Apical nor Equatorial','Somewhat Equatorial','Strongly Equatorial']
    #cut data in percentile range
    percentile_cats = pd.qcut(detection_freq, q=np_percentile_bins, labels=np_percentile_labels)
    sorted_np_estimate_rad_pos = pd.DataFrame({
        'Detection Frequency': ascending_sorted_detection_freq,
        'Estimated Radial Position': percentile_cats
    })

    np_estimate_radial_pos = pd.DataFrame({
        'Detection Frequency': detection_freq,
        'Estimated Radial Position': detection_freq.index.map(sorted_np_estimate_rad_pos['Estimated Radial Position'])
    })

    #np_estimate_radial_pos.to_csv('np_estimated_radial_pos.csv', sep=',', index=True)

    frequency_of_np_detected_each_window = np_counts_not_zero*100/(df.shape[1]-3)
    frequency_of_np_detected_each_window = frequency_of_np_detected_each_window.rename('Frequency of NP Detected Each Window')
    descending_sorted_frequency_of_np_detected_each_window = frequency_of_np_detected_each_window.sort_values(ascending=False)
    windows_percentile_bins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    windows_percentile_labels = ['10','9','8','7','6','5','4','3','2','1']
    #cut data in percentile range
    window_percentile_cats = pd.qcut(frequency_of_np_detected_each_window, q=windows_percentile_bins, labels=windows_percentile_labels)
    sorted_window_estimate_compact_rate = pd.DataFrame({
        'Frequency of NP Detected Each Window': descending_sorted_frequency_of_np_detected_each_window,
        'Compaction Rating': window_percentile_cats
    })
    window_compaction_rating = pd.DataFrame({
        'Frequency of NP Detected Each Window': frequency_of_np_detected_each_window,
        'Compaction Rating': frequency_of_np_detected_each_window.index.map(sorted_window_estimate_compact_rate['Compaction Rating'])
    })
    #window_compaction_rating.to_csv('window_compaction_rating.txt', sep='\t', index=True)

def activity3():
    df['chr'] = df.index.str.split(':').str[0] #extract chr from window
    df['start'] = df.index.str.split(':').str[1].str.split('-').str[0].astype(int) #extract start from window and convert to int
    df['end'] = df.index.str.split(':').str[1].str.split('-').str[1].astype(int) #extract end from window and convert to int

    filtered_df = df.loc[(df['chr'] == 'chr13') & (df['end'] > 21700000) & (df['start'] < 24100000)].copy()
    filtered_df.drop(columns=['chr','start','end'], inplace=True)
    filtered_df.to_excel('hist1region.xlsx')
    print(f"Number of genomic windows: {filtered_df.shape[0]}")
    hist1windows_count = (filtered_df.iloc[:,:] == 1).sum(axis=0)
    hist1np_count = (filtered_df.iloc[:,:] == 1).sum(axis=1)
    hist1windows_count = hist1windows_count[hist1windows_count != 0]
    hist1np_count = hist1np_count[hist1np_count != 0]
    print(f"Number of windows with at least one NP: {hist1windows_count.shape[0]}")


    hist1_avg_windows = hist1windows_count.mean()
    print(f"Average number of windows per NP: {hist1_avg_windows}")
    hist1_max_windows = hist1windows_count.max()
    print(f"Maximum number of windows for a NP: {hist1_max_windows}")
    hist1_min_windows = hist1windows_count.min()
    print(f"Minimum number of windows for a NP: {hist1_min_windows}")

    hist1_avg_nps = hist1np_count.mean()
    print(f"Average number of NPs per window: {hist1_avg_nps}")
    hist1_max_nps = hist1np_count.max()
    print(f"Maximum number of NPs for a window: {hist1_max_nps}")
    hist1_min_nps = hist1np_count.min()
    print(f"Minimum number of NPs for a window: {hist1_min_nps}")


    detection_freq = windows_counts*100/df.shape[0]
    ascending_sorted_detection_freq = detection_freq.sort_values(ascending=True)
    np_percentile_bins = [0,0.2,0.4,0.6,0.8,1.0]
    np_percentile_labels = ['Strongly Apical','Somewhat Apical','Neither Apical nor Equatorial','Somewhat Equatorial','Strongly Equatorial']
    #cut data in percentile range
    percentile_cats = pd.qcut(detection_freq, q=np_percentile_bins, labels=np_percentile_labels)
    sorted_np_estimate_rad_pos = pd.DataFrame({
        'Detection Frequency': ascending_sorted_detection_freq,
        'Estimated Radial Position': percentile_cats
    })

    hist1radialpos = pd.DataFrame({
        'NP': hist1windows_count.index,
        'Number of Windows': hist1windows_count,
        'Radial Position': hist1windows_count.index.map(sorted_np_estimate_rad_pos['Estimated Radial Position'])
    })
    hist1radialpos.to_csv('hist1radialpos.csv', sep=',', index=True)
    hist1radialpostcount = hist1radialpos['Radial Position'].value_counts()
    mostcommonradialpos = hist1radialpostcount.idxmax()
    mostcommoncount = hist1radialpostcount.max()
    print(f"Most common radial position: {mostcommonradialpos} ({mostcommoncount})")

    frequency_of_np_detected_each_window = np_counts_not_zero*100/(df.shape[1]-3)
    frequency_of_np_detected_each_window = frequency_of_np_detected_each_window.rename('Frequency of NP Detected Each Window')
    descending_sorted_frequency_of_np_detected_each_window = frequency_of_np_detected_each_window.sort_values(ascending=False)
    windows_percentile_bins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    windows_percentile_labels = ['10','9','8','7','6','5','4','3','2','1']
    #cut data in percentile range
    window_percentile_cats = pd.qcut(frequency_of_np_detected_each_window, q=windows_percentile_bins, labels=windows_percentile_labels)
    sorted_window_estimate_compact_rate = pd.DataFrame({
        'Frequency of NP Detected Each Window': descending_sorted_frequency_of_np_detected_each_window,
        'Compaction Rating': window_percentile_cats
    })

    hist1window_compaction_rating = pd.DataFrame({
        'Window': hist1np_count.index,
        'Number of NPs': hist1np_count,
        'Compaction Rating': hist1np_count.index.map(sorted_window_estimate_compact_rate['Compaction Rating'])
    })
    hist1window_compaction_ratingcount = hist1window_compaction_rating['Compaction Rating'].value_counts()
    mostcommoncompactionrating = hist1window_compaction_ratingcount.idxmax()
    mostcommoncount = hist1window_compaction_ratingcount.max()
    print(f"Most typical compaction rating: {mostcommoncompactionrating} ({mostcommoncount})")

def activity4():
    hist1df = pd.read_excel('hist1region.xlsx')
    hist1df.set_index('Window', inplace=True)
    #filter out np with all 0
    hist1df = hist1df.loc[:,~(hist1df==0).all(axis=0)]
    
    nps = hist1df.columns
    jaccard_index_matrix = pd.DataFrame(index=nps, columns=nps, dtype=float)

    for npi in nps:
        for npj in nps:
            if npi == npj:
                jaccard_index_matrix.loc[npi, npj] = 0
            else:
                m11 = ((hist1df[npi]==1)&(hist1df[npj]==1)).sum() #all windows where both np are present, top of jaccard index formular
                denom = ((hist1df[npi]==1)|(hist1df[npj]==1)).sum() #all windows where at least one np is present, bottom of jaccard index formular
                jacard_index = m11/denom if denom != 0 else 0
                jaccard_index_matrix.loc[npi, npj] = jacard_index

    jaccard_index_matrix.to_csv('jaccard_index_matrix.txt', sep='\t', index=True)
    jaccard_index_matrix.to_excel('jaccard_index_matrix.xlsx')

def activity5():
    jaccard_index_matrix = pd.read_excel('jaccard_index_matrix.xlsx')
    jaccard_index_matrix.set_index('Unnamed: 0', inplace=True)
    sns.heatmap(jaccard_index_matrix, cmap='plasma', annot=False)
    plt.savefig('jaccard_index_matrix.png')
    plt.clf()
    
    jaccard_distance_matrix = 1-jaccard_index_matrix
    sns.heatmap(jaccard_distance_matrix, cmap='plasma', annot=False)
    plt.savefig('jaccard_distance_matrix.png')  


if __name__ == "__main__":
    activity1()
    activity2()
    activity3()
    #activity4()
    #activity5()