import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns

#load and clean data
hist1df = pd.read_excel('hist1region.xlsx')
hist1df.set_index('Window', inplace=True)
hist1df = hist1df.loc[:,~(hist1df==0).all(axis=0)]

def get_detection_frequency():
    hist1df['Count_1s'] = (hist1df==1).sum(axis=1) # count number of 1s in each row using axis=1, which is a row wised operation
    num_cols = hist1df.shape[1] -1 # subtract 1 because count_1s is added as a column
    hist1df['Detection_Frequency'] = hist1df['Count_1s'] / num_cols #Calculate detection frequency

    detection_freq = hist1df['Detection_Frequency'] # get a detection frequency series

    #drop created columns to return hist1df to its original state
    hist1df.drop(['Count_1s'], axis=1, inplace=True) # drop the count_1s column
    hist1df.drop(['Detection_Frequency'], axis=1, inplace=True) # drop the detection frequency column 

    return detection_freq

def get_co_seg():
    num_cols = hist1df.shape[1]
    #using matrix multiplication to calculate the co-segregation matrix
    #hist1df.T is the transpose of hist1df
    #hist1df.dot(hist1df.T) is the matrix multiplication of hist1df and its transpose
    #which means each row of hist1df is multiplied by each column of hist1df.T or in other words, another row of hist1df 
    co_seg_matrix = hist1df.dot(hist1df.T) / num_cols
    return co_seg_matrix

def get_linkage(detection_freq_df,co_seg_df):
    for windowA in co_seg_df.index:
        for windowB in co_seg_df.columns:
            if windowA == windowB:
                continue
            linkage = co_seg_df.loc[windowA,windowB] - detection_freq_df.loc[windowA] * detection_freq_df.loc[windowB]
            co_seg_df.loc[windowA,windowB] = linkage
    return co_seg_df

def get_normalize_linkage(linkage_df, detection_freq_df):
    for windowA in linkage_df.index:
        for windowB in linkage_df.columns:
            if windowA == windowB:
                linkage_df.loc[windowA,windowB] = np.nan
            d = linkage_df.loc[windowA,windowB]
            if d <0:
                d_max = min((detection_freq_df.loc[windowA])*(detection_freq_df.loc[windowB]), (1-detection_freq_df.loc[windowA])*(1-detection_freq_df.loc[windowB]))
                linkage_df.loc[windowA,windowB] = d/d_max
            elif d > 0:
                d_max = min(detection_freq_df.loc[windowB]*(1-detection_freq_df.loc[windowA]), (1-detection_freq_df.loc[windowB])*detection_freq_df.loc[windowA])
                linkage_df.loc[windowA,windowB] = d/d_max
            else:
                continue
    return linkage_df

def plot_heatmap(df):
    plt.figure(figsize=(10,10))
    sns.heatmap(df, annot=False, center=0, vmin=-1, vmax=1, cmap='coolwarm')
    plt.xticks([])
    plt.yticks([])
    plt.title('Normalized Linkage Heatmap')
    plt.savefig('normalize_linkage_heatmap.png')

if __name__ == "__main__":
    detection_freq_df = get_detection_frequency()
    co_seg_df = get_co_seg()
    linkage_df = get_linkage(detection_freq_df,co_seg_df)
    normalize_linkage_df = get_normalize_linkage(linkage_df, detection_freq_df)
    normalize_linkage_df.to_csv('normalize_linkage.csv')
    plot_heatmap(normalize_linkage_df)