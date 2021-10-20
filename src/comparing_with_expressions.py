import pickle
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test

def comparing_with_expressions():

    # Designate your path
    ROOT_PATH = '.'
    DIM = 2048

    def cluster_samples_wholegenes(df_luma_metabric, cluster_number):

        # Function to cluster samples using whole gene expressions
        # Input 1: METABRIC normalized gene expressions
        # Input 2: cluster number

        # Output: cluster labels

        df_luma_metabric = df_luma_metabric.transpose()
        df_luma_metabric = df_luma_metabric.sort_index(axis=1, ascending=False)
        metabric_wholegenes = np.array(df_luma_metabric.values)

        kmeans = KMeans(n_clusters=cluster_number, random_state=42)
        cluster_metabric = kmeans.fit_predict(np.array(metabric_wholegenes))

        return cluster_metabric

    def cluster_samples_top(df_luma_metabric, cluster_number, gene_number):

        # Function to cluster samples using top variable gene expressions
        # Input 1: METABRIC normalized gene expressions
        # Input 2: cluster number

        # Output: cluster labels

        df_luma_metabric = df_luma_metabric.transpose()
        df_luma_metabric = df_luma_metabric.sort_index(axis=1, ascending=False)
        metabric_wholegenes = np.array(df_luma_metabric.values)

        mads = df_luma_metabric.mad()
        mads = pd.DataFrame(mads, columns=['mad'])
        mads = mads.sort_values(by=['mad'], ascending=False)
        genes_sorted = mads.index.to_list()[0:gene_number]

        df_top = df_luma_metabric[genes_sorted]
        metabric_top = np.array(df_top.values)

        kmeans = KMeans(n_clusters=cluster_number, random_state=42)
        cluster_metabric = kmeans.fit_predict(np.array(metabric_top))

        return cluster_metabric

    def cluster_samples_latent(metabric_latent, cluster_number):

        # Function to cluster samples using latent features
        # Input 1: METABRIC latent features
        # Input 2: cluster number

        # Output: cluster labels

        kmeans = KMeans(n_clusters=cluster_number, random_state=42)
        cluster_metabric = kmeans.fit_predict(np.array(metabric_latent))

        return cluster_metabric

    def survival_analysis_dfs(df_luma_metabric, metabric_latent, cluster_number, meta_dict, sample_list):

        true = '1:Recurred'
        false = '0:Not Recurred'

        cluster_list = cluster_samples_latent(metabric_latent, cluster_number)

        labels = ["cluster" + str(i) for i in range(1, cluster_number + 1)]
        groups = []
        events = []
        times = []

        for i in range(0, len(sample_list)):
            if meta_dict[sample_list[i]]["RFS"] == true or meta_dict[sample_list[i]]["RFS"] == false:
                rfm = float(meta_dict[sample_list[i]]["RFM"])
                if math.isnan(rfm) == False:
                    groups.append(cluster_list[i] + 1)
                    times.append(float(meta_dict[sample_list[i]]["RFM"]))
                    if meta_dict[sample_list[i]]["RFS"] == false:
                        events.append(0)
                    else:
                        events.append(1)

        print(Counter(groups))

        E = np.array(events, dtype=np.int32)
        T = np.array(times, dtype=np.float32)
        G = np.array(groups, dtype=np.int32)

        kmf = KaplanMeierFitter()

        for i in set(cluster_list):
            ix = np.array(groups) == (i + 1)
            kmf.fit(T[ix], event_observed=E[ix], label=labels[i])

        result = multivariate_logrank_test(T, G, E)
        pvalue = result.p_value
        print(pvalue)

        cluster_list = cluster_samples_wholegenes(df_luma_metabric, cluster_number)

        labels = ["cluster" + str(i) for i in range(1, cluster_number + 1)]
        groups = []
        events = []
        times = []

        for i in range(0, len(sample_list)):
            if meta_dict[sample_list[i]]["RFS"] == true or meta_dict[sample_list[i]]["RFS"] == false:
                rfm = float(meta_dict[sample_list[i]]["RFM"])
                if math.isnan(rfm) == False:
                    groups.append(cluster_list[i] + 1)
                    times.append(float(meta_dict[sample_list[i]]["RFM"]))
                    if meta_dict[sample_list[i]]["RFS"] == false:
                        events.append(0)
                    else:
                        events.append(1)

        print(Counter(groups))

        E = np.array(events, dtype=np.int32)
        T = np.array(times, dtype=np.float32)
        G = np.array(groups, dtype=np.int32)

        kmf = KaplanMeierFitter()

        for i in set(cluster_list):
            ix = np.array(groups) == (i + 1)
            kmf.fit(T[ix], event_observed=E[ix], label=labels[i])

        result = multivariate_logrank_test(T, G, E)
        pvalue = result.p_value
        print(pvalue)

        cluster_list = cluster_samples_top(df_luma_metabric, cluster_number, 5000)

        labels = ["cluster" + str(i) for i in range(1, cluster_number + 1)]
        groups = []
        events = []
        times = []

        for i in range(0, len(sample_list)):
            if meta_dict[sample_list[i]]["RFS"] == true or meta_dict[sample_list[i]]["RFS"] == false:
                rfm = float(meta_dict[sample_list[i]]["RFM"])
                if math.isnan(rfm) == False:
                    groups.append(cluster_list[i] + 1)
                    times.append(float(meta_dict[sample_list[i]]["RFM"]))
                    if meta_dict[sample_list[i]]["RFS"] == false:
                        events.append(0)
                    else:
                        events.append(1)

        print(Counter(groups))

        E = np.array(events, dtype=np.int32)
        T = np.array(times, dtype=np.float32)
        G = np.array(groups, dtype=np.int32)

        kmf = KaplanMeierFitter()

        for i in set(cluster_list):
            ix = np.array(groups) == (i + 1)
            kmf.fit(T[ix], event_observed=E[ix], label=labels[i])

        result = multivariate_logrank_test(T, G, E)
        pvalue = result.p_value
        print(pvalue)

        cluster_list = cluster_samples_top(df_luma_metabric, cluster_number, 2000)

        labels = ["cluster" + str(i) for i in range(1, cluster_number + 1)]
        groups = []
        events = []
        times = []

        for i in range(0, len(sample_list)):
            if meta_dict[sample_list[i]]["RFS"] == true or meta_dict[sample_list[i]]["RFS"] == false:
                rfm = float(meta_dict[sample_list[i]]["RFM"])
                if math.isnan(rfm) == False:
                    groups.append(cluster_list[i] + 1)
                    times.append(float(meta_dict[sample_list[i]]["RFM"]))
                    if meta_dict[sample_list[i]]["RFS"] == false:
                        events.append(0)
                    else:
                        events.append(1)

        print(Counter(groups))

        E = np.array(events, dtype=np.int32)
        T = np.array(times, dtype=np.float32)
        G = np.array(groups, dtype=np.int32)

        kmf = KaplanMeierFitter()

        for i in set(cluster_list):
            ix = np.array(groups) == (i + 1)
            kmf.fit(T[ix], event_observed=E[ix], label=labels[i])

        result = multivariate_logrank_test(T, G, E)
        pvalue = result.p_value
        print(pvalue)


    def main():
        df_luma_metabric = pd.read_csv(ROOT_PATH + "/metabric/df_luma.tsv", delimiter="\t")

        with open(ROOT_PATH + "/metabric/metabric_" + str(DIM) + ".pickle", "rb") as f:
            metabric_latent = pickle.load(f)

        with open(ROOT_PATH + "/metabric/metabric_meta_dictionary.pickle", "rb") as f:
            metabric_meta_dict = pickle.load(f)

        with open(ROOT_PATH + "/metabric/samples_list.pickle", "rb") as f:
            metabric_samples = pickle.load(f)

        survival_analysis_dfs(df_luma_metabric, metabric_latent, 2, metabric_meta_dict, metabric_samples)

    main()

comparing_with_expressions()
