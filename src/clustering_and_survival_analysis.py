import pickle
import numpy as np
import pandas as pd
import math
from sklearn.cluster import KMeans
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test
import matplotlib.pyplot as plt
import seaborn as sns

def clustering_and_survival_analysis():

    # Designate your path
    ROOT_PATH = '.'
    DIM = 2048

    def cluster_samples_latent(metabric_latent, tcga_latent, cluster_number):
        
        # Function to cluster samples 
        # Input 1: METABRIC latent features
        # Input 2: TCGA BRCA latent features
        # Input 3: cluster number
        
        # Output: cluster labels 
        
        kmeans = KMeans(n_clusters=cluster_number, random_state=42)
        cluster_metabric = kmeans.fit_predict(np.array(metabric_latent))
        cluster_tcga = kmeans.predict(tcga_latent)

        return cluster_metabric, cluster_tcga

    def survival_analysis_rfs(cluster_number, dataset, meta_dict, sample_list, cluster_list):

        # Function to implement Kaplan-Meier survival analysis
        # Input 1: cluster number
        # Input 2: dataset name ("METABRIC" or "TCGA")
        # Input 3: meta data (recurrence free survival status and months, dictionary type)
        # Input 4: sample list
        # Input 5: cluster labels
        
        # Output: multivariate log-rank test p-value

        if dataset == "METABRIC":
            true = '1:Recurred'
            false = '0:Not Recurred'
        elif dataset == "TCGA":
            true = '1:Recurred/Progressed'
            false = '0:DiseaseFree'
            sample_list = [sample[0:12] for sample in sample_list]

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

        E = np.array(events, dtype=np.int32)
        T = np.array(times, dtype=np.float32)
        G = np.array(groups, dtype=np.int32)

        kmf = KaplanMeierFitter()
        fig, ax = plt.subplots(figsize=(5, 5))
        color = ['green', 'red', 'purple', 'orange', 'blue']
        lw = 2

        for i in set(cluster_list):
            ix = np.array(groups) == (i + 1)
            kmf.fit(T[ix], event_observed=E[ix], label=labels[i])
            kmf.plot(ax=ax, ci_show=False, linewidth=lw, style="-", color=color[i])
            print(kmf.median_survival_time_)

        result = multivariate_logrank_test(T, G, E)
        pvalue = result.p_value
        print(pvalue)

        result_pair = pairwise_logrank_test(T, G, E)
        print(result_pair.p_value)

        plt.xticks([0, 50, 100, 150, 200, 250])
        plt.yticks([0.00, 0.25, 0.50, 0.75, 1.00])
        plt.xlim([0, 250])
        plt.ylim([0, 1.1])
        plt.legend(["BPS-LumA", "WPS-LumA"])
        plt.xlabel("months")
        plt.ylabel("Survival probability")
        plt.grid()
        plt.show()
        
        return pvalue


    def main():

        with open(ROOT_PATH + "/metabric/metabric_" + str(DIM) + ".pickle", "rb") as f:
            metabric_latent = pickle.load(f)

        with open(ROOT_PATH + "/tcga/tcga_" + str(DIM) + ".pickle", "rb") as f:
            tcga_latent = pickle.load(f)

        with open(ROOT_PATH + "/metabric/metabric_meta_dict.pickle", "rb") as f:
            metabric_meta_dict = pickle.load(f)

        with open(ROOT_PATH+"/tcga/tcga_meta_dictionary.pickle", "rb") as f:
            tcga_meta_dict = pickle.load(f)

        with open(ROOT_PATH + "/metabric/samples_list.pickle", "rb") as f:
            metabric_samples = pickle.load(f)

        with open(ROOT_PATH + "/tcga/samples_list.pickle", "rb") as f:
            tcga_samples = pickle.load(f)

        cluster_metabric, cluster_tcga = cluster_samples_latent(metabric_latent, tcga_latent, 2)
        p_meta = survival_analysis_rfs(2, "METABRIC", metabric_meta_dict, metabric_samples, cluster_metabric)
        p_tcga = survival_analysis_rfs(2, "TCGA", tcga_meta_dict, tcga_samples, cluster_tcga)
        
        print(p_meta)
        print(p_tcga)

    main()

clustering_and_survival_analysis()
