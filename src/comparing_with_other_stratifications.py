import pickle
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans

def comparing_previous_works():

    # Designate your path
    ROOT_PATH = '.'
    DIM = 2048

    def cluster_samples_latent(metabric_latent, tcga_latent, cluster_number):

        # Function to cluster samples using latent features
        # Input 1: METABRIC latent features
        # Input 2: cluster number

        # Output: cluster labels

        kmeans = KMeans(n_clusters=cluster_number, random_state=42)
        cluster_metabric = kmeans.fit_predict(np.array(metabric_latent))
        cluster_tcga = kmeans.predict(np.array(tcga_latent))

        return cluster_tcga

    def implement_chi_square(sample_list, cluster_list):

        df_bmc2016 = pd.read_excel(ROOT_PATH + "/dataset/bmc2016/13058_2016_724_MOESM2_ESM.xlsx")
        bmc2016_dict = {}
        for _, row in df_bmc2016.iterrows():
            bmc2016_dict[row["Sample ID"][0:12]] = row["LumA-R1/2"]

        o1 = 0
        t1 = 0
        o2 = 0
        t2 = 0

        for i in range(0, len(sample_list)):
            sample = sample_list[i]
            if cluster_list[i] == 0:
                if sample[0:12] in bmc2016_dict.keys():
                    if bmc2016_dict[sample[0:12]] == "LumA-R1":
                        o1 += 1
                    elif bmc2016_dict[sample[0:12]] == "LumA-R2":
                        t1 += 1
            else:
                if sample[0:12] in bmc2016_dict.keys():
                    if bmc2016_dict[sample[0:12]] == "LumA-R1":
                        o2 += 1
                    elif bmc2016_dict[sample[0:12]] == "LumA-R2":
                        t2 += 1

        obs = np.array([[o1, t1], [o2, t2]])
        print(chi2_contingency(obs))
        print(o1, t1, o2, t2)


        df_npj2019 = pd.read_excel(ROOT_PATH + "/dataset/npj2019/41523_2019_116_MOESM3_ESM.xlsx")
        npj2019_dict = {}
        for _, row in df_npj2019.iterrows():
            sample = str(row["Samples"]).replace(".", "-")
            npj2019_dict[sample] = row["Heterocellular Subtypes (High confidence and Mixed/Low confidence)"]

        g1 = 0
        i1 = 0
        s1 = 0
        t1 = 0
        e1 = 0
        g2 = 0
        i2 = 0
        s2 = 0
        t2 = 0
        e2 = 0

        for i in range(0, len(sample_list)):
            sample = sample_list[i]
            if cluster_list[i] == 0:
                if sample[0:12] in npj2019_dict.keys():
                    if npj2019_dict[sample[0:12]] == "Goblet.like":
                        g1 += 1
                    elif npj2019_dict[sample[0:12]] == "Inflammatory":
                        i1 += 1
                    elif npj2019_dict[sample[0:12]] == "Stem.like":
                        s1 += 1
                    elif npj2019_dict[sample[0:12]] == "TA":
                        t1 += 1
                    elif npj2019_dict[sample[0:12]] == "Enterocyte":
                        e1 += 1
            else:
                if sample[0:12] in npj2019_dict.keys():
                    if npj2019_dict[sample[0:12]] == "Goblet.like":
                        g2 += 1
                    elif npj2019_dict[sample[0:12]] == "Inflammatory":
                        i2 += 1
                    elif npj2019_dict[sample[0:12]] == "Stem.like":
                        s2 += 1
                    elif npj2019_dict[sample[0:12]] == "TA":
                        t2 += 1
                    elif npj2019_dict[sample[0:12]] == "Enterocyte":
                        e1 += 1

        obs = np.array([[g1, i1, s1, t1, e1], [g2, i2, s2, t2, e2]])
        print(chi2_contingency(obs))
        print(g1, i1, s1, t1, e1, g2, i2, s2, t2, e2)


    def main():

        with open(ROOT_PATH + "/metabric/metabric_" + str(DIM) + ".pickle", "rb") as f:
            metabric_latent = pickle.load(f)
            
        with open(ROOT_PATH + "/tcga/tcga_" + str(DIM) + ".pickle", "rb") as f:
            tcga_latent = pickle.load(f)

        with open(ROOT_PATH + "/tcga/samples_list.pickle", "rb") as f:
            tcga_samples = pickle.load(f)


        tcga_cluster = cluster_samples_latent(metabric_latent, tcga_latent, 2)
        implement_chi_square(tcga_samples, tcga_cluster)

    main()

comparing_previous_works()
