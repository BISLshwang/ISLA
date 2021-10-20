import pickle
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
import gseapy as gp

def finding_pathways_with_degs():

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

        return cluster_metabric, cluster_tcga

    def filter_degs(sample_list, cluster_list, df_luma):

        genes = df_luma["Entrez_Gene_Id"].to_list()
        
        bps_samples = []
        wps_samples = []
        for i in range(0, len(sample_list)):
            if cluster_list[i] == 0:
                bps_samples.append(sample_list[i])
            else:
                wps_samples.append(sample_list[i])

        df_bps = df_luma[bps_samples]
        df_wps = df_luma[wps_samples]

        medians = []
        median_dict = {}
        wrs_dict = {}
        sign_dict = {}
        for i in range(0, 17202):
            bps_expressions = df_bps.loc[i, :].to_list()
            wps_expressions = df_wps.loc[i, :].to_list()
            median1 = np.median(bps_expressions)
            median2 = np.median(wps_expressions)
            median_diff = median1 - median2
            medians.append(abs(median_diff))
            median_dict[genes[i]] = median_diff
            if median_diff > 0:
                sign_dict[genes[i]] = 0
            else:
                sign_dict[genes[i]] = 1
            s, p = stats.ranksums(bps_expressions, wps_expressions)
            wrs_dict[genes[i]] = p

        return median_dict, wrs_dict, sign_dict
    
    def implement_enrichr(median_dict_metabric, wrs_dict_metabric, sign_dict_metabric, median_dict_tcga, wrs_dict_tcga, sign_dict_tcga):

        degs = []
        for gene in median_dict_metabric.keys():
            if wrs_dict_metabric[gene] < 0.001 and wrs_dict_tcga[gene] < 0.001:
                if abs(median_dict_metabric[gene]) > 0.4 and abs(median_dict_tcga[gene]) > 0.353:
                    if sign_dict_metabric[gene] == sign_dict_tcga[gene]:
                        degs.append(gene)

        with open(ROOT_PATH + "/common/gsea/degs.pickle", "wb") as f:
            pickle.dump(degs, f)

        gp.enrichr(gene_list=degs, description='Enrichr',
                   gene_sets=['GO_Biological_Process_2018', 'KEGG_2019_Human', 'WikiPathways_2019_Human'])

    def main():

        df_luma_metabric = pd.read_csv(ROOT_PATH + "/metabric/df_luma.tsv", delimiter="\t")
        df_luma_tcga = pd.read_csv(ROOT_PATH + "/tcga/df_luma.tsv", delimiter="\t")

        with open(ROOT_PATH + "/metabric/metabric_" + str(DIM) + ".pickle", "rb") as f:
            metabric_latent = pickle.load(f)

        with open(ROOT_PATH + "/tcga/tcga_" + str(DIM) + ".pickle", "rb") as f:
            tcga_latent = pickle.load(f)
            
        with open(ROOT_PATH + "/metabric/samples_list.pickle", "rb") as f:
            metabric_samples = pickle.load(f)

        with open(ROOT_PATH + "/tcga/samples_list.pickle", "rb") as f:
            tcga_samples = pickle.load(f)

        cluster_metabric, cluster_tcga = cluster_samples_latent(metabric_latent, tcga_latent, 2)

        median_dict_metabric, wrs_dict_metabric, sign_dict_metabric = filter_degs(metabric_samples, cluster_metabric, df_luma_metabric)
        median_dict_tcga, wrs_dict_tcga, sign_dict_tcga = filter_degs(tcga_samples, cluster_tcga, df_luma_tcga)
        
        implement_enrichr(median_dict_metabric, wrs_dict_metabric, sign_dict_metabric, median_dict_tcga, wrs_dict_tcga, sign_dict_tcga)
        
    main()
    
finding_pathways_with_degs()
