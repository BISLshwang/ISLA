import pickle
import pandas as pd


def parsing_meta_data():

    # Designate your path
    ROOT_PATH = '.'

    def parse_metabric_meta(df_meta):

        # Function to recurrence-free survival status and months dictionary for METABRIC samples
        # Input: METABRIC meta data (recurence-free status, months, subtypes)
        # Output: meta data dictionaray ({"sample id":{"RFS":1, "RFM":10.1}})

        metabric_meta_dict = {}
        df_luma = df_meta[df_meta["Pam50 + Claudin-low subtype"]=="LumA"]

        for _, row in df_luma.iterrows():
            metabric_meta_dict[row["Patient ID"]] = {}
            metabric_meta_dict[row["Patient ID"]]["RFS"] = row["Relapse Free Status"]
            metabric_meta_dict[row["Patient ID"]]["RFM"] = row["Relapse Free Status (Months)"]

        return metabric_meta_dict

    def parse_tcga_meta(df_tcga_subtype, df_tcga_rfs):

        # Function to recurrence-free survival status and months dictionary for METABRIC samples
        # Input1: TCGA subtype data
        # Input2 :TCGA recurrence-free survival status and months
        # Output: meta data dictionaray ({"sample id":{"RFS":1, "RFM":10.1}})

        luma_samples = df_tcga_subtype[df_tcga_subtype["PAM50"]=="LumA"]["Case.ID"].to_list()
        df_luma = df_tcga_rfs[df_tcga_rfs["Patient ID"].isin(luma_samples)]

        tcga_meta_dict = {}
        for _, row in df_luma.iterrows():
            tcga_meta_dict[row["Patient ID"]] = {}
            tcga_meta_dict[row["Patient ID"]]["DFS"] = row["Disease Free Status"]
            tcga_meta_dict[row["Patient ID"]]["DFSM"] = row["Disease Free (Months)"]

        return tcga_meta_dict

    def main():

        df_meta = pd.read_csv(ROOT_PATH+"/metabric/brca_metabric_clinical_data.tsv", delimiter="\t")
        df_tcga_subtype = pd.read_csv(ROOT_PATH + "/tcga/1-s2.0-S0092867415011952-mmc2.csv")
        df_tcga_rfs = pd.read_csv(ROOT_PATH + "/tcga/brca_tcga_pub2015_clinical_data.tsv", delimiter="\t")

        metabric_meta_dict = parse_metabric_meta(df_meta)
        tcga_meta_dict = parse_tcga_meta(df_tcga_subtype, df_tcga_rfs)

        with open(ROOT_PATH + "/metabric/metabric_meta_dictionary.pickle", "wb") as f:
            pickle.dump(metabric_meta_dict, f)

        with open(ROOT_PATH+"/tcga/tcga_meta_dictionary.pickle", "wb") as f:
            pickle.dump(tcga_meta_dict, f)

    main()

parsing_meta_data()
