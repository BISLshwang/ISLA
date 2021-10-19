import pandas as pd
from sklearn.preprocessing import MinMaxScaler


def renormalizing_datasets():

    # Designate your path
    ROOT_PATH = '.'

    def renormalize_median_z_scores(df_luma_metabric, df_luma_tcga):

        # Function for renormalize median Z-scores such that ith gene expression value of jth sample to be in 0-1

        # Input1: METABRIC luminal-A breast cancer normalized gene expressions (row: Entrez ID, columns: sample ID)
        # Input2: TCGA luminal-A breast cancer normalized gene expresions (row: Entrez ID, columns: sample ID)

        # Output1: METABRIC luminal-A breast cancer renormalized gene expressions (row: sample ID, columns: Entrez ID)
        # Output2: TCGA luminal-A breast cancer renormalized gene expresions (row: sample ID, columns: Entrez ID)

        df_luma_metabric = df_luma_metabric.transpose()
        df_luma_metabric = df_luma_metabric.sort_index(axis=1, ascending=False)

        min_max_scaler = MinMaxScaler()
        min_max_scaler.fit(df_luma_metabric)

        luma_metabric_renormed = min_max_scaler.transform(df_luma_metabric)
        df_luma_metabric_renormed = pd.DataFrame(luma_metabric_renormed,
                                                 columns=df_luma_metabric.columns,
                                                 index=list(df_luma_metabric.index.values))

        df_luma_tcga = df_luma_tcga.transpose()
        df_luma_tcga = df_luma_tcga.sort_index(axis=1, ascending=False)

        luma_tcga_renormed = min_max_scaler.transform(df_luma_tcga)
        df_luma_tcga_renormed = pd.DataFrame(luma_tcga_renormed,
                                             columns=df_luma_tcga.columns,
                                             index=list(df_luma_tcga.index.values))

        print(df_luma_metabric_renormed.shape)
        print(df_luma_tcga_renormed.shape)

        return df_luma_metabric_renormed, df_luma_tcga_renormed


    def main():

        # Load inputs
        df_luma_metabric = pd.read_csv(ROOT_PATH+"/metabric/df_luma.tsv", delimiter="\t")
        df_luma_tcga = pd.read_csv(ROOT_PATH + "/tcga/df_luma.tsv", delimiter="\t")

        df_luma_metabric_renormed, df_luma_tcga_renormed = renormalize_median_z_scores(df_luma_metabric, df_luma_tcga)

        # Save outputs
        df_luma_metabric_renormed.to_csv(ROOT_PATH+"/metabric/df_renormed.tsv", sep="\t")
        df_luma_tcga_renormed.to_csv(ROOT_PATH + "/tcga/df_renormed.tsv", sep="\t")

    main()

renormalizing_datasets()
