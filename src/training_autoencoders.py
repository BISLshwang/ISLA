import pickle
import numpy as np
import pandas as pd
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model

def training_autoencoders():

    # Designate your path
    ROOT_PATH = '.'

    GENE_NUMBER = 5000
    DIM = 64
    
    def select_variable_genes():
        
        # Function to select the top 5,000 genes with the highest variability genes
        # Input: METABRIC luminal-A breast cancer normalized gene expressions (row: Entrez ID, columns: sample ID)
        # Output: The list of top 5,000 genes with the highest variability genes
        
        df_luma_metabric = pd.read_csv(ROOT_PATH+"/metabric/df_luma.tsv", delimiter="\t")
        df_luma_metabric = df_luma_metabric.transpose()
        df_luma_metabric = df_luma_metabricy.sort_index(axis=1, ascending=False)
        
        mads = df_luma_metabric.mad()
        mads = pd.DataFrame(mads, columns=['mad'])
        mads = mads.sort_values(by=['mad'], ascending=False)
        genes_sorted = mads.index.to_list()[0:5000]
        
        return genes_sorted

    def make_train_features(df_metabric_renormed, genes_sorted):

        # Function to make training input features for deep autoencoders
        # Input:  METABRIC luminal-A breast cancer renormalized gene expressions (row: sample ID, columns: Entrez ID)
        # Output: numpy array

        df_metabric_renormed = df_metabric_renormed.set_index("sample")
        df_metabric_renormed = df_metabric_renormed.sort_index(axis=1, ascending=True)
        df_metabric_renormed = df_metabric_renomred[genes_sorted]
        train_features = list(df_metabric_renormed.values)
        train_samples = list(df_metabric_renormed.index)

        return np.array(train_features), train_samples


    def make_test_features(df_tcga_renormed, genes_sorted):

        # Function to make validation input features for deep autoencoders
        # Input:  TCGA luminal-A breast cancer renormalized gene expressions (row: sample ID, columns: Entrez ID)
        # Output: numpy array

        df_tcga_renormed = df_tcga_renormed.set_index("sample")
        df_tcga_renormed = df_tcga_renormed.sort_index(axis=1, ascending=True)
        df_tcga_renormed = df_tcga_renormed[genes_sorted]
        test_features = list(df_tcga_renormed.values)
        test_samples = list(df_tcga_renormed.index)

        return np.array(test_features), test_samples

    def train_autoencoder(train_features, test_features):

        # Input: training features, validation features
        # Ouptput: Autoencoder, encoder, MSE
        # Autoencoder with 5 layers: input layer, three hidden layers (same size), output layer

        input_profile = Input(shape=(GENE_NUMBER,))
        encoded = Dense(DIM, activation="relu")(input_profile)
        encoded = Dense(DIM, activation='relu')(encoded)

        decoded = Dense(DIM, activation="relu")(encoded)
        decoded = Dense(GENE_NUMBER, activation="sigmoid")(encoded)

        autoencoder = Model(input_profile, decoded)

        encoder = Model(input_profile, encoded)

        encoded_input = Input(shape=(DIM,))
        decoder_layer = autoencoder.layers[-2]
        decoder = Model(encoded_input, decoder_layer(encoded_input))

        autoencoder.compile(optimizer="adam", loss="mean_squared_error")

        hist = autoencoder.fit(train_features, train_features,
                               epochs=100,
                               batch_size=16,
                               shuffle=True,
                               validation_data=(test_features, test_features))

        return autoencoder, encoder, hist

    def make_latent_features_metabric(encoder, train_features):

        metabric_latent = encoder.predict(train_features)
        return metabric_latent

    def make_latent_features_tcga(encoder, test_features):

        tcga_latent = encoder.predict(test_features)
        return tcga_latent


    def main():

        # Load inputs
        df_metabric_renormed = pd.read_csv(ROOT_PATH + "/metabric/df_renormed.tsv", delimiter="\t")
        df_tcga_renormed = pd.read_csv(ROOT_PATH + "/tcga/df_renormed.tsv", delimiter="\t")

        genes_sorted = select_variable_genes()
        train_features, train_samples = make_train_features(df_metabric_renormed, genes_sorted)
        test_features, test_samples = make_test_features(df_tcga_renormed, genes_sorted)
        autoencoder, encoder, hist = train_autoencoder(train_features, test_features)
        latent_features_metabric = make_latent_features_metabric(encoder, train_features)
        latent_features_tcga = make_latent_features_tcga(encoder, test_features)

        # Save outputs
        autoencoder.save(ROOT_PATH + "/autoencoder/autoencoder_" + str(DIM) + ".h5")
        autoencoder.save_weights(ROOT_PATH + "/autoencoder/autoencoder_" + str(DIM) + "_weights.h5")

        with open(ROOT_PATH + "/autoencoder/hist_loss_"+str(DIM)+".pickle", "wb") as f:
            pickle.dump(hist.history['loss'], f)

        with open(ROOT_PATH + "/autoencoder/hist_val_loss_"+str(DIM)+".pickle", "wb") as f:
            pickle.dump(hist.history['val_loss'], f)

        with open(ROOT_PATH + "/metabric/metabric_" + str(DIM) + ".pickle", "wb") as f:
            pickle.dump(latent_features_metabric, f)

        with open(ROOT_PATH + "/tcga/tcga_" + str(DIM) + ".pickle", "wb") as f:
            pickle.dump(latent_features_tcga, f)
            
        with open(ROOT_PATH + "/metabric/samples_list.pickle", "wb") as f:
            pickle.dump(train_samples, f)

        with open(ROOT_PATH + "/tcga/samples_list.pickle", "wb") as f:
            pickle.dump(test_samples, f)

    main()

training_autoencoders()
