# Identifying prognostic subgroups of luminal-A breast cancer using deep autoencoders and gene expressions
In this study, we discovered two prognostic subgroups of luminal-A breast cancer (BPS-LumA and WPS-LumA) using deep autoencoders and gene expressions. The deep autoencoders were trained using gene expression profiles of 679 luminal-A breast cancer samples in the METABRIC dataset. Then, latent features of each samples generated from the deep autoencoders were used for K-Means clustering to divide the samples into two subgroups, and Kaplan-Meier survival analysis was performed to compare prognosis (recurrence-free survival) between them. As a result, the prognosis between the two subgroups were significantly different (p-value=6.70E-05; log-rank test). This prognostic difference between two subgroups was validated using gene expression profiles of 415 luminal-A breast cancer samples in the TCGA BRCA dataset (p-value=0.004; log-rank test). Notably, the latent features were superior to the gene expression profiles and traditional dimensionality reduction method in terms of discovering the prognostic subgroups. Lastly, we discovered that ribosome-related biological functions could be potentially associated with the prognostic difference between them. Our stratification method can be contributed to understanding a complexity of luminal-A breast cancer and providing a personalized medicine.

## 1. Dataset
#### 1) METABRIC (Pereira, Bernard, et al. "The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes." Nature communications 7.1 (2016): 1-16.)
#### 2) TCGA BRCA (Ciriello, Giovanni, et al. "Comprehensive molecular portraits of invasive lobular breast cancer." Cell 163.2 (2015): 506-519.)
- All data (gene expression profiles (median Z-scores), recurrence free survival status and months, PAM50 subtype) are availble at original publications and cBioPortal (https://www.cbioportal.org/).

## 2. Source codes
- renormalizing_datasets.py
- parsing_meta_data.py
- training_autoencoder.py
- clustering_and_survival_analysis.py

## 3. Python and library versions
- python==3.7.1
- numpy==0.21.2
- pandas==1.3.3
- tensorflow==2.3.0
- scikit-learn==0.23.2
- scipy==1.7.1
- lifelines==0.24.1


