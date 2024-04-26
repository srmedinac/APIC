import numpy as np
import pandas as pd
import os

patient_feats_df = pd.read_csv("/home/smedin7/nfs/RTOG_0521/RTOG_0521_features_per_patient_immune_cells_mean_349/RTOG_0521_all_patient_features_immune_cells_349.csv")
clinical_data_df = pd.read_csv("/home/smedin7/nfs/RTOG_0521/survival_analysis/RTOG_0521_clinical_data.csv")
#rename unnamed:0 to cn_deidentified
patient_feats_df.rename(columns={'Unnamed: 0': 'cn_deidentified'}, inplace=True)
print(patient_feats_df)
print(clinical_data_df)

# merge the two dataframes
df_with_features_and_clinical = pd.merge(patient_feats_df, clinical_data_df, on='cn_deidentified')
# df_with_features_and_clinical.to_csv("/home/smedin7/nfs/RTOG_0521/RTOG_0521_features_and_clinical_data.csv")
# choose the 2793 first columns
first_2793 = df_with_features_and_clinical.iloc[:, :350]

# Select the last two columns
last_two = df_with_features_and_clinical.iloc[:, 378:]
last_two.drop('cause_of_death', axis=1, inplace=True)
rx = df_with_features_and_clinical['RX']
# Concatenate the two selections together along the columns
df_with_features_bcr_time = pd.concat([first_2793, last_two, rx], axis=1)
# Save the new dataframe to a csv file
df_with_features_bcr_time.to_csv("/home/smedin7/nfs/RTOG_0521/survival_analysis/RTOG_0521_features_immune_cells_349feats_all_events_time.csv")
print(df_with_features_bcr_time)
