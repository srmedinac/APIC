import numpy as np
import pandas as pd
import os

patient_feats_df = pd.read_csv("features.csv")
clinical_data_df = pd.read_csv("clinical_data.csv")
# rename unnamed:0 to cn_deidentified
patient_feats_df.rename(columns={"Unnamed: 0": "cn_deidentified"}, inplace=True)
print(patient_feats_df)
print(clinical_data_df)

# merge the two dataframes
df_with_features_and_clinical = pd.merge(
    patient_feats_df, clinical_data_df, on="cn_deidentified"
)
first_2793 = df_with_features_and_clinical.iloc[:, :350]  # choose where to join

# Select the last two columns
last_two = df_with_features_and_clinical.iloc[:, 378:]
last_two.drop("cause_of_death", axis=1, inplace=True)
rx = df_with_features_and_clinical["RX"]
df_with_features_bcr_time = pd.concat([first_2793, last_two, rx], axis=1)
df_with_features_bcr_time.to_csv("features_all_events_time.csv")
print(df_with_features_bcr_time)
