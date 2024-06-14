import numpy as np
import pandas as pd
import os
features_per_slide_path = "/home/smedin7/nfs/RTOG_0521/RTOG_0521_spaTIL_results_stats_per_slide_fixed_to_immune_cells"
patient_slides_df = pd.read_csv("/home/smedin7/nfs/RTOG_0521/survival_analysis/RTOG_0521_patient_slides.csv")
patient_id_to_deid = pd.read_csv("/home/smedin7/nfs/RTOG_0521/survival_analysis/RTOG_0521_deid_to_unique_patient.csv")
save_path = "/home/smedin7/nfs/RTOG_0521/RTOG_0521_features_per_patient_fixe_to_immune_cells"
features_per_slide = os.listdir(features_per_slide_path)

slides_dict = {}
slides_dict_deid = {}
unique_usis = patient_slides_df['USI'].unique()
for usi in unique_usis:
    slides_dict[usi] = patient_slides_df.loc[patient_slides_df['USI'] == usi, 'Image ID'].tolist()
    deid = patient_id_to_deid.loc[patient_id_to_deid['USI'] == usi, 'cn_deidentified'].values[0]
    slides_dict_deid[deid] = patient_slides_df.loc[patient_slides_df['USI'] == usi, 'Image ID'].tolist()
    

slides_dict_df = pd.DataFrame.from_dict(slides_dict, orient='index')
slides_dict_df_deid = pd.DataFrame.from_dict(slides_dict_deid, orient='index')
# slides_dict_df.to_csv("/home/smedin7/nfs/RTOG_0521_patient_slides_dict.csv")
# slides_dict_df_deid.to_csv("/home/smedin7/nfs/RTOG_0521_patient_slides_dict_deid.csv")

print(slides_dict_df_deid)
patient_features_for_df = []
for index, row in slides_dict_df_deid.iterrows():
    patient_features = []
    row = row.dropna()
    for slide in row:
        slide_features = np.load(os.path.join(features_per_slide_path, f"{int(slide)}.npy"))
        patient_features.append(slide_features)
    # stack the features of all slides of a patient
    patient_features = np.stack(patient_features, axis=0)
    #average the features of all slides of a patient
    patient_features = np.mean(patient_features, axis=0)
    patient_features_for_df.append(patient_features)
    print(patient_features.shape)
    np.save(os.path.join(save_path, f"{index}.npy"), patient_features)
    # print(f"Saved {index}.npy")

# make a dataframe of patient features
patient_features_for_df = np.stack(patient_features_for_df, axis=0)
patient_feats_df = pd.DataFrame(patient_features_for_df, index=slides_dict_df_deid.index)
print(patient_feats_df)
patient_feats_df.to_csv("/home/smedin7/nfs/RTOG_0521/RTOG_0521_features_per_patient_fixe_to_immune_cells/RTOG_0521_all_patient_features_fixed_to_immune_cells.csv")