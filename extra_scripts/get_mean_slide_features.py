import numpy as np
import pandas as pd
import os
from tqdm import tqdm

results_path = "/home/smedin7/nfs/RTOG_0521/RTOG_0521_spaTIL_features_fixed_to_immune_cells"
slides = os.listdir(results_path)
save_path = "/home/smedin7/nfs/RTOG_0521/RTOG_0521_spaTIL_results_mean_per_slide_immune_cells"



def get_statistics(results_path, slide):
    csv_files = os.listdir(os.path.join(results_path, slide))
    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(os.path.join(results_path, slide, csv_file))
        if df.shape[0] != 349:
            print(f"Shape of {csv_file} is {df.shape}")
            continue
        dfs.append(df)
    try:    
        df = pd.concat(dfs, axis=1, ignore_index=True)
        # calculate the mean of each row and return a df with 1 column and 349 rows
        df = df.mean(axis=1)
    except ValueError:
        print(f"Error in {slide}")
        return np.zeros((349,))

    #concatenated by statistics: e.g. sum, mean, std, min, max, skewness, variance, median of 1st feature, 2nd feature, 3rd feature, ...
    features = df.values

    return features

for slide in tqdm(slides):
    print(slide)
    stats = get_statistics(results_path, slide)
    np.save(os.path.join(save_path,f"{slide}.npy"), stats)

