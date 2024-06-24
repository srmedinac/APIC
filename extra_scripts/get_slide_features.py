import numpy as np
import pandas as pd
import os
from tqdm import tqdm

results_path = ""
slides = os.listdir(results_path)
save_path = ""


def calculate_row_statistics(row):
    std = np.nanstd(row)
    sum_ = np.nansum(row)
    min_ = np.nanmin(row)
    max_ = np.nanmax(row)
    skewness = row.skew()
    var = np.nanvar(row)
    mean = np.nanmean(row)
    median = np.nanmedian(row)
    return pd.Series(
        {
            "Standard Deviation": std,
            "Sum": sum_,
            "Min": min_,
            "Max": max_,
            "Skewness": skewness,
            "Variance": var,
            "Mean": mean,
            "Median": median,
        }
    )


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
    except ValueError:
        print(f"Error in {slide}")
        return np.zeros((2792,))
    statistics_per_row = df.apply(calculate_row_statistics, axis=1)
    features = np.concatenate(statistics_per_row.values, axis=0)

    return features


for slide in tqdm(slides):
    print(slide)
    stats = get_statistics(results_path, slide)
    print(stats.shape)
    np.save(os.path.join(save_path, f"{slide}.npy"), stats)
