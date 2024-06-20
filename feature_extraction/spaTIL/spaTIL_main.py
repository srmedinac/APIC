import os
from os.path import join, isfile
import cv2
import numpy as np
from extract_til_features import extract_til_features
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

# paths
patches_dir = r""  # patches dir
epi_stroma_masks_dir = r""  # nuclei mask
nuclei_masks_dir = r""  # nuclei mask
lymphocyte_masks_dir = r""  # nuclei mask
results_features_dir = r""  # results folder
draw_option = 0
histoqc_mask = np.ones((2048, 2048))  # replace with HistoQC mask if available
bar = tqdm(os.listdir(patches_dir))

for wsi in bar:
    bar.set_description(f"Processing: {wsi}")
    os.makedirs(join(results_features_dir, wsi), exist_ok=True)
    for tile in os.listdir(join(patches_dir, wsi, "tiles")):
        epi_stroma_mask_path = join(
            epi_stroma_masks_dir, wsi, tile.replace(".jpeg", ".png")
        )
        nuclei_mask_path = join(nuclei_masks_dir, wsi, tile.replace(".jpeg", ".png"))
        lymphocyte_mask_path = join(
            lymphocyte_masks_dir, wsi, tile.replace(".jpeg", ".png")
        )

        if (
            isfile(epi_stroma_mask_path)
            and isfile(nuclei_mask_path)
            and isfile(lymphocyte_mask_path)
        ):
            # read patches, epi/stroma mask, nuclei mask, and lymphocyte nuclei mask
            image = cv2.imread(join(patches_dir, wsi, "tiles", tile))
            image = image.astype(np.float64) / 255.0
            epi_stroma_mask = (
                cv2.imread(epi_stroma_mask_path, cv2.IMREAD_GRAYSCALE).astype(
                    np.float64
                )
            ) / 255.0
            nuclei_mask = cv2.imread(nuclei_mask_path, cv2.IMREAD_GRAYSCALE).astype(
                np.float64
            )
            lymphocyte_mask = cv2.imread(
                lymphocyte_mask_path, cv2.IMREAD_GRAYSCALE
            ).astype(
                np.float64
            )  # == 3
            epi_mask = (epi_stroma_mask * histoqc_mask).astype(np.float64)
            stroma_mask = ((1 - epi_stroma_mask) * histoqc_mask).astype(np.float64)
            nuclei_mask[(nuclei_mask != 0) & (nuclei_mask != 128)] = 1
            lymphocyte_mask = ((lymphocyte_mask == 128) * histoqc_mask).astype(
                np.float64
            )  # hovernet returns lymphocytes with values 128, confirm that this is the case for your data
            alpha = [0.56, 0.56]
            r = 0.07
            features, flag = extract_til_features(
                image,
                nuclei_mask,
                lymphocyte_mask,
                epi_mask,
                stroma_mask,
                histoqc_mask,
                draw_option,
                tile,
                alpha,
                r,
            )

            numNaN = np.isnan(features).sum()
            num0 = np.sum(features == 0)
            numf = np.sum((features != 0) & (~np.isnan(features)))
            if flag == 1:
                filename = tile.split(".jpeg")[0]
                np.savetxt(
                    join(results_features_dir, wsi, filename + ".csv"),
                    features,
                    delimiter=",",
                )
