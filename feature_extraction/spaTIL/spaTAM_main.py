"Author: Mayukhmala Jana"
import os
from os.path import join, isfile
import cv2
import numpy as np
from extract_tam_features import extract_tam_features
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
# paths
patches_dir = r"/home/smedin7/nfs/renal_project_liping/patches"  # patches dir
epi_stroma_masks_dir = (
    r"/home/smedin7/nfs/renal_project_liping/nuclei_segmentation"  # nuclei mask
)
nuclei_masks_dir = (
    r"/home/smedin7/nfs/renal_project_liping/nuclei_segmentation"  # nuclei mask
)
macrophage_masks_dir = (
    r"/home/smedin7/nfs/renal_project_liping/nuclei_segmentation"  # nuclei mask
)
results_features_dir = (
    r"/home/smedin7/nfs/renal_project_liping/spaTIL_features"  # result
)
draw_option = 0
histoqc_mask = np.ones((2048, 2048))
bar = tqdm(os.listdir(patches_dir)[1:])

for wsi in bar:
    bar.set_description(f"Processing: {wsi}")
    os.makedirs(join(results_features_dir,wsi), exist_ok=True)
    for tile in os.listdir(join(patches_dir,wsi,"tiles")):
        epi_stroma_mask_path = join(epi_stroma_masks_dir, wsi, tile.replace(".jpeg", ".png"))
        nuclei_mask_path = join(nuclei_masks_dir, wsi, tile.replace(".jpeg", ".png"))
        macrophage_mask_path = join(macrophage_masks_dir, wsi, tile.replace(".jpeg", ".png"))

        if (
            isfile(epi_stroma_mask_path)
            and isfile(nuclei_mask_path)
            and isfile(macrophage_mask_path)
        ):
            # read patches, epi/stroma mask, nuclei mask, and macrophage nuclei mask
            image = cv2.imread(join(patches_dir,wsi,"tiles",tile))
            image = image.astype(np.float64) / 255.0
            epi_stroma_mask = (
                cv2.imread(epi_stroma_mask_path, cv2.IMREAD_GRAYSCALE).astype(np.float64)
                
            ) / 255.0
            nuclei_mask = (
                cv2.imread(nuclei_mask_path, cv2.IMREAD_GRAYSCALE).astype(np.float64)
                
            )
            tam_mask = (
                cv2.imread(macrophage_mask_path, cv2.IMREAD_GRAYSCALE).astype(np.float64)
                
            )  # == 3
            epi_mask = (epi_stroma_mask * histoqc_mask).astype(np.float64)
            stroma_mask = ((1 - epi_stroma_mask) * histoqc_mask).astype(np.float64)
            #nuclei_mask = (nuclei_mask * histoqc_mask).astype(np.float64)
            nuclei_mask[(nuclei_mask != 0) & (nuclei_mask != 128)] = 1
            tam_mask = ((tam_mask == 128) * histoqc_mask).astype(np.float64)
            alpha = [0.58, 0.58]  # thresh distance betwwen nuclei
            r = 0.07
            features, flag = extract_tam_features(
                image,
                nuclei_mask,
                tam_mask,
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
            # print(numNaN, num0, numf)
            if flag == 1:
                filename = tile.split(".jpeg")[0]
                np.savetxt(
                    join(results_features_dir, wsi, filename + ".csv"), features, delimiter=","
                )
