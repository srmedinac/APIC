import numpy as np
from get_nuclei_features import get_nuclei_features
from get_spaTIL_features import get_spaTIL_features
from ESW_maker2 import ESW_maker2
from ROImaker import ROImaker
from SA_drawGraphsAndConvexHull_all import SA_drawGraphsAndConvexHull_all
from drawNucContoursByClass_SA2 import drawNucContoursByClass_SA2
import matplotlib.pyplot as plt
import scipy.io


def extract_til_features(
    image,
    nuclei_mask,
    lymphocyte_mask,
    epi_mask,
    stroma_mask,
    histoqc_mask,
    draw_option,
    filename,
    alpha,
    r,
):
    image_2 = image.copy()
    nuclei_centroids, n_1, n_2 = get_nuclei_features(image, nuclei_mask)

    nuclei_centroids_rounded = np.round(nuclei_centroids)

    epi_nuclei = np.zeros(len(nuclei_centroids_rounded), dtype=bool)
    tam_nuclei = np.zeros(len(nuclei_centroids_rounded), dtype=bool)
    features = np.array([])
    tam_nuclei_count = 0
    flag = 1

    for c in range(len(nuclei_centroids_rounded)):
        epi_nuclei[c] = epi_mask[
            int(nuclei_centroids_rounded[c, 0]), int(nuclei_centroids_rounded[c, 1])
        ]
        tam_nuclei[c] = lymphocyte_mask[
            int(nuclei_centroids_rounded[c, 0]), int(nuclei_centroids_rounded[c, 1])
        ]
        if (
            lymphocyte_mask[
                int(nuclei_centroids_rounded[c, 1]), int(nuclei_centroids_rounded[c, 0])
            ]
            == 1
        ):
            tam_nuclei_count += 1

    if len(nuclei_centroids_rounded) < 3:
        features = np.zeros(350)
    else:
        for c in range(len(nuclei_centroids_rounded)):
            epi_nuclei[c] = epi_mask[
                int(nuclei_centroids_rounded[c, 0]), int(nuclei_centroids_rounded[c, 1])
            ]
            tam_nuclei[c] = lymphocyte_mask[
                int(nuclei_centroids_rounded[c, 0]), int(nuclei_centroids_rounded[c, 1])
            ]

        coords = [nuclei_centroids[~tam_nuclei, :], nuclei_centroids[tam_nuclei, :]]
        features, feat_names = get_spaTIL_features(coords, alpha, r)
        # print(feat_names)
        if draw_option == 1:
            classes = np.zeros(len(nuclei_centroids_rounded))
            classes[~tam_nuclei] = 1
            classes[tam_nuclei] = 2

            colors = [[0.500, 0.8250, 0.00], [0.3010, 0.50, 0.330]]
            colors = [(0, 0, 1, 0.5), (0, 1, 0, 0.5)]
            V30 = ESW_maker2(epi_mask, stroma_mask, histoqc_mask)
            # TODO: hack to make it work without histoqc mask and epi/stroma masks
            V41 = ROImaker(image, histoqc_mask)

            drawNucContoursByClass_SA2(
                nuclei_mask, image_2, nuclei_centroids, classes, colors
            )
            SA_drawGraphsAndConvexHull_all(image_2, V30, V41, coords, colors, r, alpha)

    return features, flag
