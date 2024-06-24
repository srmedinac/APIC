import os
import numpy as np
import pandas as pd
import cv2
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from histoprep import SlideReader
from PIL import Image
import time

# Define paths
patch_feature_path = ""
slide_name = ""
slide_name_2 = ""
slide_path = ""
patches_path = ""
downsample = 4
patch_size = 1024
level = 1
feature = 1
patches_dir = os.path.join(patches_path, slide_name, "tiles")
patches_dir_2 = os.path.join(patches_path, slide_name_2, "tiles")


def get_feature_values(csv_dir, patches_dir, slide_name, feature_idx=1):
    feature_51_list = []
    patches = os.listdir(patches_dir)

    for csv_file in os.listdir(os.path.join(csv_dir, slide_name)):
        if csv_file.endswith(".csv"):
            csv_path = os.path.join(csv_dir, slide_name, csv_file)
            df = pd.read_csv(csv_path, header=None)
            if df.empty:
                continue
            patch_name = csv_file.replace(".csv", ".jpeg")
            if patch_name in patches:
                feature_51_list.append(df.iloc[feature_idx - 1, 0])

    return feature_51_list


def normalize_feature_values(features_img_1, features_img_2):
    # combine the two lists
    features = features_img_1 + features_img_2
    scaler = StandardScaler()
    feature_list = scaler.fit_transform(np.array(features).reshape(-1, 1))
    # transform the two separate lists
    feature_1 = scaler.transform(np.array(features_img_1).reshape(-1, 1))
    feature_2 = scaler.transform(np.array(features_img_2).reshape(-1, 1))

    return feature_1, feature_2


def make_overlay(slide_path, slide_name, patches, features, n_bins=20):
    slide = SlideReader(os.path.join(slide_path, slide_name + ".svs"))
    patch_coords = [
        (int(p.split("_")[0][1:]) // downsample, int(p.split("_")[1][1:]) // downsample)
        for p in patches
        if p.endswith(".jpeg")
    ]
    # cmap = create_binned_colormap("coolwarm", n_bins)
    cmap = cm["seismic"]
    overlay_shape = slide.level_dimensions[level][::-1] + (4,)
    overlay = np.ones(overlay_shape, dtype=np.uint8)
    print("Getting tissue mask...")
    _, tissue_mask = slide.get_tissue_mask(level=level)
    tissue_mask = tissue_mask.transpose(1, 0)
    downsampled_patch_size = patch_size // downsample
    for i, (x, y) in enumerate(patch_coords):
        if i >= len(features):
            break
        color = cmap(features[i]).squeeze()[0:3]
        # color = cmap(bin_indices[i]).squeeze()[:3]
        scaled_color = (np.array(color) * 255).astype(np.uint8)
        patch_overlay = overlay[
            x : x + downsampled_patch_size, y : y + downsampled_patch_size
        ]
        mask = tissue_mask[
            x : x + downsampled_patch_size, y : y + downsampled_patch_size
        ].astype(bool)
        patch_overlay[mask, :3] = scaled_color
        patch_overlay[mask, 3] = 128

    image = slide.read_level(level)
    image = np.array(image)[:, :, :3].transpose(1, 0, 2)
    blend = image.copy() * overlay[:, :, :3].copy()
    overlay = cv2.GaussianBlur(blend, (501, 501), 0)
    tissue_mask = cv2.medianBlur(tissue_mask, 3)
    alpha = 0.5
    array1 = np.array(image) / 255
    array2 = np.array(overlay)[:, :, :3] / 255
    output = (1 - alpha) * array1 + alpha * (array2)
    output = np.clip(output, 0, 1)
    final_image = np.zeros_like(output)
    final_image[tissue_mask == 1] = output[tissue_mask == 1]
    final_image[tissue_mask == 0] = array1[tissue_mask == 0]
    overlay = Image.fromarray((final_image * 255).astype(np.uint8).transpose(1, 0, 2))
    overlay.save(f"overlay_{slide_name}.png")

    return final_image.transpose(1, 0, 2)


def plot_overlays_side_by_side(output_1, output_2):
    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    ax[0].title.set_text("Biomarker +")
    ax[0].imshow(output_1)
    ax[0].axis("off")
    ax[1].title.set_text("Biomarker -")
    ax[1].imshow(output_2)
    ax[1].axis("off")
    plt.show()


def make_overlay_v2(
    slide_path,
    slide_name,
    patches,
    features,
    blur=True,
    tissue=False,
    output_path="",
):
    slide = SlideReader(os.path.join(slide_path, slide_name + ".svs"))
    patch_coords = [
        (int(p.split("_")[0][1:]) // downsample, int(p.split("_")[1][1:]) // downsample)
        for p in patches
        if p.endswith(".jpeg")
    ]

    cmap = cm["coolwarm"]
    overlay_shape = slide.level_dimensions[level][::-1] + (4,)
    overlay = np.ones(overlay_shape, dtype=np.uint8)

    downsampled_patch_size = patch_size // downsample

    for i, (x, y) in enumerate(patch_coords):
        if i >= len(features):
            break
        color = cmap(features[i]).squeeze()[0:3]
        scaled_color = (np.array(color) * 255).astype(np.uint8)
        patch_overlay = overlay[
            x : x + downsampled_patch_size, y : y + downsampled_patch_size
        ]
        patch_overlay[:, :, :3] = scaled_color
        patch_overlay[:, :, 3] = 128
    if blur:
        overlay_blurred = cv2.GaussianBlur(overlay, (501, 501), 0)
    else:
        overlay_blurred = overlay
    if tissue:
        _, tissue_mask = slide.get_tissue_mask(level=level)
        tissue_mask = tissue_mask.transpose(1, 0)
        tissue_mask_resized = cv2.resize(
            tissue_mask,
            (overlay_blurred.shape[1], overlay_blurred.shape[0]),
            interpolation=cv2.INTER_NEAREST,
        )
        overlay_blurred[tissue_mask_resized == 0] = [0, 0, 0, 0]
    else:
        overlay_blurred = overlay_blurred
    image = slide.read_level(level)
    image = np.array(image)[:, :, :3].transpose(1, 0, 2)
    overlay_image = overlay_blurred[:, :, :3].transpose(1, 0, 2)
    overlay_alpha = overlay_blurred[:, :, 3].transpose(1, 0)

    image_pil = Image.fromarray(image.transpose(1, 0, 2)).convert("RGBA")
    overlay_pil = Image.fromarray(overlay_image).convert("RGBA")
    overlay_alpha_pil = Image.fromarray(overlay_alpha).convert("L")

    if overlay_pil.size != image_pil.size:
        overlay_pil = overlay_pil.resize(image_pil.size, Image.NEAREST)
        overlay_alpha_pil = overlay_alpha_pil.resize(image_pil.size, Image.NEAREST)

    overlay_pil.putalpha(overlay_alpha_pil)

    # Combine images
    blended = Image.alpha_composite(image_pil, overlay_pil)

    # Save and return the final image
    final_image_pil = blended
    output_image_path = os.path.join(
        output_path,
        f"overlay_{slide_name}_feature_{feature}_blur_{blur}_tissue{tissue}.png",
    )
    final_image_pil.save(output_image_path)

    return np.array(final_image_pil).transpose(1, 0, 2)


now = time.time()
print("Getting feature values...")
features_img_1 = get_feature_values(
    patch_feature_path, patches_dir, slide_name, feature_idx=feature
)
features_img_2 = get_feature_values(
    patch_feature_path, patches_dir_2, slide_name_2, feature_idx=feature
)
print("Normalizing feature values...")
features_1, features_2 = normalize_feature_values(features_img_1, features_img_2)
# print(features_1, features_2)
print("Creating overlays...")
output_1 = make_overlay_v2(slide_path, slide_name, os.listdir(patches_dir), features_1)
output_2 = make_overlay_v2(
    slide_path, slide_name_2, os.listdir(patches_dir_2), features_2
)
print("Plotting overlays...")
plot_overlays_side_by_side(output_1, output_2)
end = time.time()
print("Time elapsed (minutes): ", (end - now) / 60)
