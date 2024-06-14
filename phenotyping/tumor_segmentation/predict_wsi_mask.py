from unet import UNet
import torch
import numpy as np
import os
from torch import nn
from torch.utils.data import DataLoader, Dataset
from PIL import Image
from tqdm import tqdm
import pytorch_lightning as pl
from sklearn.model_selection import train_test_split
from albumentations import *
from torchvision.transforms import ToTensor
import wandb
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from scipy.ndimage import binary_dilation
import cv2
import matplotlib.pyplot as plt
from torchvision.transforms.functional import to_pil_image, pil_to_tensor
from monai.losses import DiceLoss, GeneralizedDiceFocalLoss, TverskyLoss
import torchmetrics
from torchmetrics.functional.classification import dice
import openslide
import re
from histoprep import SlideReader
from skimage.morphology import (
    binary_dilation,
    disk,
    remove_small_holes,
    remove_small_objects,
    binary_erosion,
    label,
)

Image.MAX_IMAGE_PIXELS = None
wsi = ""  # name of WSI to predict
patch_locations = ""  # where to save patches for inference
path_to_wsis = ""  # location of WSI


downsample = 1
n_classes = 2
in_channels = 3
padding = True
depth = 6
wf = 2
up_mode = "upconv"
batch_norm = True
lr = 1e-3
patch_size = 1024
edge_weight = torch.tensor(0.1)
criterion = nn.CrossEntropyLoss(reduction="none", weight=torch.tensor([0.5, 0.5]))

model = UNet(
    n_classes=n_classes,
    in_channels=in_channels,
    padding=padding,
    depth=depth,
    wf=wf,
    up_mode=up_mode,
    batch_norm=batch_norm,
)

criterion = nn.CrossEntropyLoss(reduction="none")
criterion2 = GeneralizedDiceFocalLoss(softmax=True, to_onehot_y=False)
criterion3 = TverskyLoss(alpha=0.9, beta=0.1, softmax=True, to_onehot_y=False)


class TrainingModule(pl.LightningModule):
    def __init__(
        self, model, criterion, criterion2, criterion3, train_loader, val_loader
    ):
        super().__init__()
        self.model = model
        self.criterion = criterion
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.dice = torchmetrics.Dice(num_classes=n_classes)
        self.val_dice = torchmetrics.Dice(num_classes=n_classes)
        self.criterion2 = criterion2
        self.criterion3 = criterion3

    def train_dataloader(self):
        return self.train_loader

    def val_dataloader(self):
        return self.val_loader

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y, w = batch
        y_hat = self.forward(x)
        loss_matrix = self.criterion(y_hat, y)
        pred = torch.argmax(y_hat, dim=1)

        y_onehot = torch.nn.functional.one_hot(
            y,
            num_classes=2,
        ).permute(0, 3, 1, 2)
        ce_loss = (loss_matrix * (edge_weight**w)).mean()
        tversky_loss = self.criterion3(y_hat, y_onehot)
        dice_loss = self.criterion2(y_hat, y_onehot)
        loss = ce_loss + dice_loss + tversky_loss
        self.dice(pred, y)
        self.log("train_dice_loss", dice_loss, sync_dist=True, prog_bar=True)
        self.log("train_ce_loss", ce_loss, sync_dist=True, prog_bar=True)
        self.log("train_tversky_loss", tversky_loss, sync_dist=True, prog_bar=True)
        self.log("train_loss", loss, sync_dist=True, prog_bar=True)
        self.log(
            "train_dice",
            self.dice,
            on_step=True,
            on_epoch=True,
            sync_dist=True,
            prog_bar=True,
        )
        return loss

    def validation_step(self, batch, batch_idx):
        x, y, w = batch
        y_hat = self.forward(x)
        loss_matrix = self.criterion(y_hat, y)
        pred = torch.argmax(y_hat, dim=1)

        y_onehot = torch.nn.functional.one_hot(
            y,
            num_classes=2,
        ).permute(0, 3, 1, 2)
        ce_loss = (loss_matrix * (edge_weight**w)).mean()
        tversky_loss = self.criterion3(y_hat, y_onehot)
        dice_loss = self.criterion2(y_hat, y_onehot)
        loss = ce_loss + dice_loss + tversky_loss
        self.val_dice(pred, y)
        # loss = dice_loss

        self.log("val_loss", loss, sync_dist=True, prog_bar=True)
        self.log(
            "val_dice", self.val_dice, on_epoch=True, sync_dist=True, prog_bar=True
        )
        self.log("val_tversky_loss", tversky_loss, sync_dist=True, prog_bar=True)
        self.log("val_ce_loss", ce_loss, sync_dist=True, prog_bar=True)
        self.log("val_dice_loss", dice_loss, sync_dist=True, prog_bar=True)
        random_batch = torch.randint(0, len(batch), (1,)).item()
        if batch_idx == random_batch:
            random_idx = torch.randint(0, x.size(0), (1,)).item()
            self.example_image = x[random_idx].cpu()
            self.example_ground_truth = y[random_idx].cpu().type(torch.uint8) * 255
            self.example_prediction = pred[random_idx].cpu().type(torch.uint8) * 255
        return loss

    def on_validation_epoch_end(self):
        if hasattr(self, "example_image"):
            example_input_img = to_pil_image(self.example_image)
            example_ground_truth_img = to_pil_image(self.example_ground_truth)
            example_prediction_img = to_pil_image(self.example_prediction)

            # Log to wandb
            wandb.log(
                {
                    "example_input": wandb.Image(example_input_img),
                    "example_ground_truth": wandb.Image(example_ground_truth_img),
                    "example_prediction": wandb.Image(example_prediction_img),
                }
            )

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
            optimizer, T_0=10, T_mult=1, verbose=True
        )
        return {
            "optimizer": optimizer,
            "lr_scheduler": scheduler,
            "monitor": "val_dice",
        }


training_module = TrainingModule(model, criterion, criterion2, criterion3, None, None)
training_module.load_state_dict(
    torch.load(
        "tumor_segmentation_weights.ckpt",  # path to model checkpoint
    )["state_dict"]
)

tta = Compose(
    [
        VerticalFlip(p=1),
        HorizontalFlip(p=1),
        Rotate(limit=(-90, -90), p=1),
        Transpose(p=1),
        ElasticTransform(p=1),
        HueSaturationValue(
            hue_shift_limit=(-25, 0), sat_shift_limit=0, val_shift_limit=0, p=1
        ),
    ],
)

tta_deaugmentation = Compose(
    [
        Transpose(p=1),
        Rotate(limit=(90, 90), p=1),
        HorizontalFlip(p=1),
        VerticalFlip(p=1),
    ]
)


def parse_filename(filename):
    """
    Parse the filename to extract x, y, width, and height.

    :param filename: Filename in the format 'x{X}_y{Y}_w{W}_h{H}.png'.
    :return: A tuple (x, y, width, height).
    """
    match = re.search(r"x(\d+)_y(\d+)_w(\d+)_h(\d+)", filename)
    if match:
        x, y, w, h = map(int, match.groups())
        return x, y, w, h
    else:
        raise ValueError(f"Filename format is incorrect: {filename}")


def combine_predicted_patches(wsi, downsample, predicted_patches, wsi_patch_files):
    """
    Combine predicted mask patches into the original mask with overlap handled by logical OR operation.

    :param wsi: Whole Slide Image identifier.
    :param downsample: Downsample factor for the patches.
    :param mask_patch_location: Location of the mask patches.
    :param predicted_patches: List of predicted mask patches.
    :return: Combined mask.
    """
    slide = openslide.OpenSlide(os.path.join(path_to_wsis, wsi + ".svs"))
    original_shape = (slide.level_dimensions[0][0], slide.level_dimensions[0][1])
    combined_mask = np.zeros(original_shape, dtype=np.uint8)

    for idx, filename in enumerate(tqdm(wsi_patch_files)):
        x, y, w, h = parse_filename(os.path.basename(filename))
        patch = predicted_patches[idx].T
        # Rescale coordinates and dimensions
        rescaled_x = int(x / downsample)
        rescaled_y = int(y / downsample)
        rescaled_w = int(w / downsample)
        rescaled_h = int(h / downsample)
        try:
            combined_mask[
                rescaled_x : rescaled_x + rescaled_h,
                rescaled_y : rescaled_y + rescaled_w,
            ] = patch
        except:
            print("Error with patch: ", filename)
            print("Patch shape: ", patch.shape)
            print(
                "Rescaled x: ",
                rescaled_x,
                "Rescaled y: ",
                rescaled_y,
                "rescaled w: ",
                rescaled_w,
                "rescaled h: ",
                rescaled_h,
            )

    return combined_mask.T


def patch_generator(
    slide, level, downsample, patch_size, overlap, background, save_path
):
    print("Generating patches...")
    reader = SlideReader(os.path.join(path_to_wsis, slide + ".svs"))
    downsample = reader.level_downsamples[level][1]
    downsampled_patch_size = int(patch_size * downsample)
    if not os.path.exists(os.path.join(save_path, slide, "tiles")):
        threshold, tissue_mask = reader.get_tissue_mask(level=1)
        tile_coordinates = reader.get_tile_coordinates(
            tissue_mask,
            width=downsampled_patch_size,
            overlap=overlap * downsample,
            max_background=background,
        )

        tile_metadata = reader.save_regions(
            save_path,
            tile_coordinates,
            level=level,
            threshold=threshold,
            save_metrics=False,
            overwrite=True,
        )
    else:
        print("Patches already generated.")


def predict_segmentation():
    wsi_patch_files = [
        os.path.join(patch_locations, wsi, "tiles", filename)
        for filename in os.listdir(os.path.join(patch_locations, wsi, "tiles"))
        if filename.endswith(".jpeg")
    ]
    batch = []
    print("Predicting...")
    for filename in tqdm(wsi_patch_files):
        image = np.array(Image.open(filename))
        new_image = np.transpose(image, (2, 0, 1))
        new_image = new_image / 255.0
        batch.append(new_image)
    batch = np.array(batch)
    batch = torch.from_numpy(batch).float()
    batch_preds = []
    with torch.no_grad():
        y_hat = model(batch)
        print(y_hat.shape)
        pred = torch.argmax(y_hat, dim=1)
        print(pred.shape)
        batch_preds.append(pred)
    batch_preds = torch.cat(batch_preds)
    return batch_preds, wsi_patch_files


def predict_segmentation_tta():
    wsi_patch_files = [
        os.path.join(patch_locations, wsi, "tiles", filename)
        for filename in os.listdir(os.path.join(patch_locations, wsi, "tiles"))
        if filename.endswith(".jpeg")
    ]
    times_to_augment = 5
    print("Predicting with TTA...")
    batch_preds = []
    all_patches = []
    for filename in tqdm(wsi_patch_files):
        # Read and process the image

        image = np.array(Image.open(filename))

        # Apply TTA and predict for each augmented image
        all_preds = []
        for _ in range(times_to_augment):
            tta_image = tta(image=image)["image"]
            tta_image = (
                torch.from_numpy(np.transpose(tta_image, (2, 0, 1)) / 255.0)
                .unsqueeze(0)
                .float()
            )

            with torch.no_grad():
                y_hat = model(tta_image)
                pred = torch.argmax(y_hat, dim=1)
                all_preds.append(pred.squeeze(0))

        # De-augment predictions and aggregate
        deaug_preds = [
            tta_deaugmentation(image=image, mask=pred.numpy())["mask"]
            for pred in all_preds
        ]
        aggregated_pred = np.mean(np.stack(deaug_preds), axis=0)
        aggregated_pred = aggregated_pred > 0.5

        batch_preds.append(aggregated_pred)
        all_patches.append(image)

    # Combine all predictions
    all_patches = np.stack(all_patches)
    batch_preds = np.stack(batch_preds)

    return batch_preds, wsi_patch_files, all_patches


patch_generator(
    wsi,
    level=0,
    downsample=downsample,
    patch_size=patch_size,
    overlap=0,
    background=0.7,
    save_path=patch_locations,
)
batch_preds, wsi_patch_files = predict_segmentation()

generated_mask = torch.tensor(
    combine_predicted_patches(wsi, downsample, batch_preds, wsi_patch_files)
).numpy()


cleaned_mask = remove_small_objects(generated_mask.astype(bool), 5000)
filled_mask = remove_small_holes(cleaned_mask.astype(bool), 200)
selem = disk(5)
dilated_mask = binary_dilation(filled_mask, selem)
opened_mask = binary_erosion(dilated_mask, selem)
final_mask = remove_small_objects(opened_mask.astype(bool), 5000)

generated_mask = to_pil_image(generated_mask * 255)
generated_mask.save(
    os.path.join(patch_locations, wsi + f"_mask_{wsi}_1024patch3losses.png")
)


# def overlay_mask(slide, mask, level, output_path, color=(128, 255, 128, 128)):
#     # Load the image
#     reader = openslide.OpenSlide(os.path.join(path_to_wsis, slide))
#     image = reader.read_region((0, 0), level, reader.level_dimensions[level])

#     rgba_mask = np.zeros((mask.shape[0], mask.shape[1], 4), dtype=np.uint8)
#     rgba_mask[mask == 1] = color
#     rgba_mask[mask == 0] = (0, 0, 0, 0)  # Set transparent for the 0s
#     rgba_mask = Image.fromarray(rgba_mask, "RGBA")

#     # Resize mask if it's not the same size as the image
#     if rgba_mask.size != image.size:
#         rgba_mask = rgba_mask.resize(image.size, Image.NEAREST)

#     # Overlay the mask on the image
#     overlaid_image = Image.alpha_composite(image.convert("RGBA"), rgba_mask)

#     # Save or return the overlayed image
#     overlaid_image.save(output_path)
#     return overlaid_image


# path_to_wsis = ""  # location of WSI
# seg_mask = np.array(
#     Image.open(os.path.join(patch_locations, wsi + f"_mask_{wsi}.png"))
# )
# seg_mask = np.where(seg_mask >= 128, 1, 0)
# overlay = overlay_mask(
#     wsi + ".svs",
#     seg_mask,
#     1,
#     f"/Users/srmedinac/Desktop/slide_with_pred_mask_{wsi}.png",
# )
