# APIC - AI-based Pathology Image Classifier

Reference research implementation of APIC (Artificial intelligence-based Pathology Image Classifier) from:

> Medina S, Tokuyama N, Hammouda K, Pathak T, Mirtti T, Fu P, Gupta S, Lal P, Sandler HM, Correa R, Chafe S, Shah A, Efstathiou JA, Hoffman K, Straza M, Hallman MA, Jordan R, Pugh SL, Sweeney CJ, Madabhushi A. **A Computational Pathology Model to Predict Docetaxel Benefit in Localized High-Risk and Metastatic Prostate Cancer.** *Clin Cancer Res* 2025. [DOI: 10.1158/1078-0432.CCR-25-3327](https://doi.org/10.1158/1078-0432.CCR-25-3327)

and extended to enzalutamide benefit prediction in:

> Medina S, Tokuyama N, et al. **An Artificial Intelligence Pathology Image Classifier Predicts Enzalutamide Benefit in Metastatic Hormone-Sensitive Prostate Cancer: A Biomarker Analysis of the Phase III ENZAMET Trial (ANZUP 1304).** *Under review.*

## Overview

APIC analyzes **Whole Slide Images (WSI)** from prostate cancer biopsies to generate a **patient-level risk score** that predicts whether a patient with **advanced prostate cancer** will benefit from intensification with **docetaxel** or **enzalutamide**. The pipeline performs tumor segmentation, nuclei segmentation and classification, and spatial feature extraction from tumor tissue, then assigns each patient as APIC-Positive or APIC-Negative.

### Clinical context

For patients with localized high-risk or metastatic hormone-sensitive prostate cancer, adding docetaxel or an androgen receptor pathway inhibitor (ARPI) to androgen deprivation therapy (ADT) can improve outcomes, but not all patients benefit equally. APIC provides a computational biomarker from routine H&E-stained biopsy images to help identify patients most likely to benefit from treatment intensification.

### Summary of results

**Docetaxel benefit — CHAARTED and NRG/RTOG 0521 (Medina et al., Clin Cancer Res 2025).** In CHAARTED (mHSPC, n=286), APIC-positive participants (56.7%) derived significant OS benefit from docetaxel (HR 0.52, 95% CI 0.31–0.85, p=0.008); APIC-negative participants did not. In NRG/RTOG 0521 (high-risk localized, n=350), APIC-positive participants (44.7%) derived significant OS benefit from docetaxel (HR 0.49, 95% CI 0.26–0.92, p=0.023); APIC-negative participants did not. Treatment-by-APIC interactions were significant in both trials.

**Enzalutamide benefit — ENZAMET (Medina, Tokuyama et al., under review).** In the ENZAMET biomarker cohort (mHSPC, n=492), the prespecified treatment-by-APIC interaction for overall survival was significant (p=0.014). APIC-*negative* participants derived significant OS benefit from adding enzalutamide to ADT versus conventional NSAA (e.g., HR 0.54, 95% CI 0.35–0.85, p=0.008 in the biopsy subgroup), while APIC-*positive* participants did not (HR 1.06, 95% CI 0.54–2.08, p=0.865). The interaction held in a prespecified sensitivity analysis excluding participants who received early concurrent docetaxel and after adjustment for disease volume, timing of metastatic presentation, ECOG, age, Gleason score, baseline PSA, and specimen type. The APIC model and cutoff were locked from CHAARTED and applied to ENZAMET without modification.

Note the directionality across trials: APIC-*positive* tumors benefit from docetaxel, while APIC-*negative* tumors benefit from enzalutamide. Together the two studies support APIC as a candidate predictive biomarker that may inform selection between chemotherapy and ARPI intensification in advanced prostate cancer.

### Two ways to run APIC

| Goal | Use | Repository |
|---|---|---|
| Score new slides with a single command | Dockerized pipeline | [madabhushilabapic/apic](https://hub.docker.com/r/madabhushilabapic/apic) |
| Inspect, modify, or reproduce the research pipeline from raw scripts | **This repository** | [srmedinac/APIC](https://github.com/srmedinac/APIC) |

If you only want APIC risk scores and PDF reports from WSIs, use the Docker container below. The rest of this README documents the research-code pipeline in the current repository.

---

## Quick start with the Docker container

```bash
docker pull madabhushilabapic/apic:latest

docker run --gpus all \
  -v /path/to/your/slides:/data/input_slides:ro \
  -v /path/to/your/output:/data/output \
  madabhushilabapic/apic:latest \
  -i "/data/input_slides/slide_filename.svs" \
  -o /data/output
```

GPU is optional; the container falls back to CPU. The container supports single-slide, batch, multi-slide-per-patient, and research (auto-grouped) modes. See the [Docker repository README](https://hub.docker.com/r/madabhushilabapic/apic) for full run modes and flags.

---

## Research pipeline (this repository)

### Repository contents

```
APIC/
├── phenotyping/
│   ├── patch_extraction/         # Tile a WSI into 1024x1024 patches (HistoPrep)
│   ├── tumor_segmentation/       # U-Net tumor mask from patches (PyTorch)
│   └── nuclei_segmentation/      # HoverNet nuclei segmentation + classification
├── feature_extraction/
│   ├── spaTIL/                   # 350 spatial TIL features (Python)
│   └── nuclear_diversity/        # 3264 nuclear-diversity features (MATLAB)
├── extra_scripts/                # Patch -> slide -> patient feature aggregation
├── survival_analysis/            # Feature selection + APIC Cox model (R, Python)
├── models/                       # Pretrained tumor segmentation weights (git-lfs)
└── LICENSE.md
```

### System requirements

**Operating system**
- Linux (tested on Ubuntu 20.04 and 22.04). Recommended for the full pipeline.
- macOS 13+ is supported for all steps except HoverNet inference, which requires CUDA.
- Windows is supported for the MATLAB nuclear-diversity module.

**Hardware**

| Resource | Minimum | Recommended |
|---|---|---|
| CPU | 8 cores | 16+ cores |
| System RAM | 32 GB | 64 GB |
| GPU | 1x NVIDIA GPU with 12 GB VRAM (RTX 3060, T4) | 1x NVIDIA GPU with 24 GB VRAM (RTX 3090/4090, A5000, A100) |
| Disk | 1 TB (tiles and masks account for most of the footprint) | 2 TB NVMe SSD |

Peak RAM during spaTIL feature extraction scales with `num_processes` in `feature_extraction/spaTIL/spaTIL_main.py` (roughly 2 GB per worker at 2048x2048 patches). Reduce `num_processes` if out-of-memory errors occur.

**Software**
- Python 3.10 (tested with 3.9–3.11)
- PyTorch 2.0+ with CUDA 11.7 or 11.8
- R 4.2+
- MATLAB R2020a+ (Image Processing and Statistics toolboxes; needed only for the nuclear-diversity features)
- OpenSlide 3.4+ (system library required by `openslide-python`)
- Git LFS (required to pull `models/tumor_segmentation_weights.ckpt`)

### Installation

```bash
# 1. Install system dependencies (Ubuntu example)
sudo apt-get install openslide-tools git-lfs
git lfs install

# 2. Clone with LFS to pull model weights
git clone https://github.com/srmedinac/APIC.git
cd APIC
git lfs pull

# 3. Create a conda environment and install Python dependencies
conda create -n apic python=3.10 -y
conda activate apic
pip install torch==2.1.0 torchvision --index-url https://download.pytorch.org/whl/cu118
pip install histoprep openslide-python pytorch-lightning==2.1.0 monai \
            albumentations opencv-python scikit-image==0.22.0 \
            scikit-learn pandas numpy scipy tqdm networkx shapely Pillow \
            geojson wandb matplotlib

# 4. HoverNet (separate environment; follow upstream instructions)
# https://github.com/vqdang/hover_net
# PanNuke weights: https://drive.google.com/file/d/1SbSArI3KOOWHxRlxnjchO7_MbWzB4lNR/view

# 5. R packages for survival analysis
Rscript -e 'install.packages(c("survival","survminer","dplyr","caret","ggplot2","gridExtra"))'
```

---

## Input data requirements

### Whole-slide image formats

The pipeline reads WSIs via `histoprep` and `openslide`.

| Format | Extension | Tested | Notes |
|---|---|---|---|
| Aperio | `.svs` | Yes | Primary format used in CHAARTED and RTOG 0521 |
| Hamamatsu | `.ndpi` | Yes | Supported by OpenSlide |
| Leica | `.scn` | Yes | Supported by OpenSlide |
| 3DHISTECH | `.mrxs` | Partial | OpenSlide-supported; associated directory must be present |
| Generic tiled TIFF | `.tif`, `.tiff` | Yes | Must be pyramidal (BigTIFF or tiled TIFF) |
| OME-TIFF | `.ome.tif`, `.ome.tiff` | Yes | Must be pyramidal; single-resolution OME-TIFF is not supported |
| Philips iSyntax | `.isyntax` | No | Not supported by OpenSlide; convert to pyramidal TIFF first |

### Resolution and compression

- **Required base magnification:** 40x (0.25 µm/pixel). Slides digitized at 20x can be processed but patch size must be halved (512x512) and the `r` parameter in `spaTIL_main.py` doubled.
- **Compression:** JPEG (quality ≥ 80) and JPEG2000 are supported. LZW-compressed TIFF is also supported. Uncompressed pyramidal TIFF works but is disk-intensive.
- **Tile structure:** The pipeline assumes the underlying WSI contains a tiled pyramid; single-strip TIFFs will fail during random access.

### Quality control

Input slides should be passed through [HistoQC](https://github.com/choosehappy/HistoQC) (Janowczyk et al., 2019) to exclude pen marks, blur, and scanning artefacts. A HistoQC tissue mask may be supplied to `spaTIL_main.py` via `config["histoqc_mask"]`; otherwise a full-tile mask is assumed.

### Clinical file schema

The survival analysis expects a single CSV per trial with at least the following columns:

| Column | Type | Description |
|---|---|---|
| `group_uid` | string | Patient identifier |
| `os` | numeric | Overall survival time (months for CHAARTED, years for RTOG; converted internally) |
| `dead` | 0/1 | Death indicator |
| `TT_CRPC` | numeric | Time to castration-resistant progression (months) |
| `CRPC` | 0/1 | CRPC event indicator |
| `ASSIGNED_TX_ARM` / `RX` | string / integer | Treatment arm |
| `has_features` | boolean | Whether a valid feature vector exists for the patient |
| Feature columns | numeric | Six selected features per trial (listed in `survival_analysis/APIC.R`) |

A separate `patient_slides.csv` maps `USI` (patient) to `Image ID` (slide) for multi-slide aggregation in `extra_scripts/get_patient_features.py`.

---

## End-to-end pipeline

The pipeline converts a directory of WSIs into an APIC risk score. Per-step scripts expose input/output paths as module-level variables; set these before running.

### Step 1. Patch extraction
Input: WSI file. Output: folder of 1024x1024 JPEG tiles with at least 30% tissue.
```bash
# Edit paths in phenotyping/patch_extraction/extract_patches.py, then:
python phenotyping/patch_extraction/extract_patches.py
```

### Step 2. Tumor segmentation
Input: tile folder from Step 1 + `models/tumor_segmentation_weights.ckpt`. Output: binary tumor mask per slide.
```bash
# Edit wsi, patch_locations, path_to_wsis at the top of predict_wsi_mask.py
python phenotyping/tumor_segmentation/predict_wsi_mask.py
```

### Step 3. Nuclei segmentation and classification
Input: tile folder. Output: per-tile PNG masks with 6-class labels (epithelium, inflammatory, connective, neoplastic, necrosis, unknown); lymphocytes are the "inflammatory" class (label value 128).
```bash
# Edit paths inside phenotyping/nuclei_segmentation/run.sh
bash phenotyping/nuclei_segmentation/run.sh
```
Requires a working HoverNet install (https://github.com/vqdang/hover_net) with PanNuke weights. Pathologist review of a sample of lymphocyte overlays is recommended before proceeding.

### Step 4a. spaTIL features (Python)
Input: patch folder + HoverNet masks. Output: 350-dimensional feature CSV per tile.
```bash
# Edit the config dict at the top of spaTIL_main.py:
#   patches_dir, nuclei_masks_dir, results_features_dir, patch_size, r, num_processes
python feature_extraction/spaTIL/spaTIL_main.py
```

### Step 4b. Nuclear-diversity features (MATLAB)
Input: tumor-masked patches. Output: 3264-dimensional feature vector per tile.
```matlab
% In MATLAB:
cd feature_extraction/nuclear_diversity
runme
```

### Step 5. Tile -> slide -> patient aggregation
Input: per-tile feature CSVs. Output: slide-level `.npy` then patient-level `.npy` and a combined patient feature CSV. Patient-level features are averaged across all slides of a patient using `nanmean`.
```bash
# Edit results_path / save_path at the top of each script, then:
python extra_scripts/get_slide_features.py
python extra_scripts/get_patient_features.py
```

### Step 6. APIC model training, validation, and risk scoring
Input: `all_patient_features.csv` joined with clinical columns. Output: APIC-Positive / APIC-Negative assignment, Kaplan-Meier plots, Cox summaries.
```bash
Rscript survival_analysis/APIC.R                       # Trains APIC on CHAARTED and RTOG
Rscript survival_analysis/multivariable_analysis.R     # Adjusts for clinical covariates
```
For the feature-selection step used during model development:
```bash
python survival_analysis/feature_selectino.py
```
Notebooks reproducing the published figures are in `survival_analysis/survival_analysis_combined.ipynb` and `survival_analysis_individual_feats.ipynb`.

---

## Quickstart with a single example slide

The following minimal example runs Steps 1–6 on one slide. Wall-clock time on a single A5000 GPU + 16-core CPU is approximately 45 minutes per slide, with nuclei segmentation taking the majority of that time.

```bash
export WSI=/data/example/slide.svs
export OUT=/data/example/APIC_out

mkdir -p $OUT/patches $OUT/tumor $OUT/nuclei $OUT/spatil $OUT/nuclear

# Step 1
python phenotyping/patch_extraction/extract_patches.py   # edit paths to $WSI, $OUT/patches

# Step 2
python phenotyping/tumor_segmentation/predict_wsi_mask.py

# Step 3
bash phenotyping/nuclei_segmentation/run.sh

# Step 4a and 4b
python feature_extraction/spaTIL/spaTIL_main.py
matlab -batch "cd feature_extraction/nuclear_diversity; runme"

# Step 5
python extra_scripts/get_slide_features.py
python extra_scripts/get_patient_features.py

# Step 6 - score against the locked APIC model
Rscript survival_analysis/APIC.R
```

### Expected output
After Step 6, the R script prints Cox-model and treatment-APIC-interaction summaries, saves Kaplan-Meier plots to the working directory, and assigns each patient an `APIC-Positive` or `APIC-Negative` label based on the locked cutoff (33rd percentile for CHAARTED, 50th percentile for RTOG).

---

## Key findings

### CHAARTED (metastatic disease)
- APIC-positive (56.7%): significant OS benefit with docetaxel (HR 0.52, P = 0.008) and delayed CRPC (HR 0.48, P < 0.001).
- APIC-negative (43.3%): no benefit from docetaxel addition.

### NRG/RTOG 0521 (high-risk localized disease)
- APIC-positive (44.7%): significant OS benefit (HR 0.49, P = 0.023).
- APIC-negative (55.3%): no significant survival difference with docetaxel.

APIC remained independently predictive after adjustment for clinical variables in both trials.

---

## Data availability

Clinical trial data used in this study, including histology images and deidentified clinical information, were obtained from NRG Oncology and ECOG-ACRIN under data use agreements. Access is governed by the policies of the originating cooperative groups. Researchers interested in accessing these data must submit an ancillary project application through the NRG Oncology data sharing portal (https://www.nrgoncology.org/Resources/Ancillary-Projects-Data-Sharing-Application) and ECOG-ACRIN (https://ecog-acrin.org/). Processed data derived from the histology images may be made available to qualified researchers upon reasonable request to the corresponding author and execution of appropriate data transfer agreements, subject to cooperative group approval.

---

## Citation

If you use this pipeline in your research, please cite the appropriate reference(s):

**Original APIC development and docetaxel benefit prediction (CHAARTED, NRG/RTOG 0521):**

> Medina S, Tokuyama N, Hammouda K, Pathak T, Mirtti T, Fu P, Gupta S, Lal P, Sandler HM, Correa R, Chafe S, Shah A, Efstathiou JA, Hoffman K, Straza M, Hallman MA, Jordan R, Pugh SL, Sweeney CJ, Madabhushi A. **A Computational Pathology Model to Predict Docetaxel Benefit in Localized High-Risk and Metastatic Prostate Cancer.** *Clin Cancer Res* 2025. [DOI: 10.1158/1078-0432.CCR-25-3327](https://doi.org/10.1158/1078-0432.CCR-25-3327)

**Enzalutamide benefit prediction (ENZAMET):**

> Medina S, Tokuyama N, et al. **An Artificial Intelligence Pathology Image Classifier Predicts Enzalutamide Benefit in Metastatic Hormone-Sensitive Prostate Cancer: A Biomarker Analysis of the Phase III ENZAMET Trial (ANZUP 1304).** *Under review.*

```bibtex
@article{medina2025apic,
  author  = {Medina, Sebastian and Tokuyama, Naoto and Hammouda, Kamal and Pathak, Tilak
             and Mirtti, Tuomas and Fu, Pingfu and Gupta, Shilpa and Lal, Priti
             and Sandler, Howard M. and Correa, Rohann and Chafe, Susan and Shah, Amit
             and Efstathiou, Jason A. and Hoffman, Karen and Straza, Michael
             and Hallman, Mark A. and Jordan, Richard and Pugh, Stephanie L.
             and Sweeney, Christopher J. and Madabhushi, Anant},
  title   = {A Computational Pathology Model to Predict Docetaxel Benefit in Localized
             High-Risk and Metastatic Prostate Cancer},
  journal = {Clinical Cancer Research},
  year    = {2025},
  doi     = {10.1158/1078-0432.CCR-25-3327}
}

@unpublished{medina2025apicenzamet,
  author = {Medina, Sebastian and Tokuyama, Naoto and others},
  title  = {An Artificial Intelligence Pathology Image Classifier Predicts Enzalutamide
            Benefit in Metastatic Hormone-Sensitive Prostate Cancer: A Biomarker Analysis
            of the Phase III ENZAMET Trial (ANZUP 1304)},
  note   = {Under review}
}
```

## Third-party components

- [HistoPrep](https://github.com/jopo666/HistoPrep) (Pohjonen et al., 2022) for patch extraction.
- [HoverNet](https://github.com/vqdang/hover_net) (Graham et al., 2019) for nuclei segmentation and classification.
- [HistoQC](https://github.com/choosehappy/HistoQC) (Janowczyk et al., 2019) for slide-level quality control.
- Cellular diversity features adapted from Lu et al., *Lancet Digit Health* 2020 and Lu et al., *Mod Pathol* 2017 (see `feature_extraction/nuclear_diversity/README.md`).

---

## Disclaimer

**For research purposes only.** This tool is not intended for clinical use and should not be used to diagnose, treat, or make clinical decisions for any patient. The predictions generated by this pipeline have not been validated for clinical practice and are meant solely for academic and research applications.

## License

This software is licensed under the [Emory University License](LICENSE.md) for **non-commercial research purposes only**. See the LICENSE file for full terms and conditions. For commercial licensing inquiries, contact Emory University.

## Contact

Corresponding author: Anant Madabhushi (anantm@emory.edu). For questions about this repository, open a GitHub issue.
