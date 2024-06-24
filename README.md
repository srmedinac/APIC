# AI Pathology Classifier Predicts Docetaxel Benefit in High-Risk Localized Prostate Cancer: RTOG 0521 Validation

Feature extraction and computational phenotyping used in the paper "AI Pathology Classifier Predicts Docetaxel Benefit in High-Risk Localized Prostate Cancer: RTOG 0521 Validation"

## Abstract

### Background

The randomized RTOG 0521 trial demonstrated that the incorporation of adjuvant docetaxel into the standard of care (SOC) - radiotherapy (RT) and long-term androgen deprivation therapy (ADT) - enhanced overall survival (OS) in patients with high-risk localized prostate cancer (PCa). Despite this improvement, the trial did not meet its primary target hazard ratio (HR) of 0.49. This suggests that biomarker-informed stratification could have identified the risk groups more likely to benefit from adjuvant docetaxel. We aim validate a predictive artificial intelligence-based pathology imaging classifier (APIC) to identify patients likely to benefit from docetaxel.

### Methods

Whole-slide images of prostate biopsies from 350 patients enrolled in RTOG 0521 were divided into training (D1: n = 84, RT plus ADT) and validation cohorts (D2: n = 85, RT plus ADT; D3: n = 181, docetaxel with RT plus ADT). APIC was constructed by identifying and categorizing individual cells as either immune or cancer cells. Subsequently, a series of features were extracted, focusing on the architectural arrangement and diversity of the cancer and immune cells. APIC was validated for its ability to predict OS endpoints and the benefit of added docetaxel in D2 and D3.

### Results

 A significant correlation was observed between APIC and OS in D2 (HR = 2.21, 95% CI: 1.05-4.65, p = 0.0317), but not in D3 (HR = 1.05, 95% CI: 0.61-1.80, p = 0.867). The addition of docetaxel to SOC improved OS for patients with high APIC values (HR = 0.49, 95% CI: 0.26–0.92, p = 0.023), but not for those with low APIC values (HR = 1.17, 95% CI: 0.59–2.30, p = 0.656). This suggests that patients with high APIC values may experience increased OS when treated with docetaxel, indicating a potential role for APIC in treatment decision-making for localized high-risk PCa.

### Conclusion

APIC predicts benefit of docetaxel in localized high-risk PCa. To the best of our knowledge, this is the first validated AI predictive classifier for the benefit of docetaxel in high-risk, localized PCa.

## Computational phenotyping and feature extraction

### Data

Whole-slide images of core-needle prostate biopsies are used for this study, for the computational phemnotyping and feature extraction, all regions of interest of 1024x1024 pixels with at least 70% tissue at 40x magnification are considered. The data is not publicly available. The data is available upon reasonable request from the corresponding author and approval from NRG RTOG.

### Patch extraction

The whole-slide images were divided into 1024x1024 pixel patches with at least 70% tissue at 40x magnification. The patch extraction pipeline is available in Python in the `phenotyping/patch_extraction` folder, patch extraction relies heavily in the [HistoPrep](https://github.com/jopo666/HistoPrep) library.

### Tumor segmentation

The tumor region was segmented using a deep learning-based method. Pre-trained weights are available in the `models` directory. The tumor segmentation pipeline is available in Python, it uses Pytorch and CUDA GPU acceleration.

### Nuclei segmentation and classification

Nuclei segmentation and classification was performed using the state-of-the-art deep learning-based method HoverNet (citation/link), a modified HoverNet is available in the `phenotyping/nuclei_segmentation` folder, there is a `run.sh` file that extracts nuclei masks from histoprep extracted patches, if you want to use this script, please modify the paths inside it. The nuclei were classified into 6 classes: epithelium, inflammatory, connective, neoplastic, necrosis, and unknown. Lymphocytes are then identified as the nuclei classified as inflammatory. The Pannuke pretrained weights are available in the following [link](https://drive.google.com/file/d/1SbSArI3KOOWHxRlxnjchO7_MbWzB4lNR/view?usp=sharing).

### Spatial arrangement of lymphocytes

350 features were extracted from the spatial arrangement of lymphocyte and non-lymphocyte clusters. The feature extraction pipline is available in python in the `feature_extraction/spaTIL` folder.

### Tumor nuclear diversity

3264 features were extracted from nuclei in the tumor region (nuclei masks filtered by the tumor segmentation mask) to characterize the shape diversity as a proxy of molecular heterogeneity. The feature extraction pipeline is available in MATLAB in the `feature_extraction/nuclear_diversity` folder.

### Model construction and survival analysis

APIC was constructed using an elastic net penalized cox regression model. The model was trained on D1 and validated on D2 and D3 in each individual feature set first and then the top 7 features from each set were combined to construct APIC. The survival analysis was performed using Kaplan-Meier analysis and log-rank tests. The individual survival analysis and the combined survival analysis notebooks are available in the `survival_analysis` folder.
