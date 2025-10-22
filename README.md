# A Computational Pathology Model to Predict Docetaxel Benefit in Localized High-Risk and Metastatic Prostate Cancer

Computational phenotyping used in the paper "A Computational Pathology Model to Predict Docetaxel Benefit in Localized High-Risk and Metastatic Prostate Cancer"

## Abstract

Purpose: Docetaxel improves survival in metastatic hormone-sensitive prostate cancer (mHSPC) and high-risk localized disease, but benefits vary substantially among patients. Without predictive biomarkers, clinicians cannot identify patients who will benefit, exposing many to unnecessary toxicity. We developed and validated an artificial intelligence-based pathology image classifier (APIC) to predict docetaxel benefit.
Patients and Methods: We analyzed digitized H&E-stained biopsy specimens from two phase 3 trials: CHAARTED (286/790 patients with mHSPC) and NRG/RTOG 0521 (350/563 patients with high-risk localized disease). APIC used features capturing tumor-immune spatial interactions and nuclear heterogeneity. We evaluated the predictive value of APIC for docetaxel benefit on overall survival and castration-resistance using Cox proportional hazards with interaction terms.
Results: In CHAARTED, APIC-positive patients (56.7%) showed significant overall survival improvement with docetaxel (hazard ratio [HR] 0.52, 95% confidence interval [CI] 0.31–0.85; P = 0.008) and delayed castration-resistance (HR 0.48, 95% CI 0.33–0.71; P < 0.001), while APIC-negative patients (43.3%) showed no benefit (HR 1.31, 95% CI 0.71–2.44; P = 0.39). Treatment-APIC interactions were significant (P = 0.022 and P = 0.031). In NRG/RTOG 0521, APIC-positive patients (44.7%) demonstrated survival benefit (HR 0.49, 95% CI 0.26–0.92; P = 0.023), while APIC-negative patients (55.3%) showed no benefit. Treatment-APIC interaction was significant (P = 0.024). Predictive value remained significant after adjusting for clinical variables. Limitations include retrospective analysis and need for prospective validation.
Conclusions: APIC predicts docetaxel benefit in both metastatic and localized prostate cancer, independent of clinical factors. Validation in triplet therapy with androgen receptor pathway inhibitors is needed.


## Computational phenotyping and feature extraction

### Data

Whole-slide images of core-needle prostate biopsies were used for this study from two landmark phase III clinical trials (CHAARTED, NRG/RTOG 0521). Quality control assessment was performed using [HistoQC](https://github.com/choosehappy/HistoQC) (Janowczyk et al. 2019). The data is available upon request from the corresponding authors and approval from NRG Oncology RTOG and ECOG-ACRIN.

### Patch extraction

The whole-slide images were divided into 1024x1024 pixel patches with at least 70% tissue. A patch extraction script example is available in Python in the `phenotyping/patch_extraction` folder, patch extraction relies heavily on the [HistoPrep](https://github.com/jopo666/HistoPrep) library (Pohjonen et al. 2022).

### Tumor segmentation

The tumor region was segmented using a deep learning-based method. Pre-trained weights are available in the `models` directory. The tumor segmentation pipeline is available in the `phenotyping/tumor_segmentation` folder, it is written in Python, and uses Pytorch and CUDA GPU acceleration.

### Nuclei segmentation and classification

Nuclei segmentation and classification was performed using the state-of-the-art deep learning-based method [HoverNet](https://github.com/vqdang/hover_net) (Graham et al. 2019), a modified HoverNet is available in the `phenotyping/nuclei_segmentation` folder, there is a `run.sh` file that extracts nuclei masks from histoprep extracted patches, if you want to use this script, please modify the paths inside it. The nuclei were classified into 6 classes: epithelium, inflammatory, connective, neoplastic, necrosis, and unknown. Lymphocytes are then identified as the nuclei classified as inflammatory. Pathologist visual assessment is necessary to confirm lymphocyte classification is acceptable. The Pannuke pretrained weights are available in the following [link](https://drive.google.com/file/d/1SbSArI3KOOWHxRlxnjchO7_MbWzB4lNR/view?usp=sharing).

### Spatial arrangement of lymphocytes

350 features were extracted from the spatial arrangement of lymphocyte and non-lymphocyte clusters. The feature extraction pipeline is available in python in the `feature_extraction/spaTIL` folder.

### Tumor nuclear diversity

3264 features were extracted from nuclei in the tumor region (nuclei masks filtered by the tumor segmentation mask) to characterize the shape diversity as a proxy of molecular heterogeneity. The feature extraction pipeline is available in MATLAB in the `feature_extraction/nuclear_diversity` folder.

### Model construction and survival analysis

APIC was constructed using an elastic net penalized Cox regression model. The model was trained on development cohorts (50% of control arm patients from each trial) and validated on independent validation cohorts (remaining 50% of control patients plus all docetaxel-treated patients). The survival analysis was performed using Kaplan-Meier analysis and log-rank tests. The statistical analysis code is available in R in the `survival_analysis` folder.

## Key Findings

### CHAARTED Trial (Metastatic Disease)
- APIC-positive patients (56.7% of cohort) showed significant overall survival benefit with docetaxel (HR = 0.52, P=0.008)
- APIC-positive patients had delayed castration-resistant progression (HR = 0.48, P<0.001)
- APIC-negative patients showed no benefit from docetaxel addition
- 5-year overall survival: 40.2% with docetaxel vs 15.9% with ADT alone in APIC-positive patients

### NRG/RTOG 0521 Trial (High-Risk Localized Disease)
- APIC-positive patients (45% of cohort) demonstrated significant overall survival benefit (HR = 0.49, P=0.023)
- APIC-negative patients showed no significant survival difference with docetaxel
- 10-year overall survival: 74.4% with docetaxel vs 52.7% with standard care in APIC-positive patients

### Clinical Impact
- APIC remained independently predictive in multivariable analyses for both trials
- Treatment-APIC interactions were significant in both cohorts
- Results suggest that 44-55% of patients may not benefit from docetaxel addition
- First validated AI predictive classifier for docetaxel benefit in prostate cancer across disease stages

## Data Sharing

Data are available via:
- NRG Oncology: https://www.nrgoncology.org/Resources/Ancillary-Projects-Data-Sharing-Application
- ECOG-ACRIN: Contact for data sharing protocols

## Citation

If you use this work, please cite:
```
Medina S, Tokuyama N, et al. A Computational Pathology Model to Predict Docetaxel Benefit in Localized High-Risk and Metastatic Prostate Cancer.
```

## Contact

For questions regarding this work, please contact:
- Corresponding author: anantm@emory.edu
