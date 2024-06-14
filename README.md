# AI Pathology Classifier Predicts Docetaxel Benefit in High-Risk Localized Prostate Cancer: RTOG 0521 Validation

Feature extraction and computational phenotyping used in the paper "AI Pathology Classifier Predicts Docetaxel Benefit in High-Risk Localized Prostate Cancer: RTOG 0521 Validation"

## Abstract

### Background

The RTOG 0521 trial demonstrated that integrating adjuvant docetaxel into the standard of care (SOC) of radiotherapy (RT) and long-term androgen deprivation therapy (ADT) improved overall survival (OS) in patients with high-risk localized prostate cancer (PCa). However, similar trials did not confirm these findings, emphasizing the necessity for improved patient stratification to identify docetaxel beneficiaries. We aim to develop and validate an explainable artificial intelligence-based pathology imaging classifier (APIC) to identify patients likely to benefit from docetaxel.

### Methods

Data from 350 patients enrolled in RTOG 0521, treated with either RT plus ADT or docetaxel combined with RT and ADT, were divided into training (D1: n = 84, RT plus ADT) and validation sets (D2: n = 85, RT plus ADT; D3: n = 181, docetaxel with RT plus ADT). APIC was constructed by analyzing the spatial arrangement of lymphocytes and tumor nuclear diversity in digitized biopsies.

### Results

APIC showed a significant correlation with OS in D2 (HR = 2.21, 95% CI: 1.05-4.65, p = 0.0317) but not in D3 (HR = 1.05, 95% CI: 0.61-1.80, p = 0.867). Adding docetaxel to the SOC improved OS for APIC-high patients (HR = 0.49, 95% CI: 0.26–0.92, p = 0.023), but not for APIC-low patients (HR = 1.17, 95% CI: 0.59–2.30, p = 0.656).

### Conclusion

APIC was found to be predictive of chemotherapy benefit for high-risk localized PCa patients in RTOG 0521. This study represents a potentially significant contribution in personalized treatments for PCa.

## Computational phenotyping and feature extraction

### Data

Whole-slide images of core-needle prostate biopsies are used for this study, for the computational phemnotyping and feature extraction, all regions of interest of 1024x1024 pixels with at least 70% tissue at 40x magnification are considered. The data is not publicly available. The data is available upon reasonable request from the corresponding author and approval from NRG RTOG.

### Tumor segmentation

The tumor region was segmented using a deep learning-based method. Pre-trained weights are available in the `models` directory. The tumor segmentation pipeline is available in Python, it uses Pytorch and CUDA GPU acceleration.

### Nuclei segmentation and classification

Nuclei segmentation and classification was performed using the state-of-the-art deep learning-based method HoverNet, a modified HoverNet is available in the `phenotyping/nuclei_segmentation` folder. The nuclei were classified into 6 classes: epithelium, inflammatory, connective, neoplastic, necrosis, and unknown. Lymphocytes are then identified as the nuclei classified as inflammatory.

### Spatial arrangement of lymphocytes

350 features were extracted from the spatial arrangement of lymphocyte and non-lymphocyte clusters. The feature extraction pipline is available in python in the `feature_extraction/spaTIL` folder.

### Tumor nuclear diversity

3264 features were extracted from nuclei in the tumor region (nuclei masks filtered by the tumor segmentation mask) to characterize the shape diversity as a proxy of molecular heterogeneity. The feature extraction pipeline is available in MATLAB in the `feature_extraction/nuclear_diversity` folder.

### Model construction and survival analysis

APIC was constructed using an elastic net penalized cox regression model. The model was trained on D1 and validated on D2 and D3 in each individual feature set first and then the top 7 features from each set were combined to construct APIC. The survival analysis was performed using Kaplan-Meier analysis and log-rank tests. The individual survival analysis and the combined survival analysis notebooks are available in the `survival_analysis` folder.
