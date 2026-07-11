# 🧬 Code and Data for Protein-Fingerprint-Based Prediction of Cadmium Concentrations in Ark Clams

![R](https://img.shields.io/badge/R-4.3.3-276DC3?logo=r&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-Research%20code-orange)

This repository contains the R code and supporting datasets used to predict total wet-weight cadmium (Cd) concentrations in ark clam tissues from cytosolic HPLC-UV protein fingerprints.

The current workflow focuses on two closely related edible ark clam species, *Anadara kagoshimensis* and *Tegillarca granosa*. Multiple chromatographic feature sets are compared using Elastic Net regression and leave-one-species-out cross-validation.

SEC-ICP-MS profiles are used to identify Cd-associated retention-time regions and to guide penalty weighting during model training. **SEC-ICP-MS signals are not included in the predictor matrix. Predictions for held-out samples are generated from HPLC-UV fingerprints only.**

---

## 🔬 Study Overview

The analysis integrates three types of data:

- **HPLC-UV protein fingerprints** used as model predictors;
- **SEC-ICP-MS chromatographic profiles** used to define and interpret Cd-associated retention-time regions; and
- **tissue metal concentrations** used to construct the prediction target.

The current model predicts total Cd concentrations in:

- gill;
- viscera;
- muscle.

The response variable is expressed on a wet-weight basis and is log10-transformed before model fitting.

Total Cd concentration is calculated as:

```text
Cd_total_ug_g = Cd111_back_ug_g + accumulated_Cd113_ug_g
```

The concentration table also contains corresponding information for Ni, Cu, Zn, and Pb. However, the current modelling script is configured specifically for Cd.

---

## 🎯 Modelling Objective

The workflow evaluates whether cytosolic protein fingerprints measured by HPLC-UV can predict tissue Cd concentrations across ark clam species.

Model transferability is assessed by training the model on one species and testing it on the other species.

The central analytical principle is:

> HPLC-UV fingerprints provide the prediction features, whereas SEC-ICP-MS provides a Cd-specific metallomic prior for feature weighting and interpretation.

---

## 🧪 Analytical Workflow

The R script performs the following steps:

1. reads the HPLC-UV, SEC-ICP-MS, and tissue metal datasets;
2. standardizes sample identifiers and tissue names;
3. calculates total tissue metal concentrations;
4. matches HPLC-UV profiles with tissue metal measurements;
5. converts chromatographic signals into multiple feature sets;
6. applies closure normalization and centered log-ratio transformation;
7. derives Cd-associated retention-time regions from training-set SEC-ICP-MS profiles;
8. assigns reduced Elastic Net penalties to HPLC-UV features located within or near Cd-associated regions;
9. fits Elastic Net models with internal cross-validation;
10. evaluates interspecific transferability using leave-one-species-out cross-validation; and
11. generates model-performance summaries, selected-feature tables, and publication-quality figures.

---

## 🧩 HPLC-UV Feature Sets

The script constructs the following HPLC-UV feature sets:

| Feature set | Description |
|---|---|
| `fixed_auc_clr` | AUC values calculated within fixed 1-min retention-time bins, followed by closure normalization and CLR transformation |
| `quantile_clr` | AUC values calculated within quantile-based retention-time bins, followed by closure normalization and CLR transformation |
| `wavelet_clr` | Wavelet-denoised chromatograms summarized using fixed retention-time AUC bins and CLR transformation |
| `peak_clr` | AUC values calculated around consensus HPLC-UV peaks and transformed using CLR |
| `fixed_auc_clr_roi_weighted` | Fixed-bin HPLC-UV features fitted with SEC-ICP-MS-informed Elastic Net penalty factors |
| `wavelet_clr_roi_weighted` | Wavelet-denoised HPLC-UV features fitted with SEC-ICP-MS-informed Elastic Net penalty factors |

Negative UV values are set to zero before AUC calculation. A pseudo-count is added before closure normalization and CLR transformation.

---

## 🧠 SEC-ICP-MS-Informed Penalty Weighting

Cd-associated retention-time regions are derived from SEC-ICP-MS profiles using training samples only.

For each leave-one-species-out split:

- the mean training-set Cd profile is calculated;
- the profile is smoothed;
- regions exceeding a defined fraction of the maximum signal are detected;
- narrow regions are removed;
- nearby regions are merged; and
- HPLC-UV variables within or near these regions receive reduced Elastic Net penalty factors.

The default penalty factors are:

| Region | Penalty factor |
|---|---:|
| Core Cd-associated region | 0.25 |
| Flanking region | 0.50 |
| Background region | 1.00 |

Lower penalty factors allow variables within Cd-associated chromatographic regions to enter the model more readily.

The measured `Cd113` SEC-ICP-MS signal is used for ROI definition and graphical interpretation. It is not used as a predictor or response variable.

---

## ✅ Model Validation

Model transferability is evaluated using leave-one-species-out cross-validation rather than a random train-test split.

The two ark clam species are alternately used as the training and testing groups:

1. all samples of *Tegillarca granosa* are used to train the model, and the fitted model is used to predict total Cd concentrations in *Anadara kagoshimensis*;
2. all samples of *Anadara kagoshimensis* are then used to train the model, and the fitted model is used to predict total Cd concentrations in *Tegillarca granosa*.

For each validation split, the training data contain:

- HPLC-UV protein fingerprints as predictors;
- measured total Cd concentrations as the response; and
- SEC-ICP-MS profiles used only to derive Cd-associated retention-time regions and penalty factors.

For the held-out species, only the HPLC-UV fingerprints are supplied to the fitted model to generate predicted total Cd concentrations.

The measured Cd concentrations of the held-out species are not used during model training. They are retained only for comparison with the predicted values and for calculating RMSE, MAE, and R².

SEC-ICP-MS profiles from the held-out species are also excluded when defining Cd-associated retention-time regions. The ROI information is derived from the training species only, thereby preventing information leakage from the test set.

The validation procedure can therefore be summarized as:

```text
Train on T. granosa
    ↓
Predict Cd concentrations in A. kagoshimensis
    ↓
Compare predictions with measured Cd concentrations

Train on A. kagoshimensis
    ↓
Predict Cd concentrations in T. granosa
    ↓
Compare predictions with measured Cd concentrations

The model response is:

```text
log10(total Cd concentration + pseudo-count)
```

Predictions are subsequently back-transformed to the original wet-weight concentration scale for visualization.

---

## 📁 Repository Contents

```text
Ark-clam-metal-fingerprint-model/
├── README.md
├── LICENSE
├── ark_clam_cd_fingerprint_model.R
├── ark_clam_hplc_uv_profiles.xlsx
├── ark_clam_sec_icp_ms_profiles.xlsx
└── ark_clam_subcellular_metal_concentrations.xlsx
```

### File descriptions

- **`ark_clam_cd_fingerprint_model.R`**  
  Main R script for data preprocessing, HPLC-UV feature engineering, SEC-ICP-MS-informed ROI definition, Elastic Net modelling, leave-one-species-out validation, and figure generation.

- **`ark_clam_hplc_uv_profiles.xlsx`**  
  HPLC-UV chromatographic profiles used to construct the model predictors.

- **`ark_clam_sec_icp_ms_profiles.xlsx`**  
  SEC-ICP-MS chromatographic profiles used to derive Cd-associated retention-time regions and interpret selected HPLC-UV features.

- **`ark_clam_subcellular_metal_concentrations.xlsx`**  
  Tissue metal concentration data used to construct the model response.

- **`LICENSE`**  
  MIT License governing reuse of the code in this repository.

---

## 📋 Expected Input Structure

### `ark_clam_hplc_uv_profiles.xlsx`

The first column contains retention time in minutes. The remaining columns contain sample-specific HPLC-UV signals.

The script renames the first column to:

```text
time_min
```

### `ark_clam_sec_icp_ms_profiles.xlsx`

The script expects the following metadata columns:

```text
sample
Species
Tissue
rep
Time_sec
```

The current Cd workflow also requires:

```text
Cd113
```

Additional isotope-resolved metal columns may include:

```text
Cd111
Cu65
Zn68
Ni61
Pb206
```

### `ark_clam_subcellular_metal_concentrations.xlsx`

The script expects sample metadata and wet-weight metal concentration columns, including:

```text
sample
Species
Tissue
Fraction
rep
ww_g
Cd111_back_ug_g
accumulated_Cd113_ug_g
```

Only rows with:

```text
Fraction == "total"
```

and tissues classified as gill, viscera, or muscle are used for the current model.

---

## 💻 Software Requirements

The analysis was developed and tested using **R 4.3.3**.

Required R packages:

```r
install.packages(
  c(
    "readxl",
    "tidyverse",
    "glmnet",
    "pracma",
    "waveslim",
    "cowplot"
  )
)
```

The packages are used as follows:

| Package | Main purpose |
|---|---|
| `readxl` | Reading Excel input files |
| `tidyverse` | Data cleaning, transformation, plotting, and file output |
| `glmnet` | Elastic Net regression and internal cross-validation |
| `pracma` | Consensus peak detection |
| `waveslim` | Wavelet denoising |
| `cowplot` | Figure alignment and multi-panel assembly |

---

## ▶️ Running the Analysis

1. Download or clone this repository.
2. Place the R script and all three Excel files in the repository root directory.
3. Set the repository root as the working directory.
4. Run:

```r
source("Cd_SPEicpms_model_2.R")
```

Alternatively, open the R script in RStudio and run it from the repository root directory.

The input filenames must match the following names exactly:

```text
ark_clam_hplc_uv_profiles.xlsx
ark_clam_sec_icp_ms_profiles.xlsx
ark_clam_subcellular_metal_concentrations.xlsx
```

---

## 📊 Main Outputs

The script produces or stores:

- matched HPLC-UV and tissue-metal sample summaries;
- Cd concentration summaries;
- feature matrices for all HPLC-UV feature sets;
- leave-one-species-out predictions;
- RMSE, MAE, and R² summaries;
- selected Elastic Net coefficients;
- SEC-ICP-MS-derived Cd retention-time regions;
- observed-versus-predicted Cd plots;
- coefficient and SEC-ICP-MS overlay plots; and
- a combined publication-quality figure.

Key output files include:

```text
Cd_ElasticNet_selected_features_summary.csv
Cd_mode2-1.png
Cd_mode2-1.pdf
```

Additional model objects and summaries are retained in the R environment.

---

## ⚠️ Important Notes

- SEC-ICP-MS signals are **not** included in the model predictor matrix.
- SEC-ICP-MS information is used only to define penalty factors and aid interpretation.
- ROI detection is performed separately within each training split to avoid using held-out species information.
- Predictions are based only on HPLC-UV fingerprints available for the held-out samples.
- The current script is configured for Cd, although the concentration table includes additional metals.
- Input filenames, worksheet names, and column names must remain consistent with the script.

---

## 👥 Repository Maintainers

**Xue Cao**  
Shantou University  
Email: xuecao@stu.edu.cn

**Qiao-Guo Tan**  
Xiamen University  
Email: tanqg@xmu.edu.cn

---

## 📖 Citation

The associated manuscript is currently under preparation. Full citation information and the article DOI will be added after publication.

```text
https://github.com/Xuecao1005/Ark-clam-metal-fingerprint-model
```

---

## ⚖️ License

The code in this repository is distributed under the MIT License. See the `LICENSE` file for details.

Please cite the associated study and this repository when reusing the code or data.

---

## 📬 Contact

For questions about the datasets, modelling workflow, or reuse of repository materials, please contact:

**Xue Cao**  
Shantou University  
Email: xuecao@stu.edu.cn
