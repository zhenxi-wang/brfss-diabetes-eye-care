# Retinal photography and diabetes eye-care ascertainment in BRFSS

Reproducibility materials for the 2025 and 2022 Behavioral Risk Factor Surveillance System (BRFSS) analysis. The study compares self-reported annual dilation with the combined definition of annual dilation or retinal photography. The analysis covers 20 full-module states in 2025 and 12 in 2022. See [the protocol and dated amendments](PROTOCOL.md) for eligibility and analysis decisions.

## Contents

| Directory | Contents |
| --- | --- |
| [code](code) | CDC data acquisition, variable extraction, paired survey analysis, multinomial models, nonresponse analyses, tables, figures, and verification |
| [data/derived](data/derived) | Selected variables and analysis files generated from the CDC public-use records, with the survey design variables needed for domain variance |
| [results](results) | Machine-readable aggregate estimates, sensitivity analyses, Tables 1–2 and Supplementary Tables S1–S20, and [coefficients and covariance matrices](results/model_parameters) |
| [figures](figures) | Three main and three supplementary figures as PDF, PNG, and SVG |
| [sources](sources) | Official data URLs and archived checksums, extraction manifests, questionnaires, codebooks, and state-module eligibility |
| [qa](qa) | Statistical verification and artifact checks for the 17 November 2026 analysis release |

## Reproduce the estimates and figures

Use Python 3.12. Run these commands from the repository root with the included public-use extracts.

```bash
python -m pip install -r requirements.txt
python code/analyse.py
python code/figures.py
```

The procedure runs the paired descriptive analysis, joint multinomial and binary sensitivity models, inverse response weighting with and without insulin, ranked weight-influence diagnostics, MAR/MNAR imputation, and the common-state year comparison. The reference outputs include [the adjusted contrasts](results/adjusted_contrasts.csv), [nonresponse estimates](results/nonresponse_sensitivity.csv), [module response patterns](results/module_missingness.csv), and [ranked weight influence](results/ipw_ranked_influence.csv).

The model-parameter JSON files list fitted coefficients and their covariance matrices in the same parameter order. The analysis locally generates model caches used by the verification script.

## Reconstruct the extracts from official public-use data

CDC publishes the original archives on its [2025 BRFSS page](https://www.cdc.gov/brfss/annual_data/annual_2025.html) and [2022 BRFSS page](https://www.cdc.gov/brfss/annual_data/annual_2022.html). The downloader retrieves the pinned XPT archives and checks their SHA-256 hashes. The original XPT data are downloaded into the ignored `data/raw/` directory and require about 2 GB after extraction. See [source provenance and archive hashes](sources/download_manifest.json).

```bash
python code/fetch.py
python code/extract.py
python code/coding.py
python code/analyse.py
python code/verify.py
```

The Python verification independently recomputes key design-based estimates and checks extract variables against the original XPTs. The optional R verification script requires the survey, arrow, and nnet packages. See [the verification report](qa/FINAL_VERIFICATION.md) for checks completed for this release.

## Interpretation

The reported outcomes are service contacts. Their clinical adequacy and appropriateness cannot be determined from BRFSS. Estimates refer to eligible participating-state domains. Probability weighting and imputation require assumptions about unobserved responses.
