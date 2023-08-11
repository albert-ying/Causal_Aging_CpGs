# Causal Epigenetic Age Uncouple Damage and Adaptation

This repository contains the code and related functions developed for the research paper "Causal Epigenetic Age: Uncouple Damage and Adaptation".
The paper delves into the epigenetic alterations associated with aging and how these changes impact damage and adaptation at the molecular level.

## Directory Structure

```plaintext
.
├── Code: Contains all the primary scripts used for the analyses detailed in the paper.
├── Functions: Houses utility functions and scripts that support the main analysis.
└── README.md: This file.
```

### Code

- `0.2_genetic_correlation_HDL.R`: Genetic correlation analysis for HDL.
- `0.3_genetic_correlation_ldsc.R`: LD score regression for genetic correlation.
- `1_GoDMC_meQTL_aging_2SMR.R`: Two-sample Mendelian Randomization for meQTL and aging using GoDMC data.
- `2_meQTL_pwcoco.R`: Pairwise conditional concordance analysis for meQTLs.
- `3_GoDMC_enrichment_analysis.R`: Enrichment analysis using GoDMC dataset.
- `4_EWAS_data_preprocessing.R`: Pre-processing of the EWAS data.
- `5_ewas_catalog_analysis.R`: Comprehensive analysis of the EWAS catalog.
- `6-ewas_enrichment.R`: Enrichment analysis for EWAS hits.
- `7_individual_dataset_analysis.R`: Detailed analysis of individual datasets.
- `8_GoDMC_clock.R`: Analyses related to the GoDMC epigenetic clock.

### Functions

- `FUN_correlationPlot.R`: A function for creating correlation plots.
- `FUN_model_helper.R`: Helper functions for modeling purposes.
- `FUN_mol_mr.R`: Functions related to molecular Mendelian Randomization.
- `FUN_plot_wrapper.R`: Plotting wrapper functions.
- `ldsc_binary_annot_QTL.sh`: Script for LDSC binary annotation for QTLs.
- `make_ldsc_binary_annot.R`: R script to prepare binary annotations for LDSC.

## Citation

```plaintext
Ying, K. et al. Causal Epigenetic Age Uncouples Damage and Adaptation. 2022.10.07.511382 Preprint at https://doi.org/10.1101/2022.10.07.511382 (2022).

```

## License

```plaintext
This project is licensed under the MIT License - see the LICENSE.md file for details.
```
