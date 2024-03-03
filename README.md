# Bachelor's Thesis (HS2022)

## Table of Contents

- [Introduction](#introduction)
- [Contributors](#contributors)
- [Abstract](#abstract)
- [Structure](#repository-structure)
- [Usage](#usage)
- [Acknowledgements](#acknowledgements)

## Introduction

In this repository you can find the code of my Bachelor's thesis **"Portfolio Value at Risk Forecasting with GARCH-Type Models"** where I compared different GARCH-type models with respect to their ability to forecast the Value at Risk (VaR) of an equally weighted portfolio consisting of ten US large caps stocks. Most of my work revolved around implementing a *factor copula-DCC-(N)GARCH* model as described in Fortin et al. (2022).

## Contributors

* Jan Heinrich Schlegel (Author)

## Abstract:

This thesis examines the value at risk (VaR) forecasting ability of various univariate and multivariate
models for a long equity portfolio. All of the considered models involve a generalized autoregressive conditional heteroskedasticity (GARCH)-type structure. The resulting forecasts are checked for
desirable properties using violation-based backtests and compared in terms of predictive ability. We
find that the VaR forecasts of almost all univariate models are inadequate, while the multivariate
models have few problems passing these backtests. However, we do not find evidence that the multivariate models systematically outperform their univariate counterparts with regards to predictive
accuracy, or vice versa.

## Repository Structure
```bash
├── Bachelor Thesis Code
│   ├── BachelorThesis.Rproj
│   ├── Data
│   │   ├── FamaFrenchFactorsDaily.csv
│   │   ├── FFCFactors.csv
│   │   ├── MomentumFactorDaily.csv
│   │   ├── PortfolioPlrets.csv
│   │   └── StockPlrets.csv
│   ├── Results
│   │   ├── Plots
│   │   │   ├── ACF_Plots_Factors.pdf
│   │   │   ├── ACF_Plots_Shares.pdf
│   │   │   ├── Chisq_Plot_OLS.pdf
│   │   │   ├── Factors.pdf
│   │   │   ├── MD_Chisq_Plot.pdf
│   │   │   ├── PF_plots.pdf
│   │   │   └── Rplot.pdf
│   │   └── VaR
│   │       ├── COMFORT_MVG_CCC_GJR_VaR.csv
│   │       ├── COMFORT_MVG_CCC_sGARCH_VaR.csv
│   │       ├── Fortin_cop_norm_NGARCH.csv
│   │       ├── Fortin_cop_norm_sGARCH.csv
│   │       ├── Fortin_cop_skewt_sGARCH.csv
│   │       ├── Fortin_cop_t_NGARCH.csv
│   │       ├── Fortin_cop_t_sGARCH.csv
│   │       ├── Multi_Normal_DCC_GARCH_Matlab.csv
│   │       ├── Multi_Normal_DCC_GARCH_R.csv
│   │       ├── Uni_EWMA_VaR.csv
│   │       ├── Uni_MN_2_2_GARCH.csv
│   │       ├── Uni_MN_3_3_GARCH.csv
│   │       ├── Uni_Normal_GARCH_VaR.csv
│   │       ├── Uni_skewt_GJR_GARCH.csv
│   │       ├── Uni_Skewt_NGARCH.csv
│   │       └── Uni_t_GJR_GARCH.csv
│   └── Scripts
│       ├── Backtesting_Functions.R
│       ├── Backtesting.R
│       ├── Columnwise_Sum.cpp
│       ├── COMFORT_rolling_window.m
│       ├── EDA.R
│       ├── FortinCopulaGARCH.R
│       ├── Skew_t_Copula.R
│       └── Univariate_Models.R
├── Documents
│   ├── BA_Thesis_ExecutiveSummary.pdf
│   ├── BA_Thesis.pdf
│   └── BA_Thesis_Tables.pdf
└── README.md
```

## Usage

- `Data/`: Contains the datasets used.
- `Scripts/`: Includes all R and Python scripts for data processing, model fitting, and analysis.
- `Results/`: Contains all plots as well as the Value at Risk (VaR) predictions in the form of CSV files.
- `Documentation/`: Contains the thesis document and any additional notes.

## Acknowledgements

Special thanks to Prof. Dr. Marc Paolella and Soros Chitsiripanich for supervision, guidance and invaluable feedback.
