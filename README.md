# Heart_Failure_CSA

This is the repository for *Predicting survival time for critically ill patients with heart failure using conformalized survival analysis by Wang et al. (2024)*.

- `cfsurv_paper-main` contains the package from *Conformalized Survival Analysis by Candès et al. (2023)*, `R_train_predict` contains functions that we modified for our purposes.
- `Cleaned Dataset` contains csv files of the data that we filtered and cleaned from the MIMIC-IV dataset, the May10 version are the main data that we used, one is for patients with heart failure, nHF is for patients without heart failure. The June11 version is the high-dimensional dataset we used in the paper, containing 100 more medications as predictors.
- `Code` contains the code files we used to clean the data and fit the model.
- `Rshiny App` contains the files for our [online survival lower bound calculator](https://username434.shinyapps.io/Heart_failure_conformalized_survival_analysis/).
