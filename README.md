# A Blood-Test Based Predictive Model for the Cognitive Decline in Probable Alzheimer’s Disease Patients

## About
This repository include the code and result of prediction model building, training, and fine-tuning.

## Dataset
The cohort we used to train the predictive model comes from Zhong lab’s previous study, *Presymptomatic Increase of an Extracellular RNA in Blood Plasma Associates with the Development of Alzheimer’s Disease*. The data is publicly available on [GEO: GSE136243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136243).

## Directory Structure

    ├── data
    │   └── ADRC_PHGDH_with_cogscore_noAD.csv
    ├── nohup.out
    ├── README.md
    ├── result
    │   ├── algorithm_selection
    │   │   ├── PT_version1.csv
    │   │   ├── PT_version2.csv
    │   │   └── PT_version3.csv
    │   ├── modeling
    │   └── parameter_fine_tuning.csv
    └── script
        ├── algorithm
        │   ├── CF.R
        │   ├── prediction_model_conditionalCF.R
        │   ├── prediction_model_currentCF.R
        │   ├── prediction_model_current+futureCF.R
        │   ├── PT_version1.R
        │   ├── PT_version2.R
        │   └── PT_version3.R
        ├── modeling_result_illustartion.R
        ├── parameter_fine_tuning.R
        ├── PT_algorithm_selection.R
        ├── set_individual_coordinates.R
        └── storing_variables.R


