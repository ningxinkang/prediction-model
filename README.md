# A Blood-Test Based Predictive Model for the Cognitive Decline in Probable Alzheimer’s Disease Patients

## About
This repository include the code and result of prediction model building, training, and fine-tuning.

## Dataset
The cohort we used to train the predictive model comes from Zhong lab’s previous study, *Presymptomatic Increase of an Extracellular RNA in Blood Plasma Associates with the Development of Alzheimer’s Disease*. The data is publicly available on [GEO: GSE136243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136243).

## Directory Structure
'''
.
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
│   │   ├── Picture1_top_candidates.png
│   │   ├── Picture2(1)_ROC_minPT.png
│   │   ├── Picture2(2)_ROC_minPT.png
│   │   ├── Picture3_major-change_model.png
│   │   ├── Picture4_comprehensive_model.png
│   │   ├── Picture5_PTcheckers_info.png
│   │   ├── Picture6_PTcheckers_performance.png
│   │   ├── PictureS1_drop-out_span.png
│   │   ├── PictureS2_DRS_vs_MMSE.png
│   │   ├── PictureS3_MMSE_major-change_model.png
│   │   ├── PictureS4_normalization_methods.png
│   │   ├── PictureS5_strandedness.png
│   │   ├── PictureS6_age&sex_baseline.png
│   │   └── PictureS7_gender.png
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
'''

