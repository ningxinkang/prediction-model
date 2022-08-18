# A Blood-Test Based Predictive Model for the Cognitive Decline in Probable Alzheimer’s Disease Patients

## About
This repository include the code and result of prediction model building, training, and fine-tuning.

## Dataset
The cohort we used to train the predictive model comes from Zhong lab’s previous study, *Presymptomatic Increase of an Extracellular RNA in Blood Plasma Associates with the Development of Alzheimer’s Disease*. The data is publicly available on [GEO: GSE136243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136243).

## Directory Structure

    .
    ├── data                                         <- includes the sample formatted longitudinal dataset     
    │   └── ADRC_PHGDH_with_cogscore_noAD.csv
    ├── README.md
    ├── result                                       <- includes model performance, graph, etc.                                   
    │   ├── algorithm_selection                      <- performances of the model trained with three versions of PT checker
    │   │   ├── PT_version1.csv
    │   │   ├── PT_version2.csv
    │   │   └── PT_version3.csv
    │   ├── modeling                                 <- illustartions of the model, and graphs for variable comparison
    │   └── parameter_fine_tuning.csv                <- accuracy, PPV, NPV, and other specific informations of models trained with 11664 candidate parameter combinations              
    └── script                                       <- includes scripts for the prediction model algorithm and further analysis
        ├── algorithm                                <- scripts for the prediction model algorithm
        │   ├── CF.R                                 <- cognitive function checker, label cognitive score
        │   ├── prediction_model_conditionalCF.R                            <- function that combine the two labeling functions
        │   ├── prediction_model_currentCF.R                              
        │   ├── prediction_model_current+futureCF.R                        
        │   ├── PT_version1.R                                <- predictive test checker, label the exRNA PHGDH level
        │   ├── PT_version2.R                             
        │   └── PT_version3.R                   
        ├── modeling_result_illustartion.R                   <- code for modeling result illustration
        ├── parameter_fine_tuning.R                          <- code for setting loops for parameter fine tuning
        ├── PT_algorithm_selection.R                         <- code that include the candidate parameter sets for three versions of PT checkers
        ├── set_individual_coordinates.R                     <- an environment code referenced from online resources, allowing users to set individual coordinates for each facets in ggplot2
        └── storing_variables.R                              <- setup variables for information storing during the model training process



**Last update: Aug 18, 2022**

