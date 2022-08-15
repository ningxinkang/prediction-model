##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-12    
## File: modeling_result_illustartion.R          
## Functions: 
##  1. PT_checker(PT_angle, PT_span, PT_span_max, t, PHGDH, donor_id, REGTRYID)
##  2. CF_checker(CF_angle, CF_span, t, DRS, CF_start, 
##  positive_PT_end, CF_end, drop_out_span)
##  3. prediction_model(PT_span, PT_span_max, CF_span, PT_angle, CF_angle, 
##  drop_out_span, donors_PHGDH, donors_PHGDH_t, 
##  donors_DRS, donors_DRS_t, sequence)
## File Summary:
##  This file include the code of the prediction model (full test).
##  The functions allow user to input the longitudinal dataset and
##  parameters, and output with the prediction result, including
##  statistics (accuracy, positive predictive value, and negaitve
##  predictive value) and illustration.
## Adjustments:
##  Change name of parameters, add function headers.
##__________________________________________________________________

############Loading Libraries############
library(dplyr)
library(tidyr)