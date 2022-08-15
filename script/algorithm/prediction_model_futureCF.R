##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-13    
## File: prediction_model_futureCF.R          
## Functions: 
##  prediction_model(PT_span, PT_span_max, CF_span, PT_angle, CF_angle, 
##  drop_out_span, donors_PHGDH, donors_PHGDH_t, 
##  donors_DRS, donors_DRS_t, sequence)
## File Summary:
##  This file include the code of the prediction model.
##  The functions allow user to input the longitudinal dataset and
##  parameters, and output with the prediction result, including
##  statistics (accuracy, positive predictive value, and negaitve
##  predictive value) and illustration.
## Note:
##  Cognitive score phase checker: future
##__________________________________________________________________


############Loading Libraries############
library(dplyr)
library(tidyr)

#' A predictive-model-based screening test
#' 
#' We use accuracy, PPV, and NPV to evaluate the predictive model’s
#' performance. When an individual with real cognitive decline is tested
#' positive, we say it is a true positive result as the model 
#' successfully predicts the presence of the cognitive decline in 
#' that individual. If that same individual is tested negative, 
#' meaning that the model fails to predict the presence of the 
#' cognitive decline, so, we say it is a false negative result. 
#' Similarly, given another individual without real cognitive decline 
#' who is tested negative, we categorize it as a true negative result 
#' as the predictive model successfully predicts the absence of 
#' the cognitive decline in that individual; if the same individual 
#' is tested positive, meaning that the model fails to identify the 
#' absence of the cognitive decline, the test is thus a false 
#' positive test.
#' 
#' @param PT_span A parameter, defines the minimum time span between two comparing
#' PHGDH level.
#' @param PT_span_max A parameter, defines the maximum time span between two 
#' comparing PHGDH level.
#' @param CF_span A parameter, defines the minimum time span between two comparing
#' cognitive level.
#' @param PT_angle A parameter, defines the minimum angle between two comparing
#' PHGDH level.
#' @param CF_angle A parameter, defines the minimum angle between two comparing
#' cognitive score.
#' @param drop_out_span The threshold evaluated the time from the end of the 
#' positive PT period to the drop-out date
#' @param donors_PHGDH A list of lists of double, including the TPM-normalized 
#' exRNA PHGDH level being compared of all the individuals in the dataset.
#' @param donors_PHGDH_t A list of lists of double, including the years to the 
#' drop-out date for each blood draw clinical visits of all the individuals in 
#' the dataset.
#' @param donors_DRS A list of lists of int, including the DRS score being 
#' compared of all the individuals in the dataset.
#' @param donors_DRS_t A list of lists of double, including the years to the 
#' drop-out date for each cognitive test clinical visits of all the individuals
#' in the dataset.
#' @param sequence A list of donors' ids, so that the predictive model follow 
#' this sequence
#' @return A list containing:
#' full_accuracy: accuracy per full test
#' full_PPV: positive predictive value per full test
#' full_NPV: negative predictive value per full test
#' simple_accuracy: accuracy per simple test
#' simple_PPV: positive predictive value per simple test
#' simple_NPV: negative predictive value per simple test
#' simple_total_number: total number of simple tests conducted
#' simple_true_positive: number of true positives
#' simple_true_negative: number of true negatives

prediction_model <- 
  function(PT_span, PT_span_max, CF_span, PT_angle, CF_angle, drop_out_span,
           donors_PHGDH, donors_PHGDH_t, donors_DRS, donors_DRS_t, sequence){
    
    ## Setup variables for model performance evaluation of both full test 
    ## and simple test
    full_true_positive = 0
    full_false_positive = 0
    full_true_negative = 0
    full_false_negative = 0
    full_accuracy = 0
    full_PPV = 0
    full_NPV = 0
    
    simple_true_positive = 0
    simple_true_negative = 0
    simple_accuracy = 0
    simple_PPV = 0
    simple_NPV = 0
    simple_positive = 0
    simple_total_number = 0
    
    ## Start looping through individuals in the dataset
    for (x in 1:length(donors_PHGDH)) {
      ## reset to empty for each individual
      PT_end <<- PT_end[0,]
      
      ## Label exRNA PHGDH level
      PT_result = PT_checker(PT_angle, PT_span, PT_span_max, 
                             donors_PHGDH_t[[x]], donors_PHGDH[[x]], 
                             as.character(sequence[x,2]),as.integer(sequence[x,1]))
      simple_total_number = simple_total_number + as.numeric(PT_result[[1]])
      num_simple_test <<- num_simple_test %>% 
        rbind(data.frame(REGTRYID = sequence[x,1], 
                         donor_id_alias = sequence[x,2], 
                         num = as.numeric(PT_result[[1]])))
      if (length(PT_result) == 1 & nrow(PT_end) > 0) {
        ## If there is no positive PT period, meaning the individual does not
        ## have increase in exRNA PHGDH level.
        
        ## Use 'count' to keep track of the number of positive CF periods of
        ## the tested individual
        count = 0
        for (i in 1:nrow(PT_end)){
          CF_result = CF_checker(CF_angle, CF_span, 
                                 donors_DRS_t[[x]], donors_DRS[[x]],
                                 donors_PHGDH_t[[x]][[1]], NA, 
                                 PT_end[i,3], drop_out_span)
          count = count + length(CF_result)
          if ((length(CF_result)) != 0){
            ## Record the decreasing CF period and drop-out period
            for (y in 1:length(CF_result)){
              if (CF_result[[y]][[3]]== "o"){
                CF_positives <<- CF_positives %>% 
                  rbind(data.frame(REGTRYID = sequence[x,1], 
                                   donor_id_alias = sequence[x,2], 
                                   start = CF_result[[y]][[1]], 
                                   end = CF_result[[y]][[2]]))
              } else {
                drop_out_record <<- drop_out_record %>% 
                  rbind(data.frame(REGTRYID = sequence[x,1], 
                                   donor_id_alias = sequence[x,2], 
                                   start = CF_result[[y]][[1]], 
                                   end = CF_result[[y]][[2]]))
              }
            }
          }
        }
        if (count == 0) {
          ## If there is NO positive CF period, meaning the individual is not 
          ## having cognitive decline. Taken together, those are TNs.
          full_true_negative = full_true_negative + 1
          simple_true_negative = simple_true_negative + as.numeric(PT_result[[1]])
          annotation <<- annotation %>% 
            rbind(data.frame(REGTRYID = sequence[x,1], 
                             donor_id_alias = sequence[x,2], 
                             label = "TN"))
        } else {
          ## If there is positive CF period, meaning the individual is having
          ## cognitive decline. Taken together, those are FNs.
          full_false_negative = full_false_negative + 1
          annotation <<- annotation %>% 
            rbind(data.frame(REGTRYID = sequence[x,1], 
                             donor_id_alias = sequence[x,2], 
                             label = "FN"))
        }
      }
      else if (length(PT_result) > 1 & nrow(PT_end) > 0){
        ## if there is positive PT period found, meaning the individual
        ## is found to have increase in exRNA PHGDH level
        
        simple_positive = simple_positive + 1
        PT_positives <<- PT_positives %>% 
          rbind(data.frame(REGTRYID = sequence[x,1], 
                           donor_id_alias = sequence[x,2], 
                           start = PT_result[[2]][[1]], 
                           end = PT_result[[2]][[2]]))
        CF_result = CF_checker(CF_angle, CF_span, 
                               donors_DRS_t[[x]], donors_DRS[[x]],
                               donors_PHGDH_t[[x]][[1]], PT_result[[2]][[2]], 
                               NA, drop_out_span)
        if (length(CF_result) != 0){
          ## If there is decreasing CF period found, meaning the
          ## individual is having cognitive decline. Taken together,
          ## those are TPs.
          simple_true_positive = simple_true_positive + 1
          full_true_positive = full_true_positive + 1
          annotation <<- annotation %>% 
            rbind(data.frame(REGTRYID = sequence[x,1], 
                             donor_id_alias = sequence[x,2], 
                             label = "TP"))
          for (z in 1:length(CF_result)){
            ## record the positive CF period and drop-out period
            if (CF_result[[z]][[3]]== "o"){
              CF_positives <<- CF_positives %>% 
                rbind(data.frame(REGTRYID = sequence[x,1], 
                                 donor_id_alias = sequence[x,2], 
                                 start = CF_result[[z]][[1]], 
                                 end = CF_result[[z]][[2]]))
            } else {
              drop_out_record <<- drop_out_record %>% 
                rbind(data.frame(REGTRYID = sequence[x,1], 
                                 donor_id_alias = sequence[x,2], 
                                 start = CF_result[[z]][[1]], 
                                 end = CF_result[[z]][[2]]))
            }
          }
        } else {
          ## If there is NO decreasing CF period found, meaning the
          ## individual is not having cognitive decline. Taken together,
          ## those are FPs.
          full_false_positive = full_false_positive + 1
          annotation <<- annotation %>% 
            rbind(data.frame(REGTRYID = sequence[x,1], 
                             donor_id_alias = sequence[x,2], 
                             label = "FP"))
        }
      }
    }
    
    ## Calculate performance
    full_accuracy = (full_true_positive + full_true_negative)/
      (full_true_positive + full_false_positive + 
         full_true_negative + full_false_negative)
    full_PPV = full_true_positive/(full_true_positive + full_false_positive)
    full_NPV = full_true_negative/(full_true_negative + full_false_negative)
    
    
    simple_accuracy = (simple_true_positive + simple_true_negative)/
      (simple_total_number)
    simple_PPV = simple_true_positive/simple_positive
    simple_NPV = simple_true_negative/(simple_total_number - simple_positive)
    
    return(c(full_accuracy, full_PPV, full_NPV, 
             simple_accuracy, simple_PPV, simple_NPV, 
             simple_total_number, simple_true_positive, simple_true_negative))
  }
