##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-13    
## File: PT_version2.R          
## Functions: 
##  PT_checker(PT_angle, PT_span, PT_span_max, t, PHGDH, donor_id, REGTRYID)
## File Summary:
##  This file include the code of exRNA level data labeling.
##  The functions allow user to input the longitudinal dataset and
##  parameters, and output with the period identified as molecular
##  level increasing period.
## Notes:
##  Check PT history (yes)
##  Check time condition (during the test)
##__________________________________________________________________


############Loading Libraries############
library(dplyr)
library(tidyr)

#' Data labeling: Check the increase period in exRNA PHGDH level of an individual.
#' 
#' A predictive screening test is conducted after each blood test with biomarker
#'  fluctuations recorded within a specific period. In other words, the 
#'  condition for a screening test to work is that at least one other blood 
#'  draw falls within a specific time the current blood draw is taken. This 
#'  specific time range is restricted by two parameters: maximum PT span 
#'  'PT_span_max' and minimum PT span 'PT_span'. 
#' 
#' The exRNA PHGDH level is increasing when within a specific time, two exRNA 
#'  PHGDH levels l1 and l2 , sequenced from two blood tests taken at t1 and t2,
#'  are undergoing sequential increase, and the time interval and the increasing
#'  angle surpasses the conditions we defined. The increasing angle is 
#'  calculated by arctan ((l2 − l1)/(t2 − t1)), and we used the parameter
#'  minimum PT angle 'PT_angle' to set the boundary. 
#' 
#' When there is a increasing period being recognised, the algorithm will stop 
#' testing following records, which means that there will be at most one 
#' positive increasing period. 
#' 
#' @param PT_angle A parameter, defines the minimum angle between two comparing
#' PHGDH level.
#' @param PT_span A parameter, defines the minimum time span between two comparing
#'PHGDH level.
#' @param PT_span_max A parameter, defines the maximum time span between two 
#' comparing PHGDH level.
#' @param t A list of double, represents the years to the drop-out date for 
#' each blood draw clinical visits of the input individual.
#' @param PHGDH A list of double, represents the TPM-normalized exRNA PHGDH 
#' level being compared of the input individual.
#' @param donor_id The donor ID of the individual.
#' @param REGTRYID The registry ID of the individual in ADRC.
#' @return A list, with the first object an int, represents the number of 
#' simple test conducted on the input individual; the following objects are
#' lists of lists in the form of (start, end), represent the time span tested
#' to be positive for the input individual.
#' @examples
#' PT_checker(1.5, 0.4, c(-12.0, -11.0, -10.0, -9.0, -8.0, -6.5, -5.5, -4.0, -3.0, -1.5, 0.0), c(10.17161, 5.47887, 8.37505, 3.071834, 42.142075, 68.835006), "C_1_F", "14058")

PT_checker <- 
  function(PT_angle, PT_span, PT_span_max, t, PHGDH, donor_id, REGTRYID){
    ## Return variables record the number of test conducted and the time span of 
    ## the increasing period. 
    period = list()
    test = 0
    
    ## A boolean worked as a stop sign. When a increasing period is identified, 
    ## the stop sign will become TRUE, ending the testing loop.
    stop = FALSE
    
    ## Start looping through the blood test result
    ## Determine the end of comparing time span
    for (i in 2:length(PHGDH)){
      t_earliest = t[i] - PT_span_max
      t_latest = t[i] - PT_span
      
      ## Number of PT increment when there is another blood draw conducted
      ## within the time span
      if (sum(t >= t_earliest & t <= t_latest)>0){
        test = test + 1
      }
      
      ## Start comparing
      for (j in 1:(i-1)){
        if (t[j] >= t_earliest & t[j] <= t_latest){
          t_i_start = j+1
          t_i_end = i
          
          ## Evaluate if the increasing slope reach the threshold
          if ((PHGDH[[i]]-PHGDH[[j]])/(t[[i]]-t[[j]]) >= tanpi(PT_angle)){
            for (k in t_i_start:t_i_end){
              ## If the exRNA PHGDH level decrease during the time span, 
              ## PT is considered negative
              if (PHGDH[k] < PHGDH[k-1]){
                break
              } else {
                ## Check fluctuation until reach the last blood draw
                if (k == i){
                  ## Record the increasing positive PT period
                  stop = TRUE
                  period <- period %>% 
                    append(list(list(t[[j]],t[[i]])))
                  test_record <<-test_record %>% 
                    rbind(data.frame(REGTRYID = REGTRYID, 
                                     donor_id_alias = donor_id, 
                                     x1 = t[j], x2 = t[i], 
                                     y1 = PHGDH[j], y2 = PHGDH[i], 
                                     type = "positive"))
                  PT_end <<- PT_end %>% 
                    rbind(data.frame(REGTRYID = REGTRYID, 
                                     donor_id_alias = donor_id, 
                                     mark = t[i]))
                }
              }
            }
            if (stop) {break}
            ## Record the negative PT period that fail to pass 
            ## fluctuation restriction
            test_record <<-test_record %>% 
              rbind(data.frame(REGTRYID = REGTRYID, donor_id_alias = donor_id, 
                               x1 = t[j], x2 = t[i], 
                               y1 = PHGDH[j], y2 = PHGDH[i], type = "negative"))
            PT_end <<- PT_end %>% 
              rbind(data.frame(REGTRYID = REGTRYID, donor_id_alias = donor_id, 
                               mark = t[i]))
          } else {
            ## record the negative PT period that fail to pass slope restriction
            test_record <<-test_record %>% 
              rbind(data.frame(REGTRYID = REGTRYID, donor_id_alias = donor_id, 
                               x1 = t[j], x2 = t[i], 
                               y1 = PHGDH[j], y2 = PHGDH[i], type = "negative"))
            PT_end <<- PT_end %>% 
              rbind(data.frame(REGTRYID = REGTRYID, donor_id_alias = donor_id, 
                               mark = t[i]))
          }
        }
      }
      ## Stop PT test as there is already a positive PT period
      if (stop) {break}
    }
    PT_end <<- distinct(PT_end)
    return(c(test,period))
  }
