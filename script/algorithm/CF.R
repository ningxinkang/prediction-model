##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-13    
## File: CF.R          
## Functions: 
##  CF_checker(CF_angle, CF_span, t, DRS, CF_start, 
##  positive_PT_end, CF_end, drop_out_span)
## File Summary:
##  This file include the code of cognitive function data labeling.
##  The functions allow user to input the longitudinal dataset and
##  parameters, and output with the period identified as cognitive
##  decreasing period.
##__________________________________________________________________

############Loading Libraries############
library(dplyr)
library(tidyr)

#' Data Labeling: Check the decrease period in cognitive score of an individual.
#' 
#' In a given range of time, we define an individual as having cognitive 
#' decline when the two cognitive scores s1 and s2 , collected at t1 and t2
#' correspondingly, are decreasing along the time course, with the decreasing 
#' angle and the time interval passing the conditions we defined. The 
#' decreasing angle is calculated by arctan((s2 − s1)/(t2 −t1)).
#' 
#' To search for such cognitive decline, every two cognitive scores within 
#' the given time range are compared to see if all the following three 
#' conditions are met: 1) the time interval between two cognitive test dates
#' is longer than the minimum decreasing span defined; 2) cognitive scores do 
#' not increase within the time; 3) the decreasing angle between the two 
#' cognitive test scores is larger than the minimum decreasing angle defined. 
#' 
#' @param CF_angle A parameter, defines the minimum angle between two comparing
#' cognitive score.
#' @param CF_span A parameter, defines the minimum time span between two comparing
#' cognitive level.
#' @param t A list of double, represents the years to the drop-out date for 
#' each cognitive test clinical visits of the input individual.
#' @param DRS A list of int, represents the cognitive test scores being 
#' compared of the input individual.
#' @param CF_start The time indicates the start time to evaluate cognitive 
#' function. Which is always the first cognitive tests obtained during the study.
#' @param positive_PT_end The time indicates the end of positive PT period. If no 
#' positive PT period, than this variable is N/A. This variable is used to 
#' evaluate if the 'out_put_span' condition is met.
#' @param CF_end The time indicates the end time to evaluate cognitive 
#' function. If PT is positive, the 'positive_PT_end' equals to the last available 
#' date of conducting the cognitive test; if PT is negative, tbe 'positive_PT_end' 
#' equals to the last date of the PT span.
#' @param drop_out_span The threshold evaluated the time from the end of the 
#' positive PT period to the drop-out date
#' @return A list of list in the form of (start, end, category), represent the 
#' positive CF span. If category = o, that is normal output; if category = d, 
#' that is the drop out period.
#' @examples
#' CF_checker(0.4, 4.5, c(-12.0, -11.0, -10.0, -9.0, -8.0, -6.5, -5.5, -4.0, -3.0, -1.5, 0.0), c(139, 136, 142, 137, 140, 140, 137, 143, 135, 138, 132), -6.5, -3.0, 0.0, 3.0)

CF_checker <- 
  function(CF_angle, CF_span, t, DRS, 
           CF_start, positive_PT_end, CF_end, drop_out_span){
    ## return variable
    period = list()
    
    ## Figure out the CF start time index
    CF_start_temp = CF_start
    start_index = which(t==as.numeric(CF_start))
    ## if there isn't a match index, select the closest score
    while (length(start_index) == 0){
      CF_start_temp = as.numeric(CF_start_temp) + 0.5
      start_index = which(t==as.numeric(CF_start_temp))
    }
    
    ## figure out the CF end time index
    CF_end_temp = CF_end
    end_index = double()
    if (is.na(CF_end)){
      end_index = length(DRS)
    } else {
      end_index = which(t==as.numeric(CF_end))
    }
    ## if there isn't a match index, select the closest score
    while (length(end_index) == 0){
      CF_end_temp = as.numeric(CF_end_temp) - 0.5
      end_index = which(t==as.numeric(CF_end_temp))
    }
    
    ## start CF testing
    if (start_index <= end_index){
      for (x in start_index:end_index){
        start_level = DRS[[x]]
        y = x + 1
        while (y <= end_index){
          next_level = DRS[[y]]
          prev_level = DRS[[y-1]]
          ## Check fluctuation, we want the cognitive score to continuously
          ## decreasing
          if (next_level - prev_level > 0){
            y = y + 1
            break
          }
          ## Check decreasing slope
          if (t[[y]]-t[[x]] >= CF_span 
              & (DRS[[y]]-DRS[[x]])/(t[[y]]-t[[x]]) <= tanpi(CF_angle)){
            ## edge case
            if (DRS[[y]]-DRS[[x]] == 0){
              break
            } else if (length(period) == 0){
              ## direct add positive CF period if no previous periods added
              period <- period %>% append(list(list(t[[x]],t[[y]],"o")))
            } else if (period[[length(period)]][[1]] == t[[x]]){
              ## if previous positive CF period has the same start time as the 
              ## current one, delete the last CF period
              period[[length(period)]] <- NULL
              period <- period %>% append(list(list(t[[x]],t[[y]],"o")))
              if (length(period) > 1){
                if (period[[length(period)-1]][[2]] == 
                    period[[length(period)]][[2]]){
                  ## if previous period's end time equals to the newly added end 
                  ## time, delete the newly added period, because we want a 
                  ## longer CF positive period.
                  period[[length(period)]] <- NULL
                }
              }
            } else if (period[[length(period)]][[2]] == t[[y]]){
              ## if previous positive CF period end equals to current end,
              ## ignore the newly found CF period
              break
            } else{
              period <- period %>% append(list(list(t[[x]],t[[y]],"o")))
            }
          }
          y = y + 1
        }
      }
    }
    
    ## Drop-out span evaluation
    if (!is.na(positive_PT_end)){
      if (0 - positive_PT_end <= drop_out_span){
        period <- period %>% append(list(list(positive_PT_end,0,"d")))
      }
    }
    return(period)
  }