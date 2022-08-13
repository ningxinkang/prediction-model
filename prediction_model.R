##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-12    
## File: prediction_model.R          
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
library(ggplot2)
library(lemon)
library(gtools)
library(stringr)

###############FUNCTIONS################

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

##########Set coordinates for individual ggplot2 facets############
## This code block is directly copied from:
## https://gist.github.com/burchill/d780d3e8663ad15bcbda7869394a348a

UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  panel_ranges = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    layout
  },
  
  setup_panel_params = 
    function(self, scale_x, scale_y, params = list()) {
    if (!is.null(self$panel_ranges) & length(self$panel_ranges) !=
        self$num_of_panels)
      stop("Number of panel ranges does not equal the number supplied")
    
    train_cartesian <- 
      function(scale, limits, name, given_range = NULL) {
      if (is.null(given_range)) {
        expansion <- 
          ggplot2:::default_expansion(scale, expand = self$expand)
        range <- 
          ggplot2:::expand_limits_scale(scale, expansion, 
                                        coord_limits =
                                          self$limits[[name]])
      } else {
        range <- given_range
      }
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, range),
        sec = ggplot2:::view_scale_secondary(scale, limits, range),
        arrange = scale$axis_order(),
        range = range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    cur_panel_ranges <- self$panel_ranges[[self$panel_counter]]
    if (self$panel_counter < self$num_of_panels)
      self$panel_counter <- self$panel_counter + 1
    else
      self$panel_counter <- 1
    
    c(train_cartesian(scale_x, self$limits$x, "x", cur_panel_ranges$x),
      train_cartesian(scale_y, self$limits$y, "y", cur_panel_ranges$y))
  }
)

coord_panel_ranges <- 
  function(panel_ranges, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords, panel_ranges = panel_ranges, 
                   expand = expand, default = default, clip = clip)
}

################Loading datasets################
#_______________________________________________
#  Dataset format:
#    Rows: 
#       Samples
#    Columns:
#       sample_id_alias: sample id
#       x.axis: time to the drop-out date
#       TPM: TPM-normalized exRNA PHGDH level
#       REGTRYID: donor id in ADRC system
#       donor_id_alias: donor id
#       donor_group: C/N (converter/control)
#       convert_date: the year diagnosed as AD
#       MMSE: Mini-Mental State Exam (MMSE) score
#       DRS: Dementia Rating Scale-2 (DRS) score
#________________________________________________

## Prepare longitudinal list for each individual
## and their corresponding time stamp

PHGDH <- 
  read.csv("ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  drop_na(TPM) %>% arrange(x.axis) %>% group_by(donor_id_alias)
donors_PHGDH = split(PHGDH$TPM,PHGDH$donor_id_alias)
donors_PHGDH_t = split(PHGDH$x.axis,PHGDH$donor_id_alias)
PHGDH <- PHGDH[,1:5]
colnames(PHGDH)[3] <- "value"
colnames(PHGDH)[2] <- "time"

DRS <- 
  read.csv("ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  drop_na(DRS) %>% arrange(x.axis) %>% group_by(donor_id_alias)
donors_DRS = split(DRS$DRS,DRS$donor_id_alias)
donors_DRS_t = split(DRS$x.axis,DRS$donor_id_alias)
DRS <- DRS[,c(1,2,9,5,4)]
colnames(DRS)[3] <- "value"
colnames(DRS)[2] <- "time"

sequence <- distinct(PHGDH[,4:5]) %>% arrange(donor_id_alias)
convert_dates <- 
  read.csv("ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  select(REGTRYID,donor_id_alias,convert_date) %>% 
  arrange(donor_id_alias) %>% 
  distinct()

################Set up parameters################
PT_span <- 1.5
PT_span_max <- 4.5
CF_span <- 0.5
PT_angle <- 0.4
CF_angle <- -0.45
drop_out_span <- 3

################Set up storing variables################
x <- c("REGTRYID", "donor_id_alias", "start", "end")

## A dataframe, record the positive PT periods.
PT_positives <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(PT_positives) <- x

## A dataframe, record the positive CF periods.
CF_positives <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(CF_positives) <- x

## A dataframe, record the true positive current CF periods.
CF_current_TPs <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(CF_current_TPs) <- x

## A dataframe, record the true positive future CF periods.
CF_future_TPs <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(CF_future_TPs) <- x

## A dataframe, record the drop-out periods.
drop_out_record <-data.frame(matrix(ncol = 4, nrow = 0))
colnames(drop_out_record) <- x

## A dataframe, record the result of full test.
## TP/TN/FP/FN
annotation <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(annotation) <- c("REGTRYID", "donor_id_alias", "label")

## A dataframe, record the total number of simple test conducted
## for each individual
num_simple_test <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(num_simple_test) <- c("REGTRYID", "donor_id_alias", "num")

## A dataframe, record all the PT conducted
## for each individual
test_record <-data.frame(matrix(ncol = 7, nrow = 0))
colnames(test_record) <- 
  c("REGTRYID", "donor_id_alias", "x1", "x2", "y1","y2","type")

## A dataframe, record the end of PT
## for each individual
PT_end <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(PT_end) <- c("REGTRYID", "donor_id_alias", "mark")

## A dataframe, record the category of the TP CF decline period
## current/future/drop-out
decline_categorization <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(decline_categorization) <- 
  c("REGTRYID", "donor_id_alias", "label")

## A dataframe, record the time and PHGDH level of two ends of 
## positive PT periods.
PT_highlight <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(PT_highlight) <- c("sample_id_alias", "time", "value",
                            "REGTRYID","donor_id_alias")

## A dataframe, record the time and PHGDH level of two ends of 
## positive CF current periods.
CF_current_highlight <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(CF_current_highlight) <- c("sample_id_alias", "time", "value",
                                    "REGTRYID","donor_id_alias", "index")

## A dataframe, record the time and PHGDH level of two ends of 
## positive CF future periods.
CF_future_highlight <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(CF_future_highlight) <- c("sample_id_alias", "time", "value",
                                   "REGTRYID","donor_id_alias", "index")

################PREDICTION TEST#####################
performance <- prediction_model(PT_span, PT_span_max, CF_span, PT_angle,
                                CF_angle, drop_out_span, donors_PHGDH,
                                donors_PHGDH_t, donors_DRS, 
                                donors_DRS_t, sequence)

####################Manage result for illustration#################

##________________simple PT result categorization_________________##
## In the figure illustration, we want to color-code the curve that
## represents the result of every single simple PT conducted.
test_record <- 
  merge(test_record, annotation, by = "donor_id_alias", all.x = TRUE)
test_record <- test_record[-8]
## 'label' is the prediction result of the full test. 
## 'type' is the PT labeling result.
test_record <- test_record %>% 
  mutate(category = case_when(type == "positive" & label == "TP" ~ "TP",
                              type == "positive" & label == "FP" ~ "FP",
                              type == "negative" & label == "TN" ~ "TN",
                              type == "negative" & label == "TP" ~ "FN",
                              type == "negative" & label == "FN" ~ "FN",
                              type == "negative" & label == "FP" ~ "TN"))
## 'category' is the simple test prediction result
## 'final' is the simple test prediction category, 
## which will be used for color-code.
test_record <- test_record %>% 
  mutate(final = case_when(category == "TP" | category == "TN" ~ "right",
                           category == "FP" | category == "FN" ~ "wrong"))

##________________full test categorization_________________##
## 'label' is the full test prediction result
## 'for_color' is the prediction category, 
## which will be used for color-code.
annotation <- annotation %>% 
  mutate(for_color = case_when(label == "TP" | label == "TN" ~ "right", 
                               label == "FP" | label == "FN" ~ "wrong"))


##________________PT highlight area_________________##
for (i in 1:nrow(sequence)){
  index = as.list(which(PT_positives$donor_id_alias == 
                          as.character(sequence[i,2])))
  if (length(index) != 0){
    for (m in 1:length(index)){
      line = index[[m]]
      ## fetch the value and time for points within the positive PT period.
      new = PHGDH %>% 
        dplyr::filter(donor_id_alias == as.character(sequence[i,2]) & 
                        time >= PT_positives[line,3] & 
                        time <= PT_positives[line,4])
      PT_highlight <- rbind(PT_highlight,new)
    }
  }
}

##________________CF highlight area_________________##
## Step 1: split CF positive periods into future and current
for (i in 1:nrow(CF_positives)){
  PT_positive_index = as.list(which(PT_positives$donor_id_alias == 
                           as.character(CF_positives[i,2])))
  PT_index = as.list(which(test_record$donor_id_alias == 
                                as.character(CF_positives[i,2])))
  if (length(PT_positive_index) != 0){
    ## if there is positive PT period for the individual
    if (PT_positives[PT_positive_index[[1]],4] <= CF_positives[i,3]){
      CF_future_TPs <<- CF_future_TPs %>% 
        rbind(data.frame(REGTRYID = CF_positives[i,1], 
                         donor_id_alias = CF_positives[i,2], 
                         start = CF_positives[i,3], 
                         end = CF_positives[i,4]))
    } else {
      CF_current_TPs <<- CF_current_TPs %>% 
        rbind(data.frame(REGTRYID = CF_positives[i,1], 
                         donor_id_alias = CF_positives[i,2], 
                         start = CF_positives[i,3], 
                         end = CF_positives[i,4]))
    }
  } else if (as.character(test_record[PT_index[[1]],8]) == "FN"){
    ## if there is NO positive PT period for the individual, but do
    ## have positive CF period.
    end = -Inf
    for (j in 1:length(PT_index)){
      if (test_record[PT_index[[j]],4] > end){
        end = test_record[PT_index[[j]],4]
      }
    }
    if (end <= CF_positives[i,3]){
      CF_future_TPs <<- CF_future_TPs %>% 
        rbind(data.frame(REGTRYID = CF_positives[i,1], 
                         donor_id_alias = CF_positives[i,2], 
                         start = CF_positives[i,3], 
                         end = CF_positives[i,4]))
    } else {
      CF_current_TPs <<- CF_current_TPs %>% 
        rbind(data.frame(REGTRYID = CF_positives[i,1], 
                         donor_id_alias = CF_positives[i,2], 
                         start = CF_positives[i,3], 
                         end = CF_positives[i,4]))
    }
  }
}

## Step 2: Manage the highlighted area of current CFs
## Use 'trace' to trace the rows that is covered within one highlited area
trace = 1
for (i in 1:nrow(sequence)){
  index = as.list(which(CF_current_TPs$donor_id_alias == 
                           as.character(sequence[i,2])))
  if (length(index) != 0){
    for (m in 1:length(index)){
      line = index[[m]]
      new = DRS %>% 
        dplyr::filter(donor_id_alias == as.character(sequence[i,2]) & 
                        time >= CF_current_TPs[line,3] & 
                        time <= CF_current_TPs[line,4])
      new['index'] = trace
      CF_current_highlight <- rbind(CF_current_highlight,new)
      trace = trace + 1
    }
  }
}

## Step 3: Manage the highlighted area of future CFs
## Use 'trace' to trace the rows that is covered within one highlited area
trace = 1
for (i in 1:nrow(sequence)){
  index = as.list(which(CF_future_TPs$donor_id_alias == 
                           as.character(sequence[i,2])))
  if (length(index) != 0){
    for (m in 1:length(index)){
      line = index[[m]]
      new = DRS %>% 
        dplyr::filter(donor_id_alias == as.character(sequence[i,2]) & 
                        time >= CF_future_TPs[line,3] & 
                        time <= CF_future_TPs[line,4])
      new['index'] = trace
      CF_future_highlight <- rbind(CF_future_highlight,new)
      trace = trace + 1
    }
  }
}

##________________Number of specific test conducted_________________##
full_current_TP = 0
full_future_TP = 0
full_dropout_TP = 0
full_current_FN = 0
full_future_FN = 0

simple_current_TP = 0
simple_future_TP = 0
simple_dropout_TP = 0
simple_current_FN = 0
simple_future_FN = 0

## full test decline categorization

## For one individual, there may be multiple positive CF period
## So a summary label are created here.
## current > future > drop-out
for (i in 1:nrow(annotation)){
  if (annotation[i,3] == "TP" | annotation[i,3] == "FN"){
    donor = annotation[i,2]
    current = as.list(which(CF_current_TPs$donor_id_alias == donor))
    if (length(current) == 0){
      future = as.list(which(CF_future_TPs$donor_id_alias == donor))
      if (length(future) == 0){
        ## CASE 1: drop-out
        decline_categorization <<- decline_categorization %>% 
          rbind(data.frame(REGTRYID = annotation[i,1], 
                           donor_id_alias = annotation[i,2], 
                           label = "drop-out span"))
        if(annotation[i,3] == "TP"){
          full_dropout_TP = full_dropout_TP + 1
        }
      } else {
        ## CASE 2: future decline
        decline_categorization <<- decline_categorization %>% 
          rbind(data.frame(REGTRYID = annotation[i,1], 
                           donor_id_alias = annotation[i,2], 
                           label = "future"))
        if(annotation[i,3] == "TP"){
          full_future_TP = full_future_TP + 1
        } else {
          full_future_FN = full_future_FN + 1
        }
      }
    } else {
      ## CASE 3: current decline
      decline_categorization <<- decline_categorization %>% 
        rbind(data.frame(REGTRYID = annotation[i,1], 
                         donor_id_alias = annotation[i,2], 
                         label = "current"))
      if(annotation[i,3] == "TP"){
        full_current_TP = full_current_TP + 1
      } else {
        full_current_FN = full_current_FN + 1
      }
    }
  }
}

## simple test number summary
for (i in 1:nrow(test_record)){
  ## test_record[i,9] is 'category',which is the simple test prediction result
  if (test_record[i,9] == "TP" | test_record[i,9] == "FN"){
    donor = test_record[i,1]
    current = as.list(which(CF_current_TPs$donor_id_alias == donor))
    if (length(current) == 0){
      future = as.list(which(CF_future_TPs$donor_id_alias == donor))
      if (length(future) == 0){
        ## CASE 1: drop-out
        if(test_record[i,9] == "TP"){
          simple_dropout_TP = simple_dropout_TP + 1
        }
      } else {
        ## CASE 2: future decline
        if(test_record[i,9] == "TP"){
          simple_future_TP = simple_future_TP + 1
        } else {
          simple_future_FN = simple_future_FN + 1
        }
      }
    } else {
      ## CASE 3: current decline
      if(test_record[i,9] == "TP"){
        simple_current_TP = simple_current_TP + 1
      } else {
        simple_current_FN = simple_current_FN + 1
      }
    }
  }
}

####################Result summary#################
caption <- paste("Details: accuracy_per_full_test =", performance[1],
                 "; PPV_per_full_test =", performance[2],
                 "; NPV_per_full_test =", performance[3],
                 "; num_of_current_TP_full_test =", full_current_TP,
                 "; num_of_future_TP_full_test =", full_future_TP,
                 "; num_of_drop-out_TP_full_test =", full_dropout_TP,
                 "; num_of_current_FN_full_test =", full_current_FN,
                 "; num_of_future_FN_full_test =", full_future_FN,
                 "; accuracy_per_simple_test =", performance[4],
                 "; PPV_per_simple_test =", performance[5],
                 "; NPV_per_simple_test =", performance[6],
                 "; total_num_of_simple_tests =", performance[7],
                 "; avg_simple_tests_per_person =", performance[7]/20,
                 "; num_of_TP_simple_test = ", performance[8],
                 "; num_of_TN_simple_test = ", performance[9],
                 "; num_of_current_TP_simple_test =", simple_current_TP,
                 "; num_of_future_TP_simple_test =", simple_future_TP,
                 "; num_of_drop-out_TP_simple_test =", simple_dropout_TP,
                 "; num_of_current_FN_simple_test =", simple_current_FN,
                 "; num_of_future_FN_simple_test =", simple_future_FN,
                 "; min_PT_span =", PT_span,
                 "years; max_PT_span =", PT_span_max,
                 "years; min_CF_span =", CF_span, 
                 "years; min_PT_angle =", PT_angle,
                 "pi; min_CF_angle =", CF_angle,
                 "pi; drop_out_span =", drop_out_span, "years")

###############Merge datasets for data illustartion##############
PHGDH$score_type <- "TPM"
DRS$score_type <- "DRS"
overall <- bind_rows(PHGDH,DRS) %>% select(-REGTRYID)

## Annotate variables for facet graphing
PT_highlight$score_type <- "TPM"
CF_current_highlight$score_type <- "DRS"
CF_future_highlight$score_type <- "DRS"
test_record$score_type <- "TPM"
annotation$score_type <- "DRS"
decline_categorization$score_type <- "DRS"
num_simple_test$score_type <- "DRS"

#################Start draw graph###############
donor_list <- sequence$donor_id_alias
donor_plots <- list()
## draw multi-facet graph for each individual and
## store in a list 'donor_plots'
for (donor in donor_list){
  donor_plots[[donor]] <- 
    ggplot(overall %>% dplyr::filter(donor_id_alias==donor), 
           aes(time,value))+
    scale_colour_manual(values=c("#E41A1C", "#808B96"))+
    scale_fill_manual(values=c("#E41A1C", "#808B96"))+
    {if(nrow(PT_highlight %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_line(data= PT_highlight %>% dplyr::filter(donor_id_alias==donor), 
                col = "grey30", alpha = 0.4, size = 5)}+
    {if(nrow(CF_current_highlight %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_line(data= CF_current_highlight %>% 
                  dplyr::filter(donor_id_alias==donor), 
                aes(group = index), col = "grey30", 
                alpha = 0.4, size = 5)}+
    {if(nrow(CF_future_highlight %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_line(data= CF_future_highlight %>% 
                  dplyr::filter(donor_id_alias==donor), 
                aes(group = index), col = "#101097", 
                alpha = 0.4, size = 5)}+
    geom_line() +
    geom_point(size = 0.7)+
    {if(nrow(test_record %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_curve(data = test_record %>% 
                   dplyr::filter(donor_id_alias==donor), 
                 aes(x = x1, y = y1, xend = x2, yend = y2, group = final),
                 color = ifelse(
                   test_record[which(test_record$donor_id_alias == donor),10] ==
                     "wrong","#808B96","#E41A1C"),alpha = 0.8, curvature = -1)}+
    {if (nrow(drop_out_record %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_rect(data=drop_out_record %>% dplyr::filter(donor_id_alias==donor), 
                inherit.aes=FALSE, aes(xmin = start, xmax = end, 
                                       ymin = -Inf, ymax = Inf), 
                color="transparent", fill="#FFFF33", alpha=0.3)}+
    
    ## annotation
    geom_vline(data = convert_dates %>% dplyr::filter(donor_id_alias==donor), 
               aes(xintercept = convert_date), 
               linetype="dotted", color = "#E41A1C", size=0.5)+
    {if(nrow(annotation %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_label(data = annotation %>% dplyr::filter(donor_id_alias==donor), 
                 aes(x = -21, y = 120, label = label),
                 fill = ifelse(
                   annotation[which(annotation$donor_id_alias == donor),4] == 
                     "wrong","#808B96","#E41A1C"),
                 color = "white", fontface = "bold",label.size = NA)}+
    {if(nrow(num_simple_test %>% dplyr::filter(donor_id_alias==donor))!=0)
      geom_label(aes(x = -25, y = 120, label = num), 
                 data = num_simple_test %>% 
                   dplyr::filter(donor_id_alias==donor), 
                 fill = "#808B96", color = "white",label.size = NA,
                 fontface = "bold")}+
    {if(nrow(decline_categorization %>% 
             dplyr::filter(donor_id_alias==donor))!=0)
      geom_label(data = decline_categorization %>% 
                   dplyr::filter(donor_id_alias==donor), 
                 aes(x = -21, y = 100, label = label),
                 fill = "transparent",color = "black",
                 label.size = NA,size = 3)}+
    facet_wrap( ~ score_type, ncol = 1, 
                strip.position = "left", scales = "free_y") +
    
    ## theme
    theme_bw()+
    ggtitle(donor)+
    theme(plot.title = element_text(hjust = -0.1, size = 11,face = "bold"),
      legend.position = "none", 
      aspect.ratio=1/3, 
      axis.text.x = element_text(size = 6.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(2,5,2,5, "points"))+
    
    ## scale
    scale_x_continuous(limits=c(-27,0), n.breaks=14)+
    coord_panel_ranges(panel_ranges = list(
      list(x=c(-27,0.5), y=c(80,150)), # Panel 1
      list(x=c(-27,0.5), y=c(-1,70)) # Panel 2
    ))
}

## Combine individual plots
all <- ggarrange(plotlist=donor_plots,ncol = 5, nrow = 4)

## Add x-axis title
all <- annotate_figure(all,
                       bottom = text_grob("Time to drop out the study (years)",
                                          hjust = 0.5))

## Add caption
annotate_figure(all,bottom = textbox_grob(
  caption,
  hjust = 0.2, vjust = 3, halign = 0, valign = 0,
  width = unit(15, "inch"), height = unit(1.5, "inch"),
  orientation = "upright",
  margin = unit(c(5, 5, 5, 5), "pt")
))

