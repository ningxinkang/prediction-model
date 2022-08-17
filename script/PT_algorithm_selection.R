##__________________________________________________________________
## Author: Ningxin Kang (nik010@ucsd.edu)       
## Last update: 2022-08-17    
## File: PT_algorithm_selection.R          
## File Summary:
##  This file include the code that allow user to compare
##  between three versions of PT algorithm using the candidate
##  parameter combinations.
## Note:
##  parameter: 24 * 7 candidate combinations
##  CF: conditional
##__________________________________________________________________

## Change according to need: PT_version1/2/3
PT_version <- "PT_version3"
source_dir <- paste("./script/algorithm/", PT_version, ".R", sep = "")
output_dir <- 
  paste("./result/algorithm_selection/", PT_version, ".csv", sep = "")

############Loading Libraries############
library(dplyr)
library(tidyr)

############Set up loop range#############
PT_span_list <- seq(0.5, 2, by = 0.5)
PT_span_max <- 4.5
CF_span_list <- seq(0.5, 1.5, by = 0.5)
PT_angle_list <- seq(0.4, 0.45, by = 0.05)
CF_angle <- -0.45
drop_out_span_list <- seq(0, 3, by = 0.5)
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
  read.csv("./data/ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  drop_na(TPM) %>% arrange(x.axis) %>% group_by(donor_id_alias)
donors_PHGDH = split(PHGDH$TPM,PHGDH$donor_id_alias)
donors_PHGDH_t = split(PHGDH$x.axis,PHGDH$donor_id_alias)
PHGDH <- PHGDH[,1:5]
colnames(PHGDH)[3] <- "value"
colnames(PHGDH)[2] <- "time"

DRS <- 
  read.csv("./data/ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  drop_na(DRS) %>% arrange(x.axis) %>% group_by(donor_id_alias)
donors_DRS = split(DRS$DRS,DRS$donor_id_alias)
donors_DRS_t = split(DRS$x.axis,DRS$donor_id_alias)
DRS <- DRS[,c(1,2,9,5,4)]
colnames(DRS)[3] <- "value"
colnames(DRS)[2] <- "time"

sequence <- distinct(PHGDH[,4:5]) %>% arrange(donor_id_alias)
convert_dates <- 
  read.csv("./data/ADRC_PHGDH_with_cogscore_noAD.csv", header = TRUE, sep = ',') %>% 
  select(REGTRYID,donor_id_alias,convert_date) %>% 
  arrange(donor_id_alias) %>% 
  distinct()


############Sourcing function############
source(source_dir)
source('./script/algorithm/CF.R')
source('./script/algorithm/prediction_model_conditionalCF.R')

############Set up output table#############
contigency_table <- data.frame(matrix(ncol = 18, nrow = 0))
z <- c("accuracy_per_full_test", "PPV_per_full_test", "NPV_per_full_test", 
       "accuracy_per_simple_test", "PPV_per_simple_test", "NPV_per_simple_test",
       "total_simple_test", "avg_simple_test_per_full_test",
       "positive_simple","negative_simple","TP_simple","TN_simple",
       "PT_span", "PT_max_span","CF_span", 
       "PT_angle", "CF_angle", "drop_out_span")
colnames(contigency_table) <- z

for (i in 1:length(PT_span_list)){
  for (j in 1:length(CF_span_list)){
    for (k in 1:length(PT_angle_list)){
      for (l in 1:length(drop_out_span_list)){
        source('./script/storing_variables.R')
        info <- prediction_model(PT_span_list[[i]], 
                                 PT_span_max, 
                                 CF_span_list[[j]],
                                 PT_angle_list[[k]], 
                                 CF_angle, 
                                 drop_out_span_list[[l]], 
                                 donors_PHGDH, 
                                 donors_PHGDH_t, 
                                 donors_DRS, 
                                 donors_DRS_t, 
                                 sequence)
        contigency_table <<- contigency_table %>% 
          rbind(data.frame(accuracy_per_full_test= info[[1]], 
                           PPV_per_full_test = info[[2]], 
                           NPV_per_full_test = info[[3]], 
                           accuracy_per_simple_test= info[[4]], 
                           PPV_per_simple_test = info[[5]], 
                           NPV_per_simple_test = info[[6]], 
                           total_simple_test = info[[7]], 
                           avg_simple_test_per_full_test =info[[7]]/20, 
                           positive_simple = info[[8]], 
                           negative_simple = info[[9]], 
                           TP_simple = info[[10]], 
                           TN_simple= info[[11]], 
                           PT_span=PT_span_list[[i]], 
                           PT_max_span = PT_span_max, 
                           CF_span = CF_span_list[[j]], 
                           PT_angle = PT_angle_list[[k]], 
                           CF_angle = CF_angle, 
                           drop_out_span = drop_out_span_list[[l]]))
      }
    }
  }
}
write.table(contigency_table, output_dir,
            row.names=FALSE, col.names=TRUE, sep=",", quote = FALSE)