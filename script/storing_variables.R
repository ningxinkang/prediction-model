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