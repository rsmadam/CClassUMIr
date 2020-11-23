# Libraries needed
library(MASS)
library(tidyverse)
library(randomForest)
library(caret)
library(nnet)
library(pROC)
#install.packages('e1071', dependencies=TRUE)
set.seed(67)

##################################################################################################################################      
# Open the csv file
rpm <- read.csv("data/cleaned/RPM_prob_p005.csv")
rc <- read.csv("data/cleaned/RC_prob_p005.csv")

# Make all miRNA data numeric
for (i in 3:ncol(rpm)){
  rpm[, i] <- as.numeric(rpm[, i])
}

# Only use the features with read count mean > 1
col_means <- colMeans(rc[-c(1,2)])
feature_bigger1 <- min(which(sort(col_means)>1)) # First feature with a mean > 1
feature_names_bigger1 <- names(sort(col_means)[feature_bigger1 :length(sort(col_means))]) 
feature_names <- c("Sample.ID", "CMS", feature_names_bigger1)

rpm <- rpm %>% select(feature_names)


# Divide into train validation and test sets
#nrow(rpm)*.5 #170
#nrow(rpm)*.3 #102
#nrow(rpm)*.2 #67

##################################################################################################################################      
####################### Make a for loop for cross validation with oversampling ####################### 
# Make an empty variable
cv_comparison = NULL
for (i in 1:10){
  
  rpm_split <- rpm %>% 
    mutate(split = sample(rep(c("train", "valid", "test"), times = c(170, 102, 67))))

  rpm_train_pre_over <- rpm_split %>% filter(split == "train") %>% select(-split) 
  rpm_valid <- rpm_split %>% filter(split == "valid") %>% select(-split)
  rpm_test  <- rpm_split %>% filter(split == "test") %>% select(-split)
  #rpm  <- rpm %>% select(-split)
  
  # Do OVERSAMPLING, duplicate all samples from CMS 1, 3 and 4 from the training set
  rpm_train_CMS134 <- rpm_train_pre_over %>% filter(CMS == "CMS1" | CMS == "CMS3" | CMS == "CMS4")
  rpm_train <- bind_rows(rpm_train_pre_over, rpm_train_CMS134)

  # Check class imbalance
  #rpm_train %>% ggplot(aes(x = CMS, fill = CMS)) + 
  #  geom_bar() + 
  #  theme_minimal() +
  #  ggtitle("CMS counts")
  
  # Make another df to get rid of redundant variables (those that are 0 for all samples). 
  # For example needed for LDA
  remove_0_variables <- function(df, df_2, df_3){
    df_no0 <- df[c(1, 2)]
    for (i in 3:(ncol(df)-1)){ #ncol(df)-1 because of the added split column
      if ((sum(df[i]) != 0) & (sum(df_2[i]) != 0) & (sum(df_3[i]) != 0)){
        df_no0 <- bind_cols(df_no0, df[i]) 
      }
    }
    df_no0 <- bind_cols(df_no0, df[ncol(df)]) 
    return(df_no0)
  } 

  # Use this function to get rid of variables that are equal for all samples, within one of the 3 groups.
  rpm_train_no0 <- remove_0_variables(rpm_train, rpm_valid, rpm_test)
  rpm_valid_no0 <- remove_0_variables(rpm_valid, rpm_train, rpm_test)
  rpm_test_no0 <- remove_0_variables(rpm_test, rpm_train, rpm_valid)

  ##################################################################################################################################      
  ##################################################################################################################################      
  ####################### Make the RF model ####################### 
  model_RF <- randomForest(CMS ~ . - Sample.ID, data = rpm_train)
  #model_RF
  
  # Predict CMS from the validation data
  pred_rpm_valid_RF <- predict(model_RF, newdata = rpm_valid)
  
  # Print the confusion matrix
  cmat_RF <- confusionMatrix(pred_rpm_valid_RF, rpm_valid$CMS)
  #cmat_RF
  
  ####################### Make the RF_no0 model #######################
  model_RF_no0 <- randomForest(CMS ~ . - Sample.ID, data = rpm_train_no0)
  #model_RF_no0
  
  # Predict CMS from the validation data
  pred_rpm_valid_RF_no0 <- predict(model_RF_no0, newdata = rpm_valid_no0)
  
  # Print the confusion matrix
  cmat_RF_no0 <- confusionMatrix(pred_rpm_valid_RF_no0, rpm_valid_no0$CMS)
  #cmat_RF_no0
  
  # Make a row with the accuracy of the models, how well CMS 3 is predicted, 
  # and the amount of CMS3 in the training set (to check whether that matters, doesnt seem so)
  new_row <- data_frame("model_RF acc" = cmat_RF$overall[1], 
                        "model_RF_no0 acc" = cmat_RF_no0$overall[1], 
                        "specifiticy CMS3" = cmat_RF$byClass[3, 1], 
                        "prevalence CMS3" =  cmat_RF$byClass[3, 8])
  
  # Add this row to cv_comparison
  cv_comparison <- rbind(cv_comparison, new_row)
  ##################################################################################################################################      
}
cv_comparison 

cv_comparison[order(cv_comparison$`model_RF acc`, decreasing = TRUE),]
colMeans(cv_comparison)

# ROC metric (DOES NOT WORK YET)
#pred_rpm_valid_RF_no0_prob <- predict(model_RF_no0, newdata = rpm_valid_no0, type = "prob")
#str(pred_rpm_valid_RF_no0_prob)
#multiclass.roc(rpm_valid_no0$CMS, pred_rpm_valid_RF_no0_prob)
#multiclass.roc(rpm_valid_no0$CMS, pred_rpm_valid_RF_no0_prob)
#length(rpm_valid_no0$CMS) == length(as.numeric(pred_rpm_valid_RF_no0_prob))

##################################################################################################################################      
####################### Check importance of features ####################### 
imp_RF <- importance(model_RF)
imp_RF_df <- tibble(variable = rownames(imp_RF), importance = c(imp_RF)) %>% arrange(-importance) 

imp_RF_df %>%
  ggplot(aes(x = reorder(variable, -importance), y = importance, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")





