# Libraries needed
library(MASS)
library(tidyverse)
library(randomForest)
library(caret)
library(nnet)
library(pROC)
#install.packages('e1071', dependencies=TRUE)
set.seed(67)

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
nrow(rpm)*.5 #170
nrow(rpm)*.3 #102
nrow(rpm)*.2 #67

rpm <- rpm %>% 
  mutate(split = sample(rep(c("train", "valid", "test"), times = c(170, 102, 67))))

rpm_train_pre_over <- rpm %>% filter(split == "train") %>% select(-split) 
rpm_valid <- rpm %>% filter(split == "valid") %>% select(-split)
rpm_test  <- rpm %>% filter(split == "test") %>% select(-split)
#rpm  <- rpm %>% select(-split)

# Do OVERSAMPLING, duplicate all samples from CMS 1, 3 and 4 from the training set
rpm_train_CMS134 <- rpm_train_pre_over %>% filter(CMS == "CMS1" | CMS == "CMS3" | CMS == "CMS4")
rpm_train <- bind_rows(rpm_train_pre_over, rpm_train_CMS134)

# Check class imbalance
rpm_train %>% ggplot(aes(x = CMS, fill = CMS)) + 
  geom_bar() + 
  theme_minimal() +
  ggtitle("CMS counts")

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
####################### Make the LDA model ###################### 
# LDA cannot use redundant variables
# Therefore we use the no0 dfs
model_LDA <- lda(CMS ~ . - Sample.ID, data = rpm_train_no0)
model_LDA

# Predict CMS from the validation data
pred_rpm_valid_LDA <- predict(model_LDA, newdata = rpm_valid_no0)
pred_rpm_valid_LDA

# Print the confusion matrix
cmat_LDA <- confusionMatrix(pred_rpm_valid_LDA$class, rpm_valid_no0$CMS)
cmat_LDA

##################################################################################################################################      
##################################################################################################################################      
##################################################################################################################################      
####################### Make the RF model ####################### 
model_RF <- randomForest(CMS ~ . - Sample.ID, data = rpm_train)
model_RF

# Predict CMS from the validation data
pred_rpm_valid_RF <- predict(model_RF, newdata = rpm_valid)

# Print the confusion matrix
cmat_RF <- confusionMatrix(pred_rpm_valid_RF, rpm_valid$CMS)
cmat_RF

####################### Make the RF_no0 model #######################
model_RF_no0 <- randomForest(CMS ~ . - Sample.ID, data = rpm_train_no0)
model_RF_no0

# Predict CMS from the validation data
pred_rpm_valid_RF_no0 <- predict(model_RF_no0, newdata = rpm_valid_no0)

# Print the confusion matrix
cmat_RF_no0 <- confusionMatrix(pred_rpm_valid_RF_no0, rpm_valid_no0$CMS)
cmat_RF_no0

##################################################################################################################################      
####################### Make the LR model ####################### 

rpm_train2 <- rpm_train
rpm_train2$CMS_relevel <- relevel(rpm_train$CMS, ref = "CMS2") 

#model_LR <- multinom(CMS_relevel ~ . - CMS - Sample.ID, data = rpm_train2, MaxNWts = 8000, maxit = 180)
model_LR <- multinom(CMS ~ . - Sample.ID, data = rpm_train, MaxNWts = 8000, maxit = 120)
#model_LR

# Predict CMS from the validation data
pred_rpm_valid_LR <- predict(model_LR, newdata = rpm_valid)

# Print the confusion matrix
#cmat_LR <- confusionMatrix(pred_rpm_valid_LR, relevel(rpm_valid$CMS, ref = "CMS2"))
cmat_LR <- confusionMatrix(pred_rpm_valid_LR, rpm_valid$CMS)
cmat_LR

####################### Make the LR_no0 model ####################### 
rpm_train_no0_2 <- rpm_train_no0
rpm_train_no0_2$CMS_relevel <- relevel(rpm_train_no0$CMS, ref = "CMS2") 

#model_LR_no0 <- multinom(CMS_relevel ~ . - CMS - Sample.ID, data = rpm_train_no0_2, MaxNWts = 8000, maxit = 180)
model_LR_no0 <- multinom(CMS ~ . - Sample.ID, data = rpm_train_no0, MaxNWts = 8000, maxit = 120)
#model_LR_no0

# Predict CMS from the validation data
pred_rpm_valid_LR_no0 <- predict(model_LR_no0, newdata = rpm_valid_no0)

# Print the confusion matrix
#cmat_LR <- confusionMatrix(pred_rpm_valid_LR, relevel(rpm_valid$CMS, ref = "CMS2"))
cmat_LR_no0 <- confusionMatrix(pred_rpm_valid_LR_no0, rpm_valid_no0$CMS)
cmat_LR_no0

# ROC metric (DOES NOT WORK YET)
pred_rpm_valid_RF_no0_prob <- predict(model_RF_no0, newdata = rpm_valid_no0, type = "prob")
str(pred_rpm_valid_RF_no0_prob)
multiclass.roc(rpm_valid_no0$CMS, pred_rpm_valid_RF_no0_prob)
multiclass.roc(rpm_valid_no0$CMS, pred_rpm_valid_RF_no0_prob)
length(rpm_valid_no0$CMS) == length(as.numeric(pred_rpm_valid_RF_no0_prob))

##################################################################################################################################      
####################### Compare models ####################### 
cmat_LDA$table
cmat_LDA$overall[1]
cmat_RF$table
cmat_RF$overall[1]
cmat_RF_no0$table
cmat_RF_no0$overall[1]
cmat_LR$table
cmat_LR$overall[1]
cmat_LR_no0$table
cmat_LR_no0$overall[1]

##################################################################################################################################      
####################### Check importance of features ####################### 

imp_RF <- importance(model_RF)
imp_RF_df <- tibble(variable = rownames(imp_RF), importance = c(imp_RF)) %>% arrange(-importance) 

imp_RF_df %>%
  ggplot(aes(x = reorder(variable, -importance), y = importance, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")





