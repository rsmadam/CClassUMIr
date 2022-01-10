# CMS-miRaCl
Training a microRNA-assigned CMS-classifier (miRaCl).

Scripts for gathering and cleaning data for training and validation sets, for training a random forest classifier, and for reproducing figures.

Compiled classifiers are in classifiers.RData, `miRaCl` is the full version, `miRaCl20` is the parsimonous version, `miRaCl20A` is the microarray version. 
You can load the classifiers into RStudio after download using: 
```{r}
require(caret)
load("your-download-location/classifiers.RData")
```

You could use the classifier as such: 

```{r}
### prepare your dataset using variance stabilizing transformation ### 
# (after potentially removing outliers and summarizing miRNA isoforms, depending on your dataset) #
require(DESeq2)
miR_DATASET_vst <- as.data.frame( t( varianceStabilizingTransformation(
  as.matrix( round( miR_DATASET, 0 ) ) ) ) )

### prediction of posterior probability per class using miRaCl-20 ###
pred_prob_cms <- predict(miRaCl20, newdata = miR_DATASET_vst, type = "prob")

### method to get an estimate of the confidence ###
diffSecond <- function(x, output) {
  ordered = sort(x, decreasing=T)
  confid = abs(ordered[1] - ordered[2])
  confid
} #gets the absolute difference between first and second class prob.
pred_read_RF$d2ndProb <- apply(pred_prob_cms, 1, diffSecond)

### prediction of most likely class label ###
pred_prob_cms$CMS.miRaCl20 <- predict(miRaCl20, newdata = miR_DATASET_vst, type = "raw")
summary(pred_prob_cms$CMS.miRaCl20)

```        
                          