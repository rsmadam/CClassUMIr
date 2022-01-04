# CMS-miRaCl
Training a microRNA-assigned CMS-classifier (miRaCl).
Scripts for gathering and cleaning data for training and validation sets, for training a random forest classifier, and for reproducing figures.
Compiled classifiers are in classifiers.RData, miRaCl is the full version, miRaCl20 is the parsimonous version, miRaCl20A is the microarray version. 
You can load the classifiers into RStudio after download using 
load("your-download-location/classifiers.RData").