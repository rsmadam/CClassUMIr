library(tidyverse)

# Open the csv file
rpm <- read.csv("data/cleaned/RPM_prob_p005.csv")
rc <- read.csv("data/cleaned/RC_prob_p005.csv")

# Make count data numeric
for (i in 3:ncol(rpm)){
  rpm[, i] <- as.numeric(rpm[, i])
} 

rpm %>% ggplot(aes(x = CMS, fill = CMS)) + 
  geom_bar() + 
  theme_minimal() +
  ggtitle("CMS counts")
#ggsave(filename = "CMS_counts.png", path = "Data/Images/")

# PCA
pca <- prcomp(rpm[, -c(1, 2)]) 

# Plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

# Make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per[1:5], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

# now make a fancy looking plot that shows the PCs and the variation:
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], CMS = rpm$CMS)
head(pca_data)

pca_data %>% ggplot(aes(x = pca$x[,1], y = pca$x[,2], color = CMS)) +
  #geom_point(alpha = 0.7) +
  geom_density_2d() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA") +
  theme_minimal()
#ggsave(filename = "PCA.png", path = "Data/Images/")
