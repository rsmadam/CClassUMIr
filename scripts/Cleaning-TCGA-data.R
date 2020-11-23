# Libraries needed
library(tidyverse)

# Import data
miRNA_samples <- as_tibble(read.table(file = 'data/raw/gdc_sample_sheet.2019-07-01.tsv', sep = '\t', header = TRUE))
CMS_samples <- as_tibble(read.table(file = 'data/raw/cms_labels.txt', sep = '\t', header = TRUE))
# Adjust SampleID so it is the same as the Case.ID
CMS_samples$SampleId <- gsub("\\.", "-", CMS_samples$SampleId) 

# There can be multiple miRNA expression files per Sample.ID
strange_double <- filter(miRNA_samples, Sample.ID == "TCGA-A6-6780-01A")

# For loop to count the amount of overlap
count <- 0
for (i in unique(miRNA_samples$Case.ID)){
  if (i %in% unique(CMS_samples$SampleId)){
    count <- count + 1
  }
}
print(count)

#Combine datasets
inner_miRNA_CMS <- merge(filter(miRNA_samples, miRNA_samples$Sample.Type != "Solid Tissue Normal"), CMS_samples, by.x = "Case.ID", by.y = "SampleId")

# List of the File.IDs
File.ID_list <- list.files("data/raw/miRNA_TCGA/")

# Check whether the File.IDs are indeed the same
count2 <- 0
for (i in inner_miRNA_CMS$File.ID){
  if (i %in% File.ID_list){
    count2 <- count2 + 1
  }
}
print(count2)

# For checking whether the values are correct
# File.ID
inner_miRNA_CMS$File.ID[3]
# File.Name
inner_miRNA_CMS$File.Name[3]
# Sample.ID
inner_miRNA_CMS$Sample.ID[3]

# Load an example to get the miRNA names (first column)
example <- as_tibble(read.table(file = 'data/raw/example.txt', sep = '\t', header = TRUE))

## RPM
# Loop over the files and add the reads per milion (rpm) miRNA mapped columns
rpm_table <- example[1]
for (i in 1:length(inner_miRNA_CMS$Case.ID)){
  # Get the location of the file
  file <- paste0('data/raw/miRNA_TCGA/', inner_miRNA_CMS$File.ID[i], "/", inner_miRNA_CMS$File.Name[i])
  # Open the file and add the right column to to df
  rpm_table <- add_column(rpm_table, as_tibble(read.table(file = file, sep = '\t', header = TRUE))$reads_per_million_miRNA_mapped)
  # Add the column name, which is the Sample.ID
  colnames(rpm_table)[ncol(rpm_table)] <- toString(inner_miRNA_CMS$Sample.ID[i])
}

## RC
# Loop over the files and add the read counts (rc) columns
rc_table <- example[1]
for (i in 1:length(inner_miRNA_CMS$Case.ID)){
  file <- paste0('data/raw/miRNA_TCGA/', inner_miRNA_CMS$File.ID[i], "/", inner_miRNA_CMS$File.Name[i])
  rc_table <- add_column(rc_table, as_tibble(read.table(file = file, sep = '\t', header = TRUE))$read_count)
  colnames(rc_table)[ncol(rc_table)] <- toString(inner_miRNA_CMS$Sample.ID[i])
}


# Function to transpose a dataframe
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(df = .) %>%
    tibble::as_data_frame(x = .)
  return(t_df)
}

## RPM
# Transpose the rpm_table
rpm_t <- transpose_df(rpm_table)
# Set the headers right
colnames(rpm_t) <- c("Sample.ID", t(rpm_t[1, -1]))
rpm_t <- rpm_t[-1, ]
# Add the CMS lables, with "CMS" colname
rpm <- add_column(rpm_t, inner_miRNA_CMS$CMS, .after = 1)
colnames(rpm)[2] <- "CMS"
# Save the cleaned RPM data
write_csv(rpm, "data/cleaned/RPM.csv")

## RC
# Transpose the rc_table and set the headers right
rc_t <- transpose_df(rc_table)
colnames(rc_t) <- c("Sample.ID", t(rc_t[1, -1]))
rc_t <- rc_t[-1, ]
# Add the CMS lables, with "CMS" colname
rc <- add_column(rc_t, inner_miRNA_CMS$CMS, .after = 1)
colnames(rc)[2] <- "CMS"
# Save the cleaned RC data
write_csv(rc, "data/cleaned/RC.csv")

