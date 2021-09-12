#!/usr/bin/env Rscript

# library to read xls files
#library(gdata)

# read subject information files:
subject.info = read.table("outd01_subject_info.txt", header=TRUE)

print(head(subject.info))
print(dim(subject.info))
print(table(subject.info$ptsd))
print(table(subject.info[subject.info$ptsd==0,]$Sex))
print(table(subject.info[subject.info$ptsd==1,]$Sex))
print(table(subject.info[subject.info$ptsd==2,]$Sex))
# -------------------------- fsl gray matter volume ---------------------

fc = read.table("outd01_fc_features.txt", header=TRUE)

