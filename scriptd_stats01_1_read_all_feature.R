#!/usr/bin/env Rscript

# library to read xls files
#library(gdata)

data_dir = "data_d01_features"
slash = .Platform$file.sep

# read subject information files:
subject.info = read.table(paste(data_dir, "pnc_subject_info.csv", sep = slash), header=TRUE, sep = ',')

print(head(subject.info))
print(dim(subject.info))
#print(table(subject.info$ptsd))
#print(table(subject.info[subject.info$ptsd==0,]$Sex))
#print(table(subject.info[subject.info$ptsd==1,]$Sex))
#print(table(subject.info[subject.info$ptsd==2,]$Sex))

#print("mean age of 3 ptsd groups:")
#ag <- aggregate(. ~ ptsd, subject.info[, c("ptsd", "age")], function(x) c(mean = mean(x), sd = sd(x)))
#print(ag)

# recode gender:
subject.info$Gender = NaN;
subject.info[subject.info$Sex=='F',]$Gender = 0;
subject.info[subject.info$Sex=='M',]$Gender = 1;
print(head(subject.info))
print(table(subject.info$Sex))

print("mean age of males and females:")
ag <- aggregate(. ~ Sex, subject.info[, c("Sex", "age_at_cnb")], function(x) c(mean = mean(x), sd = sd(x)))
print(ag)

# -------------------------- cat12 gray matter volume ---------------------

#cat.vbm.hammers = read.table(paste(data_dir, "ROI_catROI_hammers_Vgm.csv", sep = slash), header=TRUE, sep = ",")[,-1]

# -------------------------- cat12 gray matter volume ---------------------

cat.vbm.neuromorph = read.table(paste(data_dir, "ROI_catROI_neuromorphometrics_Vgm.csv", sep = slash), header=TRUE, sep = ",")
cat.vbm.subid = data.frame(SUBJID = cat.vbm.neuromorph$names)

# remove columns with more than 70% zeros
percent.zeros = lapply(cat.vbm.neuromorph, function(x){ length(which(x==0))/length(x)})
#print(percent.zeros)
cat.vbm.neuromorph = cat.vbm.neuromorph[, percent.zeros<.3]
cat.vbm.neuromorph = cat.vbm.neuromorph[, -1]
print(head(cat.vbm.neuromorph))


# -------------------------- fsl gray matter volume ---------------------

#fsl.vbm = read.table(paste(data_dir, "outd01_all_fslvbm_features.txt", sep = slash), header=TRUE)

# ------------------------- spm gray matter volume ----------------------

#spm.vbm = read.table(paste(data_dir, "outd01_all_spmvbm_features.txt", sep = slash), header=TRUE)
#spm.vbm.aal = read.table(paste(data_dir, "outd01_all_spm.vbm.aal_features.txt", sep = slash), header=TRUE)
#spm.vbm.bn246 = read.table(paste(data_dir, "outd01_all_spm.vbm.bn246_features.txt", sep = slash), header=TRUE, )

# ------------------------------- alff ---------------------------------

#alff = read.table(paste(data_dir, "outd01_all_alff_features.txt", sep = slash), header=TRUE)
#alff.aal = read.table(paste(data_dir, "outd01_all_alff.aal_features.txt", sep = slash), header=TRUE, )
alff.bn246 = read.table(paste(data_dir, "ALFF_AvgExtract.txt", sep = slash), header=TRUE, sep = '\t')
print(head(alff.bn246[,c(1,ncol(alff.bn246))]))
alff.bn246 = alff.bn246[,-c(1,ncol(alff.bn246))]
print(dim(alff.bn246))

# ------------------------------- falff --------------------------------

#falff = read.table(paste(data_dir, "outd01_all_falff_features.txt", sep = slash), header=TRUE)
#falff.aal = read.table(paste(data_dir, "outd01_all_falff.aal_features.txt", sep = slash), header=TRUE, )
falff.bn246 = read.table(paste(data_dir, "fALFF_AvgExtract.txt", sep = slash), header=TRUE, sep = '\t')
print(head(falff.bn246[,c(1,ncol(alff.bn246))]))
falff.bn246 = falff.bn246[,-c(1,ncol(falff.bn246))]
print(dim(falff.bn246))
# --------------------------------- reho ----------------------------------

#reho = read.table(paste(data_dir, "outd01_all_reho_features.txt", sep = slash), header=TRUE)
#reho.aal = read.table(paste(data_dir, "outd01_all_reho.aal_features.txt", sep = slash), header=TRUE, )
reho.bn246 = read.table(paste(data_dir, "Reho_AvgExtract.txt", sep = slash), header=TRUE, sep = '\t')
print(head(reho.bn246[,c(1,ncol(reho.bn246))]))
reho.bn246 = reho.bn246[,-c(1,ncol(reho.bn246))]
print(dim(reho.bn246))
# ------------------------------- label.fa --------------------------------

label.fa = read.table(paste(data_dir, "WMlabelResults_FA_all.csv", sep = slash), header = TRUE, sep = ',')
print(head(label.fa[, c(1,2,ncol(label.fa))]))
label.fa = label.fa[, -c(1,2, ncol(label.fa))]
print(dim(label.fa))

# ------------------------------- label.md --------------------------------
label.md = read.table(paste(data_dir, "WMlabelResults_MD_all.csv", sep = slash), header = TRUE, sep = ',')
print(head(label.md[, c(1,2,ncol(label.md))]))
label.md = label.md[, -c(1,2, ncol(label.md))]
print(dim(label.md))

# ------------------------------- tract.fa --------------------------------

tract.fa = read.table(paste(data_dir, "WMtractResults_FA_all.csv", sep = slash), header = TRUE, sep = ',')
print(head(tract.fa[, c(1,2,ncol(tract.fa))]))
tract.fa = tract.fa[, -c(1,2, ncol(tract.fa))]
print(dim(tract.fa))

# ------------------------------- tract.md --------------------------------

tract.md = read.table(paste(data_dir, "WMtractResults_MD_all.csv", sep = slash), header = TRUE, sep = ',')
print(head(tract.md[, c(1,2,ncol(tract.md))]))
tract.md = tract.md[, -c(1,2, ncol(tract.md))]
print(dim(tract.md))

