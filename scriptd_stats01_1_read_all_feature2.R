#!/usr/bin/env Rscript

# library to read xls files
#library(gdata)

data_dir = "data_d01_features"
slash = .Platform$file.sep

# read subject information files:
subject.info = read.table(paste(data_dir, "outd01_all_subject_info.txt", sep = slash), header=TRUE)


# udpates 2/1/2021/ add trauma_age as features:
# read subject information files:
subject.info = read.table("data_d01_features/outd01_all_subject_info2.txt", sep = ' ', 
                          header=TRUE, stringsAsFactors = F)
colnames(subject.info) = c('SUBJID','Sex', 'age',  'trauma_age', 'ptsd')
# code sex: 1 as female, 2 as male
subject.info[subject.info$Sex=='F',]$Sex=1
subject.info[subject.info$Sex=='M',]$Sex=2

# to avoid problem when converting df to matrix:
subject.info <- type.convert(subject.info, as.is = TRUE)




print(head(subject.info))
print(dim(subject.info))
print(table(subject.info$ptsd))
print(table(subject.info[subject.info$ptsd==0,]$Sex))
print(table(subject.info[subject.info$ptsd==1,]$Sex))
print(table(subject.info[subject.info$ptsd==2,]$Sex))

print("mean age of 3 ptsd groups:")
ag <- aggregate(. ~ ptsd, subject.info[, c("ptsd", "age")], function(x) c(mean = mean(x), sd = sd(x)))
print(ag)
# -------------------------- cat12 gray matter volume ---------------------

cat.vbm.hammers = read.table(paste(data_dir, "ROI_catROI_hammers_Vgm.csv", sep = slash), header=TRUE, sep = ",")[,-1]

# -------------------------- cat12 gray matter volume ---------------------

cat.vbm.neuromorph = read.table(paste(data_dir, "ROI_catROI_neuromorphometrics_Vgm.csv", sep = slash), header=TRUE, sep = ",")[,-1]

# -------------------------- fsl gray matter volume ---------------------

fsl.vbm = read.table(paste(data_dir, "outd01_all_fslvbm_features.txt", sep = slash), header=TRUE)

# ------------------------- spm gray matter volume ----------------------

spm.vbm = read.table(paste(data_dir, "outd01_all_spmvbm_features.txt", sep = slash), header=TRUE)
spm.vbm.aal = read.table(paste(data_dir, "outd01_all_spm.vbm.aal_features.txt", sep = slash), header=TRUE)
spm.vbm.bn246 = read.table(paste(data_dir, "outd01_all_spm.vbm.bn246_features.txt", sep = slash), header=TRUE, )

# ------------------------------- alff ---------------------------------

alff = read.table(paste(data_dir, "outd01_all_alff_features.txt", sep = slash), header=TRUE)
alff.aal = read.table(paste(data_dir, "outd01_all_alff.aal_features.txt", sep = slash), header=TRUE, )
alff.bn246 = read.table(paste(data_dir, "outd01_all_alff.bn246_features.txt", sep = slash), header=TRUE, )

# ------------------------------- falff --------------------------------

falff = read.table(paste(data_dir, "outd01_all_falff_features.txt", sep = slash), header=TRUE)
falff.aal = read.table(paste(data_dir, "outd01_all_falff.aal_features.txt", sep = slash), header=TRUE, )
falff.bn246 = read.table(paste(data_dir, "outd01_all_falff.bn246_features.txt", sep = slash), header=TRUE, )

# --------------------------------- reho ----------------------------------

reho = read.table(paste(data_dir, "outd01_all_reho_features.txt", sep = slash), header=TRUE)
reho.aal = read.table(paste(data_dir, "outd01_all_reho.aal_features.txt", sep = slash), header=TRUE, )
reho.bn246 = read.table(paste(data_dir, "outd01_all_reho.bn246_features.txt", sep = slash), header=TRUE, )

# ------------------------------- label.fa --------------------------------

label.fa = read.table(paste(data_dir, "outd01_all_label.fa_features.txt", sep = slash), header = TRUE)

# ------------------------------- label.md --------------------------------
label.md = read.table(paste(data_dir, "outd01_all_label.md_features.txt", sep = slash), header = TRUE)

# ------------------------------- tract.fa --------------------------------

tract.fa = read.table(paste(data_dir, "outd01_all_tract.fa_features.txt", sep = slash), header = TRUE)

# ------------------------------- tract.md --------------------------------

tract.md = read.table(paste(data_dir, "outd01_all_tract.md_features.txt", sep = slash), header = TRUE)

