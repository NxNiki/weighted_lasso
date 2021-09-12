
reproducibility.index = function(x, k){
  r =  (abs(x - k / 2) - (k / 2) %% 1) / floor(k / 2)
}


file.name.lasso = 'coefs_ptsd_1_vs_2_folds_10_5_alpha_1_none_logtransform_1_mf.vbm.aal.dti.csv'
file.name.elasticnet = 'coefs_ptsd_1_vs_2_folds_10_5_alpha_0.5_none_logtransform_1_mf.vbm.aal.dti.csv'

file.name.wlasso = 'coefs_ptsd_1_vs_2_folds_10_5_alpha_1_mean.diff.boot_logtransform_1_mf.vbm.aal.dti.csv'
file.name.welasticnet = 'coefs_ptsd_1_vs_2_folds_10_5_alpha_0.5_mean.diff.boot_logtransform_1_mf.vbm.aal.dti.csv'



#file_date = '2019-06-10'
#file_date = '2019-07-17'# reproducibility contains zero, nfold outer=10
#file_date = '2019-07-11'# reproducibility does not contain zero, nfold outer=10
#file_date = '2020-03-03' # udpated penalty calculation
file_date = '2021-02-26' # add trauma age as features.
nfolds = "10"

# elasticnet:
file_dir1 = paste0("result_d03_trauma_age_cv_glmnet_folds_", nfolds, "_5_alpha_0.5_none_logtransform_1_", file_date)

# alpha: 1 lasso, type.measure: class
file_dir2 = paste0("result_d03_trauma_age_cv_glmnet_folds_", nfolds, "_5_alpha_1_none_logtransform_1_", file_date)

# weighted elasticnet:
file_dir3 = paste0("result_d03_trauma_age_cv_glmnet_folds_", nfolds, "_5_alpha_0.5_mean.diff.boot_logtransform_1_", file_date)

# weighted alpha: 1 lasso, type.measure: class
file_dir4 = paste0("result_d03_trauma_age_cv_glmnet_folds_", nfolds, "_5_alpha_1_mean.diff.boot_logtransform_1_", file_date)


#  elastic net:
file.elasticnet = read.table(paste0(file_dir1, .Platform$file.sep, file.name.elasticnet), sep = ",", header = T)
#mydf$task <- factor(mydf$task, levels = c("up", "down", "left", "right", "front", "back"))

#  lasso:
file.lasso = read.table(paste0(file_dir2, .Platform$file.sep, file.name.lasso), sep = ",", header = T)

# weighted elastic net:
file.welasticnet = read.table(paste0(file_dir3, .Platform$file.sep, file.name.welasticnet), sep = ",", header = T)

# weighted lasso:
file.wlasso = read.table(paste0(file_dir4, .Platform$file.sep, file.name.wlasso), sep = ",", header = T)

# make ROI index for each modality:
roi.list=file.lasso[,1]
aal.length = length(grep('aal', roi.list))
fa.label.length = length(grep('label.fa', roi.list))
md.label.length = length(grep('label.md', roi.list))
fa.tract.length = length(grep('tract.fa', roi.list))
md.tract.length = length(grep('tract.md', roi.list))

roi.index=c(rep(NA, 4), seq(1, aal.length), seq(1, fa.label.length), seq(1, fa.tract.length), seq(1, md.label.length), seq(1, md.tract.length))


# find common rows with nonzero count >=7 for at least one method:
nonzerocount.elasticnet = apply(file.elasticnet[, 2:11] != 0, 1, sum)
nonzerocount.lasso = apply(file.lasso[, 2:11] != 0, 1, sum)
nonzerocount.weighted.elasticnet = apply(file.welasticnet[, 2:11] != 0, 1, sum)
nonzerocount.weighted.lasso = apply(file.wlasso[, 2:11] != 0, 1, sum)

common.rows=which(nonzerocount.elasticnet>=7|nonzerocount.lasso>=7|nonzerocount.weighted.elasticnet>=7|nonzerocount.weighted.lasso>=7)
common.rows=common.rows[common.rows>=2]

# compute mean coefs and select rows:
mean.coefs = apply(file.elasticnet[, 2:11], 1, mean)
file.elasticnet.count = data.frame(brain.region = file.elasticnet[common.rows, 1], 
                                   roi.index = roi.index[common.rows],
                                   elasticnet =  nonzerocount.elasticnet[common.rows],
                                   elasticnet.coefs = mean.coefs[common.rows])

nonzerocount = apply(file.lasso[, 2:11] != 0, 1, sum)
mean.coefs = apply(file.lasso[, 2:11], 1, mean)
file.lasso.count = data.frame(brain.region = file.lasso[common.rows, 1], 
                                   lasso =  nonzerocount.lasso[common.rows],
                                   lasso.coefs = mean.coefs[common.rows])

nonzerocount = apply(file.welasticnet[, 2:11] != 0, 1, sum)
mean.coefs = apply(file.welasticnet[, 2:11], 1, mean)
file.welasticnet.count = data.frame(brain.region = file.welasticnet[common.rows, 1], 
                                   weighted.elasticnet =  nonzerocount.weighted.elasticnet[common.rows],
                                   weighted.elasticnet.coefs = mean.coefs[common.rows])

nonzerocount = apply(file.wlasso[, 2:11] != 0, 1, sum)
mean.coefs = apply(file.wlasso[, 2:11], 1, mean)
file.wlasso.count = data.frame(brain.region = file.wlasso[common.rows, 1], 
                                   weighted.lasso =  nonzerocount.weighted.lasso[common.rows],
                                   weighted.lasso.coefs = mean.coefs[common.rows])

file.merge = Reduce(function(x,y) merge(x = x, y = y, by = "brain.region", all = T), 
                    list(file.elasticnet.count[,c(1,2,3)], file.welasticnet.count[,c(1,2)], file.lasso.count[,c(1,2)], file.wlasso.count[,c(1,2)]))

file.merge.coefs = Reduce(function(x,y) merge(x = x, y = y, by = "brain.region", all = T), 
                    list(file.elasticnet.count[,c(1,2,4)], file.welasticnet.count[,c(1,3)], file.lasso.count[,c(1,3)], file.wlasso.count[,c(1,3)]))

file.merge[,3:6] = reproducibility.index(file.merge[,3:6], k=10)

write.table(file.merge, 'coefs_count_merge_20210404.csv', sep = ',', row.names = F)
write.table(file.merge.coefs, 'coefs_merge_20210404.csv', sep = ',', row.names = F)

