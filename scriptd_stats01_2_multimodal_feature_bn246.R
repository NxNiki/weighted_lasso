#!/usr/bin/env Rscript


#multimodal.feature.bn246 = scale(cbind(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, label.fa, tract.fa, label.md, tract.md))
multimodal.feature.bn246 = scale(cbind(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, label.fa, tract.fa))

mf.dti = scale(cbind(label.fa, label.md, tract.fa, tract.md))
dti.fa = scale(cbind(label.fa, tract.fa))
dti.md = scale(cbind(label.md, tract.md))

mf.vbm.bn246.dti = scale(cbind(spm.vbm.bn246, label.fa, label.md, tract.fa, tract.md))

mf.resting.bn246 = scale(cbind(alff.bn246, reho.bn246, falff.bn246))

mf.vbm.bn246.fa = scale(cbind(spm.vbm.bn246, label.fa, tract.fa))

mf.vbm.resting.bn246 = scale(cbind(spm.vbm.bn246, alff.aal, reho.aal))

#feature.name = c("spm.vbm.bn246", "alff.bn246", "falff.bn246", "reho.bn246", "dti.fa", "dti.md", "multimodal.feature.bn246", "mf.dti", "mf.vbm.bn246.dti", "mf.resting.bn246", "mf.vbm.bn246.fa", "mf.vbm.resting.bn246")
#feature.list = list(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, dti.fa, dti.md, multimodal.feature.bn246, mf.dti, mf.vbm.bn246.dti, mf.resting.bn246, mf.vbm.bn246.fa, mf.vbm.resting.bn246)

#feature.name = c("spm.vbm.bn246", "alff.bn246", "falff.bn246", "reho.bn246", "dti.fa", "dti.md", "multimodal.feature.bn246", "FC")
#feature.list = list(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, dti.fa, dti.md, multimodal.feature.bn246, fc)

feature.name = c("spm.vbm.bn246", "alff.bn246", "falff.bn246", "reho.bn246", "dti.fa", "dti.md", "multimodal.feature.bn246")
feature.list = list(spm.vbm.bn246, alff.bn246, falff.bn246, reho.bn246, dti.fa, dti.md, multimodal.feature.bn246)


