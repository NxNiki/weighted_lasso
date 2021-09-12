#!/usr/bin/env Rscript


multimodal.feature = scale(cbind(cat.vbm.neuromorph, alff, reho, label.fa, tract.fa, label.md, tract.md))
multimodal.feature.aal = scale(cbind(spm.vbm.aal, alff.aal, reho.aal, label.fa, tract.fa, label.md, tract.md))
multimodal.feature.bn246 = scale(cbind(spm.vbm.bn246, alff.bn246, reho.bn246, label.fa, tract.fa, label.md, tract.md))

mf.dti = scale(cbind(label.fa, label.md, tract.fa, tract.md))
dti.fa = scale(cbind(label.fa, tract.fa))

mf.vbm.dti = scale(cbind(cat.vbm.neuromorph, label.fa, label.md, tract.fa, tract.md))
mf.vbm.aal.dti = scale(cbind(spm.vbm.aal, label.fa, label.md, tract.fa, tract.md))
mf.vbm.bn246.dti = scale(cbind(spm.vbm.bn246, label.fa, label.md, tract.fa, tract.md))

mf.resting = scale(cbind(alff, reho))
mf.resting.aal = scale(cbind(alff.aal, reho.aal))
mf.resting.bn246 = scale(cbind(alff.bn246, reho.bn246))

mf.vbm.fa = scale(cbind(cat.vbm.neuromorph, label.fa, tract.fa))
mf.vbm.aal.fa = scale(cbind(spm.vbm.aal, label.fa, tract.fa))
mf.vbm.bn246.fa = scale(cbind(spm.vbm.bn246, label.fa, tract.fa))

mf.vbm.resting = scale(cbind(cat.vbm.neuromorph, alff, reho))
mf.vbm.resting.aal = scale(cbind(spm.vbm.aal, alff.aal, reho.aal))
mf.vbm.resting.bn246 = scale(cbind(spm.vbm.bn246, alff.aal, reho.aal))

#feature.name = c("no.brain","cat.vbm.hammers", "cat.vbm.neuromorph", "spm.vbm", "spm.vbm.aal", "spm.vbm.bn246", "alff", "alff.aal", "alff.bn246", "falff", "falff.aal", "falff.bn246", "reho", "reho.aal", "reho.bn246", "label.fa", "tract.fa", "dti.fa", "multimodal.feature", "multimodal.feature.aal", "multimodal.feature.bn246", "mf.dti", "mf.vbm.dti", "mf.vbm.aal.dti", "fm.vbm.bn246.dti", "mf.resting", "mf.resting.aal", "mf.resting.bn246", "mf.vbm.fa", "mf.vbm.aal.fa", "mf.vbm.bn246.fa", "mf.vbm.resting", "mf.vbm.resting.aal", "mf.vbm.resting.bn246")
#
#feature.list = list(NULL, cat.vbm.hammers, cat.vbm.neuromorph, spm.vbm, spm.vbm.aal, spm.vbm.bn246, alff, alff.aal, alff.bn246, falff, falff.aal, falff.bn246, reho, reho.aal, reho.bn246, label.fa, tract.fa, dti.fa, multimodal.feature, multimodal.feature.aal, multimodal.feature.bn246, mf.dti, mf.vbm.dti, mf.vbm.aal.dti, mf.vbm.bn246.dti, mf.resting, mf.resting.aal, mf.resting.bn246, mf.vbm.fa, mf.vbm.aal.fa, mf.vbm.bn246.fa, mf.vbm.resting, mf.vbm.resting.aal, mf.vbm.resting.bn246)

feature.name = c("cat.vbm.hammers", "cat.vbm.neuromorph", "spm.vbm", "spm.vbm.aal", 
                 "spm.vbm.bn246", "alff", "alff.aal", "alff.bn246", "falff", "falff.aal", 
                 "falff.bn246", "reho", "reho.aal", "reho.bn246", "label.fa", "tract.fa",
                 "dti.fa", "multimodal.feature", "multimodal.feature.aal", "multimodal.feature.bn246", 
                 "mf.dti", "mf.vbm.dti", "mf.vbm.aal.dti", "fm.vbm.bn246.dti", "mf.resting", "mf.resting.aal", 
                 "mf.resting.bn246", "mf.vbm.fa", "mf.vbm.aal.fa", "mf.vbm.bn246.fa", "mf.vbm.resting", 
                 "mf.vbm.resting.aal", "mf.vbm.resting.bn246")

feature.list = list(cat.vbm.hammers , cat.vbm.neuromorph, spm.vbm, spm.vbm.aal, 
                    spm.vbm.bn246, alff, alff.aal, alff.bn246, falff, falff.aal, 
                    falff.bn246, reho, reho.aal, reho.bn246, label.fa, tract.fa,
                    dti.fa, multimodal.feature, multimodal.feature.aal, multimodal.feature.bn246, 
                    mf.dti, mf.vbm.dti, mf.vbm.aal.dti, mf.vbm.bn246.dti, mf.resting, mf.resting.aal, 
                    mf.resting.bn246, mf.vbm.fa, mf.vbm.aal.fa, mf.vbm.bn246.fa, mf.vbm.resting, 
                    mf.vbm.resting.aal, mf.vbm.resting.bn246)

feature.list2 = list()
for (i in 1:length(feature.name)){
    feature.list2 = c(feature.list2, list(list(feature.list[[i]], feature.name[i])))
}

feature.list = feature.list2

#feature.name = feature.name[c(14,16,20)]
#feature.list = feature.list[c(14,16,20)]

#feature.name = c("fc")
#feature.list = list(fc)


