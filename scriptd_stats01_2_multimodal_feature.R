#!/usr/bin/env Rscript


multimodal.feature = scale(cbind(spm.vbm, alff, reho, label.fa, tract.fa, label.md, tract.md))
multimodal.feature.aal = scale(cbind(spm.vbm.aal, alff.aal, reho.aal, label.fa, tract.fa, label.md, tract.md))
multimodal.feature.bn246 = scale(cbind(spm.vbm.bn246, alff.bn246, reho.bn246, label.fa, tract.fa, label.md, tract.md))

mf.dti = scale(cbind(label.fa, label.md, tract.fa, tract.md))
dti.fa = scale(cbind(label.fa, tract.fa))
dti.md = scale(cbind(label.md, tract.md))

mf.vbm.dti = scale(cbind(spm.vbm, label.fa, label.md, tract.fa, tract.md))
mf.vbm.aal.dti = scale(cbind(spm.vbm.aal, label.fa, label.md, tract.fa, tract.md))
mf.vbm.bn246.dti = scale(cbind(spm.vbm.bn246, label.fa, label.md, tract.fa, tract.md))

mf.resting = scale(cbind(alff, reho))
mf.resting.aal = scale(cbind(alff.aal, reho.aal, falff.aal))
mf.resting.bn246 = scale(cbind(alff.bn246, reho.bn246, falff.bn246))

mf.vbm.fa = scale(cbind(spm.vbm, label.fa, tract.fa))
mf.vbm.aal.fa = scale(cbind(spm.vbm.aal, label.fa, tract.fa))
mf.vbm.bn246.fa = scale(cbind(spm.vbm.bn246, label.fa, tract.fa))

mf.vbm.resting = scale(cbind(spm.vbm, alff, reho))
mf.vbm.resting.aal = scale(cbind(spm.vbm.aal, alff.aal, reho.aal))
mf.vbm.resting.bn246 = scale(cbind(spm.vbm.bn246, alff.aal, reho.aal))

#feature.name = c("no.brain","cat.vbm.hammers", "cat.vbm.neuromorph", "spm.vbm", "spm.vbm.aal", "spm.vbm.bn246", "alff", "alff.aal", "alff.bn246", "falff", "falff.aal", "falff.bn246", "reho", "reho.aal", "reho.bn246", "label.fa", "tract.fa", "dti.fa", "multimodal.feature", "multimodal.feature.aal", "multimodal.feature.bn246", "mf.dti", "mf.vbm.dti", "mf.vbm.aal.dti", "fm.vbm.bn246.dti", "mf.resting", "mf.resting.aal", "mf.resting.bn246", "mf.vbm.fa", "mf.vbm.aal.fa", "mf.vbm.bn246.fa", "mf.vbm.resting", "mf.vbm.resting.aal", "mf.vbm.resting.bn246")
#
#feature.list = list(NULL, cat.vbm.hammers, cat.vbm.neuromorph, spm.vbm, spm.vbm.aal, spm.vbm.bn246, alff, alff.aal, alff.bn246, falff, falff.aal, falff.bn246, reho, reho.aal, reho.bn246, label.fa, tract.fa, dti.fa, multimodal.feature, multimodal.feature.aal, multimodal.feature.bn246, mf.dti, mf.vbm.dti, mf.vbm.aal.dti, mf.vbm.bn246.dti, mf.resting, mf.resting.aal, mf.resting.bn246, mf.vbm.fa, mf.vbm.aal.fa, mf.vbm.bn246.fa, mf.vbm.resting, mf.vbm.resting.aal, mf.vbm.resting.bn246)

#feature.name = c("cat.vbm.hammers", "cat.vbm.neuromorph", "spm.vbm", "spm.vbm.aal", "spm.vbm.bn246", "alff", "alff.aal", "alff.bn246", "falff", "falff.aal", "falff.bn246", "reho", "reho.aal", "reho.bn246", "label.fa", "tract.fa", "dti.fa", "multimodal.feature", "multimodal.feature.aal", "multimodal.feature.bn246", "mf.dti", "mf.vbm.dti", "mf.vbm.aal.dti", "fm.vbm.bn246.dti", "mf.resting", "mf.resting.aal", "mf.resting.bn246", "mf.vbm.fa", "mf.vbm.aal.fa", "mf.vbm.bn246.fa", "mf.vbm.resting", "mf.vbm.resting.aal", "mf.vbm.resting.bn246")

#feature.list = list(cat.vbm.hammers, cat.vbm.neuromorph, spm.vbm, spm.vbm.aal, spm.vbm.bn246, alff, alff.aal, alff.bn246, falff, falff.aal, falff.bn246, reho, reho.aal, reho.bn246, label.fa, tract.fa, dti.fa, multimodal.feature, multimodal.feature.aal, multimodal.feature.bn246, mf.dti, mf.vbm.dti, mf.vbm.aal.dti, mf.vbm.bn246.dti, mf.resting, mf.resting.aal, mf.resting.bn246, mf.vbm.fa, mf.vbm.aal.fa, mf.vbm.bn246.fa, mf.vbm.resting, mf.vbm.resting.aal, mf.vbm.resting.bn246)

feature.list = list(list(spm.vbm, "spm.vbm"),
                    list(cat.vbm.neuromorph, "cat.vbm.neuromorph"),
                    list(alff.bn246, "alff.bn246"),
                    list(reho.bn246, "reho.bn246"),
                    list(dti.fa, "FA"),
                    list(dti.md, "MD"),
                    list(multimodal.feature.bn246, "multimodal.alff.bn246")
)
#feature.name = feature.name[c(14,16,20)]
#feature.list = feature.list[c(14,16,20)]

#feature.name = c("fc")
#feature.list = list(fc)


