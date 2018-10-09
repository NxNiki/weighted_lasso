#!/usr/bin/env Rscript


mf.catvbm_alff_reho_dti.bn246 = scale(cbind(cat.vbm.neuromorph, alff.bn246, reho.bn246,
                                                                 label.fa, tract.fa, label.md, tract.md))

mf.catvbm_falff_reho_dti.bn246 = scale(cbind(cat.vbm.neuromorph, falff.bn246, reho.bn246, 
                                                                  label.fa, tract.fa, label.md, tract.md))

mf.dti = scale(cbind(label.fa, label.md, tract.fa, tract.md))
dti.fa = scale(cbind(label.fa, tract.fa))
dti.md = scale(cbind(label.md, tract.md))

mf.catvbm.dti = scale(cbind(cat.vbm.neuromorph, label.fa, label.md, tract.fa, tract.md))

mf.reho.alff.bn246 = scale(cbind(alff.bn246, reho.bn246))
mf.reho.falff.bn246 = scale(cbind(reho.bn246, falff.bn246))

mf.catvbm.fa = scale(cbind(cat.vbm.neuromorph, label.fa, tract.fa))
mf.catvbm.md = scale(cbind(cat.vbm.neuromorph, label.md, tract.md))

mf.catvbm.alff.reho.bn246 = scale(cbind(cat.vbm.neuromorph, alff.bn246, reho.bn246))
mf.catvbm.falff.reho.bn246 = scale(cbind(cat.vbm.neuromorph, falff.bn246, reho.bn246))

#feature.list = list(list(cat.vbm.neuromorph, "cat.vbm.neuromorph"),
#                    list(alff.bn246, "alff.bn246"),
#                    list(falff.bn246, "falff.bn246"),
#                    list(reho.bn246, "reho.bn246"),
#                    list(label.fa, "label.fa"),
#                    list(tract.fa, "tract.fa"),
#                    list(label.md, "label.md"),
#                    list(tract.md, "tract.md"),
#                    list(dti.fa, "FA"),
#                    list(dti.md, "MD"),
#                    list(mf.dti,"FA_MD"),
#                    list(mf.catvbm.dti, "catvbm_fa_md"),
#                    list(mf.reho.alff.bn246, "reho_alff_bn246"),
#                    list(mf.reho.falff.bn246, "reho_falff_bn246"),
#                    list(mf.catvbm.fa, "catvbm_fa"),
#                    list(mf.catvbm.md, "catvbm_md"),
#                    list(mf.catvbm.alff.reho.bn246, "catvbm_alff_reho"),
#                    list(mf.catvbm.falff.reho.bn246, "catvbm_falff_reho"),
#                    list(mf.catvbm_alff_reho_dti.bn246, "multimodal.alff.bn246"),
#                    list(mf.catvbm_falff_reho_dti.bn246, "multimodal.falff.bn246")
#                    )

feature.list = list(list(cat.vbm.neuromorph, "cat.vbm.neuromorph"),
                    list(alff.bn246, "alff.bn246"),
                    list(reho.bn246, "reho.bn246"),
                    list(dti.fa, "FA"),
                    list(dti.md, "MD"),
                    list(mf.catvbm_alff_reho_dti.bn246, "multimodal.alff.bn246")
)

