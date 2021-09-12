#!/usr/bin/env Rscript
library(ggpubr)
# Create some data format
# :::::::::::::::::::::::::::::::::::::::::::::::::::

set.seed(1234)
wdata = data.frame(
   sex = factor(rep(c("F", "M"), each=200)),
      weight = c(rnorm(200, 55), rnorm(200, 58)))
	  head(wdata, 4)


data("ToothGrowth")
df <- ToothGrowth
head(df, 4)
for (i in 1:1){

plot.name=(paste("plot.test", toString(i), ".png", sep = "_"))
print(plot.name)
png(plot.name)
p <- ggboxplot(df, x = "dose", y = "len",
                 color = "dose", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
				                 add = "jitter", shape = "dose")

my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
p = p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)                   # Add global p-value

print(p)
dev.off()
}
