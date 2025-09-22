###### PLOT VARIANT COUNTS USING R ######

#Load libraries
library(dplyr)
library(tidyverse)

#Create new column containing variants
modified <- modified %>%
mutate(variants = paste(modified$REF, ">", modified$ALT ))

#Counts the variants
modified <- modified %>%
group_by(variants) %>%
summarize(count = n())

#Plot boxplot 
barplot_2 <- barplot(counts$count, las=2, ylim = c(0, max(counts$count)+50), names.arg = counts$variants, xlab="Variants", ylab="Counts")
text(x = barplot, y = counts$count, label = counts$count, pos = 3, cex = 0.8, col = "blue")
