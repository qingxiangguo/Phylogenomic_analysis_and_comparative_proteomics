library(UpSetR)


data <- read.csv("eight_matrix_strict.csv",header=T)

pdf("eight_orthomcl.pdf",  width=11, height=6)

upset(data, nsets = 8, order.by="freq", sets.bar.color = "#56B4E9", mainbar.y.label = "OG Intersections", sets.x.label = "Proteins per Species",nintersects = 70)


dev.off()
