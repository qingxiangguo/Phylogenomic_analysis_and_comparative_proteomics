library(VennDiagram)
a<-read.csv("./mhon_OG",header=F,sep="\t")
b<-read.csv("./mwul_OG",header=F,sep="\t")
c<-read.csv("./tkit_OG",header=F,sep="\t")

a <- a[[1]]
b <- b[[1]]
c <- c[[1]]

venn.plot <- venn.diagram(x=list(M.honghuensis=a, M.wulii=b, T.kitauei=c), filename = "Venn_3set_2.png", fill=c("dodgerblue","darkorange1","orchid3"), imagetype = "png", alpha = 0.50, cex = 1.5, cat.col = c("white","white","white"), margin = 0.05,  cat.cex = 0.1, col = "black", fontfamily = "sans", overrideTriple = 1, euler.d = TRUE, scaled = TRUE);

