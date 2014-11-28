#!~/bin/Rscript

# 06_EDA.R

options(echo = TRUE)

load("../data/05_gexnormed_cpv_matched.rda")

pdf("../results/06_EDA.pdf", height = 10, width = 10)

library(MASS)

sel.hisd = apply(gex, 1, sd) > quantile(apply(gex, 1, sd), 0.9)
gex.hisd = gex[sel.hisd,]
gex.hisd.std = (gex.hisd - rowMeans(gex.hisd)) / apply(gex.hisd, 1, sd)
dist.hisd.std = dist(t(gex.hisd.std))
mds.hisd.std = isoMDS(dist.hisd.std)

library(ggplot2)
plot_data = data.frame(mds_x = mds.hisd.std$points[,1], mds_y = mds.hisd.std$points[,2], histo_type = abbreviate(cpvs$Path.HistoType), purity = samples$purity_qpure)
ggplot(plot_data, aes(x = mds_x, y = mds_y, colour = histo_type)) + geom_point() + ggtitle("Histological subtype clusters")
ggplot(plot_data[!is.na(plot_data$purity),], aes(x = mds_x, y = mds_y, colour = purity)) + geom_point() + ggtitle("Purity effect")
ggplot(plot_data, aes(x = mds_x, y = mds_y, colour = histo_type, size = purity)) + geom_point() + ggtitle("Histological subtype clusters and purity effect")

dev.off()

sessionInfo()
