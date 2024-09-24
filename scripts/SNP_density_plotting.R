library(tidyverse)
library(karyoploteR)
library(BiocManager)
library(plyranges)
args = commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], sep = "\t")
data <- data %>% select(-4) 
data <- data %>% separate(., V5, into = c("DR2", "AF", "IMP", "AC", "AN"), sep = ";")
data <- data %>% select(1,2,3,6)
data <- data %>% mutate(IMP = if_else(IMP == "IMP", "IMP", "TYPED"))

data_typed <- data %>% filter(data$IMP == "TYPED")
data_imputed <- data %>% filter(data$IMP == "IMP")

data_typed$start <- as.numeric(data_typed$V3, -1)
colnames(data_typed)[3] <- "end"
colnames(data_typed)[1] <- "snp"
colnames(data_typed)[2] <- "chr"
data_typed$chr <- gsub("^", "chr", data_typed$chr)


data_imputed$start <- as.numeric(data_imputed$V3, -1)
colnames(data_imputed)[3] <- "end"
colnames(data_imputed)[1] <- "snp"
colnames(data_imputed)[2] <- "chr"
data_imputed$chr <- gsub("^", "chr", data_imputed$chr)









data_imputedGR <- makeGRangesFromDataFrame(data_imputed, start.field = "start", end.field = "end" )
data_typedGR <- makeGRangesFromDataFrame(data_typed, start.field = "start", end.field = "end", )




pdf(args[2], width = 15, height = 12)

# All chromosomes, 500kb bins
kp <- plotKaryotype(genome = "hg38", plot.type = 1, chromosomes = c(paste0("chr", c(1:5,6,7:22))), labels.plotter = NULL)
kpAddChromosomeNames(kp, chr.names = c(paste0("chr", c(1:5,6,7:22))), xoffset = -0.02, yoffset = 50)
kp <- kpPlotDensity(kp, data=data_imputedGR, window.size = 5e5, border="#0e87eb", col="#0e87eb", r1=0.95)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.95, cex=0.525)
kp <- kpPlotDensity(kp, data=data_typedGR, window.size = 5e5, border="black", col="#027148", apha = 0.5, r1=0.95)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.95, cex=0.525, side = 2)

dev.off()

pdf(args[3], width = 15, height = 12)
# chromosome 12 through 16, 500 kb bins  
kp_2 <- plotKaryotype(genome = "hg38", plot.type=1, chromosomes = c(paste0("chr", c(11,12,13,14,15))))
kpRect(kp_2, chr="chr13", x0=55e6, x1=81e6, y0=-0.1, y1=0.5, col = "grey")
kp_2 <- kpPlotDensity(kp_2, data=data_imputedGR, r0=0, r1=0.95, window.size = 5e5, border="#0e87eb", col="#0e87eb")
kpAxis(kp_2, ymax=kp_2$latest.plot$computed.values$max.density, r0=0, r1=0.95, cex=0.55, numticks = 6)
kp_2 <- kpPlotDensity(kp_2, data=data_typedGR, r0=0, r1=0.95, window.size = 5e5, border="black", col="#027148")
kpAxis(kp_2, ymax=kp_2$latest.plot$computed.values$max.density, r0=0, r1=0.95, cex=0.55, side = 2, numticks = 6)
kpAddBaseNumbers(kp_2)

dev.off()
# Chromosome 3 interesting region
data_typedGR <- data_typedGR %>% plyranges::filter(start >= 55e6 ) %>% plyranges::filter(end <= 81e6) %>% plyranges::filter(seqnames == "chr13")
data_imputedGR <- data_imputedGR %>% plyranges::filter(start >= 55e6 ) %>% plyranges::filter(end <= 81e6) %>% plyranges::filter(seqnames == "chr13")

zoom<-toGRanges("chr13",55e6,82e6)


pdf(args[4], width = 15, height = 12)
kp_3 <- plotKaryotype(genome = "hg38", plot.type=1, zoom = zoom)
kp_3 <- kpPlotDensity(kp_3, data=data_imputedGR, r0=0, r1=1, window.size = 1e5, border="#0e87eb", col="#0e87eb")
kpAxis(kp_3, ymax=kp_3$latest.plot$computed.values$max.density, r0=0, r1=1, cex=0.7, numticks = 6, labels = c("0", "34", "68", "103", "137", "171"))
kp_3 <- kpPlotDensity(kp_3, data=data_typedGR, r0=0, r1=1, window.size = 1e5, border="black", col="#027148")
kpAxis(kp_3, ymax=kp_3$latest.plot$computed.values$max.density, r0=0, r1=1, cex=0.7, side = 2, numticks = 6)
kpAddBaseNumbers(kp_3 ,cex = 1, add.units = T, minor.ticks = T, clipping = FALSE)
dev.off()
