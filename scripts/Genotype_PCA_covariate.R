library(SNPRelate)
args <- commandArgs(trailingOnly = TRUE)
#5. Genotype PCA

vcf.fn = args[1]
snpgdsVCF2GDS(vcf.fn, args[2], method = "biallelic.only")
genofile <- snpgdsOpen(args[2])
ccm_pca <- snpgdsPCA(genofile,num.thread=4)

pca_genotype <- ccm_pca$eigenvect[, 1:20]
colnames(pca_genotype) <- paste0("pc", 1:20)
rownames(pca_genotype) <- ccm_pca$sample.id
pca_genotype0 = data.frame(SampleID=ccm_pca$sample.id,pca_genotype)
pca_var0 = data.frame(pc=1:20,eigenval=ccm_pca$eigenval[1:20],varprop=ccm_pca$varprop[1:20])

# output
write.table(pca_genotype0, args[3], sep = "\t", row.names = F, quote = FALSE)
write.table(pca_var0, args[4], sep = "\t", row.names = F, quote = FALSE)
