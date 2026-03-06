library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(tidyverse)
#----------------------------------------------------------------------------
### functions
"%&%" = function(a, b) { paste0(a, b) }
# Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
  qnorm(rank(x) / (length(x)+1))
}


# Function to order by sampleID
custom_sort <- function(s) {
  parts <- strsplit(s, "_")[[1]]
  as.integer(parts[2])
}

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# Input data for TMM calculation
# Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id

Counts = read.table(args[1], row.names = 1, header = T)

# TPM matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
TPM = read.table(args[2], row.names = 1, header = T)

## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number

# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds
count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2
sample_count_threshold = 10

#keep the genes with >=0.1 tpm and >=6 read counts in >=20% samples.
expr_tpm = TPM[rownames(expr_counts),samids]
tpm_th = rowSums(expr_tpm >= tpm_threshold)
count_th = rowSums(expr_counts >= count_threshold)
ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl1 & ctrl2
TMM_pass = TMM[mask,] ##row is gene; column is sample

###expression values (TMM) were inverse normal transformed across samples.
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
#----------------------------------------------------------------------------
### 2. prepare bed file

annotation <- read.table(gzfile(args[3]), header = F, sep = "\t")
colnames(annotation) <- c("Chr", "Start", "End", "gene_id")

# Modify annotation file
# Need to extract gene information etc if want to convert to symbols
geneid = annotation$gene_id

expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,] # expr_matrix TMM_inv

# prepare bed file for tensorQTL
bed_annot = annotation[annotation$gene_id %in% rownames(expr_matrix),]
bed = data.frame(bed_annot,expr_matrix[bed_annot$gene_id,])
bed = bed[bed[,1] %in% as.character(paste0("chr",1:22)),]
bed[,1] = gsub("chr", "", bed[,1])
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
bed[,1] = gsub("^", "chr", bed[,1])
colnames(bed)[1] = "#Chr"
# output bed file
write.table(bed,file = args[4], sep = "\t", row.names = F, quote = F)

#--------------------------------------------------------------------------------
### 3. Identify hidden covariates

expr <- t(bed[,-(1:4)]) # CHR , START, END, GeneID - all removed before transposing

print(dim(expr))


# Run PCA
prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE) #This should take less than a minute.
PCs<-prcompResult$x #368*368. The columns are PCs.

# Choose K
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)


RNGkind("L'Ecuyer-CMRG")
set.seed(1)
resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05)
print(resultRunBE$numOfPCsChosen)


K_elbow<-resultRunElbow #12.
K_BE<-resultRunBE$numOfPCsChosen #29.
K_Max <- 12 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE","MAX"),values=c(K_elbow,K_BE,K_Max),
                         titleText="")



ggplot2::ggsave(args[5],width=16,height=11,unit="cm", dpi = 600)

#------------------------------------------------------------------------------------
#### 4. Integrate with known covariates


cov <- read.table(args[6], sep = "\t", header = T) %>% select(-c(Group, subgroup)) # uninformative
# convert to factor and numeric
print(head(cov))
cov$gender <- factor(cov$gender, levels = c("male", "female"))

# Specify the levels for the "smear_results" column
cov$smear_results <- factor(cov$smear_results, levels = c("Positive", "Negative"))


cov$subgroup_att <- factor(cov$subgroup_att, levels = c("Outbreak TB strain", "Short_ATT", "Long_ATT", "Difficult TB Cases", "TB Drug Resistance"))

cov$days_from_att <- scale(as.numeric(cov$days_from_att), center = TRUE)

cov$Age <- scale(as.numeric(cov$Age), center = TRUE)
cov$response_group <- factor(cov$response_group)


pca <- read.table(args[7], sep = "\t", header = T) %>% select(1:4) # top 3 genotype PCS
known <- left_join(pca, cov, by = c("SampleID" = "Patient_ID")) # join
known_sorted <- known %>% arrange(sapply(SampleID, custom_sort)) # sort so in same order as PC
rownames(known_sorted) <- known_sorted$SampleID
known_sorted <- known_sorted %>% select(-1)

known_sorted <- known_sorted[rownames(expr),] # get in same order to avoid any doubt.

dim(known_sorted) # check dimensions

identical(rownames(known_sorted),rownames(expr)) #TRUE is good.
print(rownames(known_sorted)) # check order
print(rownames(expr)) # check order

colnames(known_sorted)[1:3]<-paste0("genotypePC",1:3) 


# Select PCS from BE analysis to use (upper bound)

PCsTop<-PCs[,1:K_elbow] #48*K_elbow.

print(dim(PCsTop))

# Remove highly correlated features
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(known_sorted,PCsTop,unadjustedR2_cutoff=0.9)

PCsTop<-scale(PCsTop) #Optional. Could be helpful for avoiding numerical inaccuracies.
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)

# transpose these
write.table(t(covariatesToUse), args[8], sep = "\t", row.names = T, col.names = T, quote = F)
