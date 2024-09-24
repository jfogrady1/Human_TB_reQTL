## Haeder, Schaeuble et al., 2023, Nature Communications
## perform response eQTL anaysis
# input: results from QTLtools - nominal path
# output: reQTLs for all stimuli and timepoints

library(data.table)
library(tidyverse)
args<-commandArgs(T)



# output files
outputfile<-paste0("/home/workspace/jogrady/heQTL/results/reQTLs/",args[3],"_",args[4],"_","reQTL_results.txt")
outputfile_sig<-paste0("/home/workspace/jogrady/heQTL/results/reQTLs/",args[3],"_",args[4],"_","reQTL_sig_results.txt")
#file for gene annotation
genes_anno<- args[5]


# get nominal associations for reference set which we will compare top eQTLs in our data set against
base_name<-args[1]
a_file<-paste0(args[1]) # baseline

stim_name<-args[2]
b_file<-paste0(args[2])



# read eQTL files
a<-fread(a_file)
# need beta and se

a_info<-data.frame(gene_id=a$phenotype_id,snp=a$variant_id,beta=a$slope,var=abs(a$slope_se^2))
a_info$association <- paste0(a_info$snp, "-", a_info$gene_id)
b <- fread(b_file) %>% filter(is_eGene == TRUE)
head(b)
b_info<-data.frame(gene_id=b$phenotype_id,snp=b$variant_id,beta=b$slope,var=abs(b$slope_se^2))
b_info$association <- paste0(b_info$snp,"-",b_info$gene_id)
a_info <- a_info %>% filter(association %in% b_info$association)

# cross reference
b_info <- b_info %>% filter(association %in% a_info$association)
all(b_info$association == a_info$association)


# Now sort

a_info <- a_info[order(a_info$gene_id),]
b_info <- b_info[order(b_info$gene_id),]

all(a_info$gene_id == b_info$gene_id)

dim(a_info)
dim(b_info)
all(b_info$snp == a_info$snp)
# calculate z-score
a_info$beta

zscore<-(a_info$beta - b_info$beta)/(sqrt(a_info$var+b_info$var))

hist(zscore)
# calculate p-value for response eQTL
pval<-2*pnorm(-abs(zscore))
pval_corr<-p.adjust(pval,method = "BH")
min(pval_corr)
hist(pval_corr, breaks = 10)
0.05/500

# get all information in one file
symbols <- fread(args[5]) %>% filter(V3 == "gene")
symbols <- symbols %>% separate(V9, into = c("gene_id","gene_type","gene_name"), sep = ";")
symbols$gene_id <- gsub('gene_id "', '', symbols$gene_id)
symbols$gene_id <- gsub('"', '', symbols$gene_id)
symbols$gene_id <- gsub(' ', '', symbols$gene_id)
symbols$gene_name <- gsub('gene_name "', '', symbols$gene_name)
symbols$gene_name <- gsub('"', '', symbols$gene_name)
symbols$gene_name <- gsub(' ', '', symbols$gene_name)
head(symbols)
matched_genes<-symbols[match(a_info$gene_id,symbols$gene_id),]
matched_genes
response<-data.frame(matched_genes,snp_id=a_info$snp,a_info[,3:4],b_info[,3:4],zscore,pval,pval_corr)
response <- response %>% select(-c(V6, V8, V2, V3))
colnames(response)<-c("chr", "gene_start", "gene_end", "strand", "ensembl_gene_id", "gene_type", "gene_name", "snp_id", "beta_baseline", "var_baseline", "beta_stimulated", "var_stimulated", "zscore", "pval", "Padj")
response$comparison <- paste0(args[3],"_V_",args[4])
response_sort<-response[order(response$Padj),]
response_sort_sig<-response_sort[response_sort$Padj<=0.05,]
response_sort_sig$gene_name[response_sort_sig$hgnc==""]<-NA
dim(response_sort_sig)
head(response_sort_sig)
write.table(response_sort,outputfile, quote = F, row.names = F)
write.table(response_sort_sig,outputfile_sig, quote = F, row.names = F)




