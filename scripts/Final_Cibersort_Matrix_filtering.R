library(data.table)
markers <- read.table(file='/home/workspace/jogrady/heQTL/work/scRNA_seq/Conserved_Marker_genes.txt', sep='\t')

expression_matrix <- fread(file="/home/workspace/jogrady/heQTL/work/scRNA_seq/preSigMatrix.tsv")      

