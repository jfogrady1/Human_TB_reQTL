# Script to transpose interaction file and get interaction file in right format

library(data.table)
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



args = commandArgs(trailingOnly=TRUE)
data <- fread(args[1])


head(data)
head(data)
data <- data %>% select(1:16)
head(data)
data <- as.data.frame(data)
hist(data$`Classical Monocyte`)
dim(data)

colMeans(data[2:16])
data_1 <- data[1:48,]
data_2 <- data[49:96,]
data_3 <- data[97:144,]
data_4 <- data[145:192,]
data_5 <- data[193:240,]


data_1_medians <- apply(data_1[,-1], 2,median)
data_2_medians <- apply(data_2[,-1], 2, median)
data_3_medians <- apply(data_3[,-1], 2, median)
data_4_medians <- apply(data_4[,-1], 2, median)
data_5_medians <- apply(data_5[,-1], 2, median)


data_1_medians <- data_1_medians[data_1_medians > 0.01]
data_2_medians <- data_2_medians[data_2_medians > 0.01]
data_3_medians <- data_3_medians[data_3_medians > 0.01]
data_4_medians <- data_4_medians[data_4_medians > 0.01]
data_5_medians <- data_5_medians[data_5_medians > 0.01]

median(data_3$`Classical Monocyte`)

data_1_medians
data_2_medians 
data_3_medians
data_4_medians
data_5_medians

data_1$Mixture <- gsub("_T0", "", data_1$Mixture)
data_2$Mixture <- gsub("_T1", "", data_2$Mixture)
data_3$Mixture <- gsub("_T2", "", data_3$Mixture)
data_4$Mixture <- gsub("_T3", "", data_4$Mixture)
data_5$Mixture <- gsub("_T4", "", data_5$Mixture)
data_1$Mixture
final_cell_types = intersect(names(data_1_medians), intersect(names(data_2_medians), intersect(names(data_3_medians), intersect(names(data_4_medians), names(data_5_medians)))))
final_cell_types <- final_cell_types[1:6] # 6 cell types
final_cell_types


data_2$Mixture <- data_1$Mixture
data_3$Mixture <- data_1$Mixture
data_4$Mixture <- data_1$Mixture
data_5$Mixture <- data_1$Mixture

celltypes <- c(args[2], args[3], args[4], args[5], args[6])

for (celltype in celltypes) {
  
  vec <- c("Mixture", celltype)
  vec <- gsub("_", " ", vec)
  data_1_temp <- data_1[,vec]
  data_2_temp <- data_2[,vec]
  data_3_temp <- data_3[,vec]
  data_4_temp <- data_4[,vec]
  data_5_temp <- data_5[,vec]
  
  colnames(data_1_temp)[2] <- gsub(" ", "_", colnames(data_1_temp)[2])
  colnames(data_2_temp)[2] <- gsub(" ", "_", colnames(data_2_temp)[2])
  colnames(data_3_temp)[2] <- gsub(" ", "_", colnames(data_3_temp)[2])
  colnames(data_4_temp)[2] <- gsub(" ", "_", colnames(data_4_temp)[2])
  colnames(data_5_temp)[2] <- gsub(" ", "_", colnames(data_5_temp)[2])
  
  vec[2] <- gsub(" ", "_", vec[2])
  print(vec[2])
  celltype <- gsub(" ", "_", celltype)
  data_1_temp[,vec[2]] <- inverse_normal_transform(data_1_temp[,vec[2]])
  data_2_temp[,vec[2]] <- inverse_normal_transform(data_2_temp[,vec[2]])
  data_3_temp[,vec[2]] <- inverse_normal_transform(data_3_temp[,vec[2]])
  data_4_temp[,vec[2]] <- inverse_normal_transform(data_4_temp[,vec[2]])
  data_5_temp[,vec[2]] <- inverse_normal_transform(data_5_temp[,vec[2]])
  
  
  write.table(data_1_temp, file = paste0("/home/workspace/jogrady/heQTL/work/ieQTL/T0_",celltype,"_interaction_input_from_cibersort.txt"), sep = "\t", row.names = F, quote = F)
  write.table(data_2_temp, file = paste0("/home/workspace/jogrady/heQTL/work/ieQTL/T1_",celltype,"_interaction_input_from_cibersort.txt"), sep = "\t", row.names = F, quote = F)
  write.table(data_3_temp, file = paste0("/home/workspace/jogrady/heQTL/work/ieQTL/T2_",celltype,"_interaction_input_from_cibersort.txt"), sep = "\t", row.names = F, quote = F)
  write.table(data_4_temp, file = paste0("/home/workspace/jogrady/heQTL/work/ieQTL/T3_",celltype,"_interaction_input_from_cibersort.txt"), sep = "\t", row.names = F, quote = F)
  write.table(data_5_temp, file = paste0("/home/workspace/jogrady/heQTL/work/ieQTL/T4_",celltype,"_interaction_input_from_cibersort.txt"), sep = "\t", row.names = F, quote = F)
  
}

write.table(data_1, file = args[7], sep = "\t", row.names = F, quote = F)
write.table(data_2, file = args[8], sep = "\t", row.names = F, quote = F)
write.table(data_3, file = args[9], sep = "\t", row.names = F, quote = F)
write.table(data_4, file = args[10], sep = "\t", row.names = F, quote = F)
write.table(data_5, file = args[11], sep = "\t", row.names = F, quote = F)

