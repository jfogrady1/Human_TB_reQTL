args = commandArgs(trailingOnly=TRUE)
timepoint <- args[2]

data <- read.table(args[1], sep = "\t", header = T)

# Ensure that names match
colnames(data) <- c("Geneid","P_1", "P_3", "P_13", "P_29", "P_36", "P_75", "P_76", "P_83", "P_84", "P_85", "P_87", "P_92", "P_99", "P_105", "P_117", "P_124", "P_127", "P_136", "P_137", "P_146", "P_157", "P_175", "P_176", "P_201", "P_203", "P_221", "P_227", "P_254", "P_256", "P_257", "P_258", "P_265", "P_266", "P_268", "P_278", "P_294", "P_297", "P_302", "P_327", "P_335", "P_348", "P_367", "P_371", "P_389", "P_395", "P_400", "P_432", "P_461")

write.table(data, paste0("/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/count_matrix_ordered_",args[2],".txt"), row.names = F, quote = F, col.names = T, sep = "\t")
