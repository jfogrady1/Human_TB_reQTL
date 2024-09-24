# read the count data
# Read in all of the files
counts0 <- as.matrix(read.table("count_matrix_ordered_T0.txt", sep = "\t", header = T, row.names = 1))
counts1 <- as.matrix(read.table("count_matrix_ordered_T1.txt", sep = "\t", header = T, row.names = 1))
counts2 <- as.matrix(read.table("count_matrix_ordered_T2.txt", sep = "\t", header = T, row.names = 1))
counts3 <- as.matrix(read.table("count_matrix_ordered_T3.txt", sep = "\t", header = T, row.names = 1))
counts4 <- as.matrix(read.table("count_matrix_ordered_T4.txt", sep = "\t", header = T, row.names = 1))

colnames(counts0) <- paste0(colnames(counts0),"_T0")
colnames(counts1) <- paste0(colnames(counts1),"_T1")
colnames(counts2) <- paste0(colnames(counts2),"_T2")
colnames(counts3) <- paste0(colnames(counts3),"_T3")
colnames(counts4) <- paste0(colnames(counts4),"_T4")


counts.table = cbind(counts0, counts1, counts2, counts3, counts4)
