# Script to add cell tyeps to covariates file
# This is requried for ieqtl mapping
# y ~ g + i + g:i
# TensorQTL does not do this automatically


library(tidyverse)
library(data.table)
# Read in interaction files

# NK cell
T0_NK_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_NK_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("NK_cell")) %>% as.data.frame() %>% magrittr::set_rownames("NK_cell")
T1_NK_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_NK_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("NK_cell")) %>% as.data.frame() %>% magrittr::set_rownames("NK_cell")
T2_NK_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_NK_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("NK_cell")) %>% as.data.frame() %>% magrittr::set_rownames("NK_cell")
T3_NK_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_NK_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("NK_cell")) %>% as.data.frame() %>% magrittr::set_rownames("NK_cell")
T4_NK_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_NK_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("NK_cell")) %>% as.data.frame() %>% magrittr::set_rownames("NK_cell")


read.table('/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T0_NK_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T0_ieqtlcovs_NK_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T1_NK_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T1_ieqtlcovs_NK_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T2_NK_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T2_ieqtlcovs_NK_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T3_NK_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T3_ieqtlcovs_NK_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T4_NK_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T4_ieqtlcovs_NK_cell.txt',  sep = "\t", row.names = T, quote = F)



# B cell
T0_Naïve_B_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_Naïve_B_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Naïve_B_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Naïve_B_cell")
T1_Naïve_B_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_Naïve_B_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Naïve_B_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Naïve_B_cell")
T2_Naïve_B_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_Naïve_B_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Naïve_B_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Naïve_B_cell")
T3_Naïve_B_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_Naïve_B_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Naïve_B_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Naïve_B_cell")
T4_Naïve_B_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_Naïve_B_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Naïve_B_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Naïve_B_cell")


read.table('/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T0_Naïve_B_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T0_ieqtlcovs_Naïve_B_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T1_Naïve_B_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T1_ieqtlcovs_Naïve_B_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T2_Naïve_B_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T2_ieqtlcovs_Naïve_B_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T3_Naïve_B_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T3_ieqtlcovs_Naïve_B_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T4_Naïve_B_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T4_ieqtlcovs_Naïve_B_cell.txt',  sep = "\t", row.names = T, quote = F)




# Classical Monocyte
T0_Classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_Classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Classical_Monocyte")
T1_Classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_Classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Classical_Monocyte")
T2_Classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_Classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Classical_Monocyte")
T3_Classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_Classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Classical_Monocyte")
T4_Classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_Classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Classical_Monocyte")


read.table('/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T0_Classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T0_ieqtlcovs_Classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T1_Classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T1_ieqtlcovs_Classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T2_Classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T2_ieqtlcovs_Classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T3_Classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T3_ieqtlcovs_Classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T4_Classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T4_ieqtlcovs_Classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)



# Memory_CD8+_T_cell
T0_Memory_CD8_T_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Memory_CD8+_T_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Memory_CD8+_T_cell")
T1_Memory_CD8_T_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Memory_CD8+_T_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Memory_CD8+_T_cell")
T2_Memory_CD8_T_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Memory_CD8+_T_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Memory_CD8+_T_cell")
T3_Memory_CD8_T_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Memory_CD8+_T_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Memory_CD8+_T_cell")
T4_Memory_CD8_T_cell_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_Memory_CD8+_T_cell_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Memory_CD8+_T_cell")) %>% as.data.frame() %>% magrittr::set_rownames("Memory_CD8+_T_cell")


read.table('/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T0_Memory_CD8_T_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T0_ieqtlcovs_Memory_CD8+_T_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T1_Memory_CD8_T_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T1_ieqtlcovs_Memory_CD8+_T_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T2_Memory_CD8_T_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T2_ieqtlcovs_Memory_CD8+_T_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T3_Memory_CD8_T_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T3_ieqtlcovs_Memory_CD8+_T_cell.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T4_Memory_CD8_T_cell_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T4_ieqtlcovs_Memory_CD8+_T_cell.txt',  sep = "\t", row.names = T, quote = F)



# B cell
T0_Non_classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T0_Non-classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Non-classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Non-classical_Monocyte")
T1_Non_classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T1_Non-classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Non-classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Non-classical_Monocyte")
T2_Non_classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T2_Non-classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Non-classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Non-classical_Monocyte")
T3_Non_classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T3_Non-classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Non-classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Non-classical_Monocyte")
T4_Non_classical_Monocyte_interaction <- fread('/home/workspace/jogrady/heQTL/work/ieQTL/T4_Non-classical_Monocyte_interaction_input_from_cibersort.txt') %>% pivot_wider(names_from = c("Mixture"), values_from = c("Non-classical_Monocyte")) %>% as.data.frame() %>% magrittr::set_rownames("Non-classical_Monocyte")


read.table('/home/workspace/jogrady/heQTL/data/covariate/T0_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T0_Non_classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T0_ieqtlcovs_Non-classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T1_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T1_Non_classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T1_ieqtlcovs_Non-classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T2_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T2_Non_classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T2_ieqtlcovs_Non-classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T3_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T3_Non_classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T3_ieqtlcovs_Non-classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)
read.table('/home/workspace/jogrady/heQTL/data/covariate/T4_eqtlcovs.txt') %>% as.data.frame() %>% rbind(., T4_Non_classical_Monocyte_interaction) %>% write.table(., '/home/workspace/jogrady/heQTL/data/covariate/T4_ieqtlcovs_Non-classical_Monocyte.txt',  sep = "\t", row.names = T, quote = F)

