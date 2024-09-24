# Script to check allele concordance of sites called common to merged_HG001 timepoints
# needed to make sure that deep variant is calling alleles correctly
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
#args[1] <- "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/DV_merged_raw_calls.txt"
#args[2] <- "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/HG001_calls.txt"
#args[3] <- "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/HG005_calls.txt"
#args[4] <- "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/HG006_calls.txt"
#args[5] <- "/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/HG007_calls.txt"

merged <- read.table(args[1], sep ="\t", header = F)
HG001 <- read.table(args[2], sep = "\t", header = F)
HG005 <- read.table(args[3], sep = "\t", header = F)
HG006 <- read.table(args[4], sep = "\t", header = F)
HG007 <- read.table(args[5], sep = "\t", header = F)


head(HG001)
merged <-merged %>% unite(., POS, c(V1, V2), sep = "_")
HG001 <- HG001 %>% unite(., POS, c(V1, V2), sep = "_")
HG005 <- HG005 %>% unite(., POS, c(V1, V2), sep = "_")
HG006 <- HG006 %>% unite(., POS, c(V1, V2), sep = "_")
HG007 <- HG007 %>% unite(., POS, c(V1, V2), sep = "_")



colnames(merged) <- c("POS", "mergedREF", "mergedALT")
colnames(HG001) <- c("POS", "HG001REF", "HG001ALT")
colnames(HG005) <- c("POS", "HG005REF", "HG005ALT")
colnames(HG006) <- c("POS", "HG006REF", "HG006ALT")
colnames(HG007) <- c("POS", "HG007REF", "HG007ALT")

merged_HG001<- inner_join(merged, HG001)
merged_HG005<- inner_join(merged, HG005)
merged_HG006<- inner_join(merged, HG006)
merged_HG007<- inner_join(merged, HG007)

merged_HG001$MatchALT <- ifelse(merged_HG001$mergedALT == merged_HG001$HG001ALT & merged_HG001$mergedREF == merged_HG001$HG001REF, TRUE, FALSE)
merged_HG005$MatchALT <- ifelse(merged_HG005$mergedALT == merged_HG005$HG005ALT & merged_HG005$mergedREF == merged_HG005$HG005REF, TRUE, FALSE)
merged_HG006$MatchALT <- ifelse(merged_HG006$mergedALT == merged_HG006$HG006ALT & merged_HG006$mergedREF == merged_HG006$HG006REF, TRUE, FALSE)
merged_HG007$MatchALT <- ifelse(merged_HG007$mergedALT == merged_HG007$HG007ALT & merged_HG007$mergedREF == merged_HG007$HG007REF, TRUE, FALSE)


dim(merged_HG001)
dim(merged_HG005)
dim(merged_HG006)
dim(merged_HG007)


table(merged_HG001$MatchALT)
table(merged_HG005$MatchALT)
table(merged_HG006$MatchALT)
table(merged_HG007$MatchALT)


merged_HG001_spurious <- merged_HG001 %>% filter(MatchALT == FALSE)
merged_HG005_spurious <- merged_HG005 %>% filter(MatchALT == FALSE)
merged_HG006_spurious <- merged_HG006 %>% filter(MatchALT == FALSE)
merged_HG007_spurious <- merged_HG007 %>% filter(MatchALT == FALSE)

dim(HG001_spurious_triallelic)
HG001_spurious_triallelic <- merged_HG001_spurious %>% filter(str_detect(HG001ALT, ",") | str_detect(mergedALT, ","))
HG005_spurious_triallelic <- merged_HG005_spurious %>% filter(str_detect(HG005ALT, ",") | str_detect(mergedALT, ","))
HG006_spurious_triallelic <- merged_HG006_spurious %>% filter(str_detect(HG006ALT, ",") | str_detect(mergedALT, ","))
HG007_spurious_triallelic <- merged_HG007_spurious %>% filter(str_detect(HG007ALT, ",") | str_detect(mergedALT, ","))



dim(HG001_spurious_triallelic)

dim(HG005_spurious_triallelic)

dim(HG006_spurious_triallelic)

dim(HG007_spurious_triallelic)

proportions_1 <- table(merged_HG001$MatchALT) / length(merged_HG001$mergedALT)
proportions_1 <- as.data.frame(proportions_1)
colnames(proportions_1)[1] <- "category"
proportions_1$Group <- "HG001"

proportions_5 <- table(merged_HG005$MatchALT) / length(merged_HG005$mergedALT)
proportions_5 <- as.data.frame(proportions_5)
colnames(proportions_5)[1] <- "category"
proportions_5$Group <- "HG005"

proportions_6 <- table(merged_HG006$MatchALT) / length(merged_HG006$mergedALT)
proportions_6 <- as.data.frame(proportions_6)
colnames(proportions_6)[1] <- "category"
proportions_6$Group <- "HG006"

proportions_7 <- table(merged_HG007$MatchALT) / length(merged_HG007$mergedALT)
proportions_7 <- as.data.frame(proportions_7)
colnames(proportions_7)[1] <- "category"
proportions_7$Group <- "HG007"

proportions <- rbind(proportions_1, proportions_5, proportions_6, proportions_7)

proportions <- proportions %>%
  mutate(Percentage = Freq * 100)

meantrue <- proportions %>% filter(category==TRUE)
meantrue <- mean(meantrue$Percentage)

meanfalse <- proportions %>% filter(category==FALSE)
meanfalse <- mean(meanfalse$Percentage)
meanfalse

ggplot(proportions, aes(x = category, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.8) +
  scale_fill_brewer(palette="Paired") + theme_bw() +
  labs(y = "Percentage (%)",
       x = "Overlapping variants") +
  geom_segment(y = meanfalse, x = 0, yend = meanfalse, xend = 1.55, linetype = "dashed", size = 1)+
  geom_segment(y = meantrue, x = 0, yend = meantrue, xend = 3.1, linetype = "dashed", size = 1) +
  theme(legend.position = "right",
        title = element_text(size = 16),
        legend.background = element_rect(color = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size= 14, colour = "black"),
        legend.key.size = unit(2, "cm"),
        legend.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14)) +
  scale_x_discrete(labels=c("FALSE" = "Discordant", "TRUE" = "Concordant")) +
  annotate("text", x = 1, y = 5, size = 8, label = paste(round(meanfalse,3),"%")) +
  annotate("text", x = 2, y = 100, size = 8, label = paste(round(meantrue,3),"%"))
  
  

ggsave("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/DV_GIAB_concordance.pdf", width = 12, height = 12, dpi = 600)

#ggsave(filename = args[6], width = 15, height = 12, dpi = 600, device = "pdf")

# Remove triallelic sites

merged_HG001 <- merged_HG001 %>% filter(!(POS %in% HG001_spurious_triallelic$POS))
merged_HG005 <- merged_HG005 %>% filter(!(POS %in% HG005_spurious_triallelic$POS))
merged_HG006 <- merged_HG006 %>% filter(!(POS %in% HG006_spurious_triallelic$POS))
merged_HG007 <- merged_HG007 %>% filter(!(POS %in% HG007_spurious_triallelic$POS))





dim(merged_HG001)

proportions_1 <- table(merged_HG001$MatchALT) / length(merged_HG001$mergedALT)
proportions_1 <- as.data.frame(proportions_1)
colnames(proportions_1)[1] <- "category"
proportions_1$Group <- "HG001"

proportions_5 <- table(merged_HG005$MatchALT) / length(merged_HG005$mergedALT)
proportions_5 <- as.data.frame(proportions_5)
colnames(proportions_5)[1] <- "category"
proportions_5$Group <- "HG005"

proportions_6 <- table(merged_HG006$MatchALT) / length(merged_HG006$mergedALT)
proportions_6 <- as.data.frame(proportions_6)
colnames(proportions_6)[1] <- "category"
proportions_6$Group <- "HG006"


proportions_7 <- table(merged_HG007$MatchALT) / length(merged_HG007$mergedALT)
proportions_7 <- as.data.frame(proportions_7)
colnames(proportions_7)[1] <- "category"
proportions_7$Group <- "HG007"

proportions <- rbind(proportions_1, proportions_5, proportions_6, proportions_7)

proportions <- proportions %>%
  mutate(Percentage = Freq * 100)

meantrue <- proportions %>% filter(category==TRUE)
meantrue <- mean(meantrue$Percentage)
proportions
meanfalse <- proportions %>% filter(category==FALSE)
meanfalse <- mean(meanfalse$Percentage)


ggplot(proportions, aes(x = category, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position=position_dodge(), alpha = 0.8) +
  scale_fill_brewer(palette="Paired") + theme_bw() +
  labs(y = "Percentage (%)",
       x = "Overlapping variants") +
  geom_segment(y = meanfalse, x = 0, yend = meanfalse, xend = 1.55, linetype = "dashed", size = 1)+
  geom_segment(y = meantrue, x = 0, yend = meantrue, xend = 3.1, linetype = "dashed", size = 1) +
  theme(legend.position = "right",
        title = element_text(size = 16),
        legend.background = element_rect(color = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size= 14, colour = "black"),
        legend.key.size = unit(2, "cm"),
        legend.title = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 14)) +
  scale_x_discrete(labels=c("FALSE" = "Discordant", "TRUE" = "Concordant")) +
  annotate("text", x = 1, y = 5, size = 8, label = paste(round(meanfalse,3),"%")) +
  annotate("text", x = 2, y = 100, size = 8, label = paste(round(meantrue,3),"%"))


ggsave("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/variantcomparison/DV_no_triallelic_GIAB_concordance.pdf", width = 12, height = 12, dpi = 600)
