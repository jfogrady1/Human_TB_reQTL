# VEP calculation

library(vcfR)
library(data.table)
library(ggplot2)
library(tidyverse)
install.packages("ggsankey")
library(networkD3)
# install.packages("remotes")

library(ggsankey)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("dplyr")
library(dplyr) # Also needed

data = read.table("/home/workspace/jogrady/heQTL/software/ensembl-vep/VEP_human_analysis.txt")
vcf = read.vcfR("/home/workspace/jogrady/heQTL/work/DNA_seq/imputation/final/DV_IMPUTED_R0.6_filtered.vcf.gz")

variant_info = data.frame(vcf@fix)

head(variant_info)
variant_info$Imputed = str_split_i(variant_info$INFO, ";",i = 3)
variant_info$Imputed = if_else(variant_info$Imputed == "IMP", "IMPUTED", "CALLED")
variant_info = variant_info %>% select(ID, Imputed)
colnames(variant_info)[1] <- "Uploaded_variation"

table(variant_info$Imputed) 

# CALLED IMPUTED 
#428027 1078921 

colnames(data) <- c("Uploaded_variation", "SYMBOL","Feature", "Feture_type","Consequence","cDNA_position","CANONICAL")
head(data)
data$Consequence <- gsub(",.*", "", data$Consequence)

data <- left_join(data, variant_info)
data <-  data %>% mutate(New_ID = case_when(Consequence == "splice_donor_5th_base_variant" ~ "Splice variant",
                                       Consequence == "splice_donor_region_variant" ~ "Splice variant",
                                       Consequence == "splice_donor_variant" ~ "Splice variant",
                                       Consequence == "splice_region_variant" ~ "Splice variant",
                                       Consequence == "splice_acceptor_variant" ~ "Splice variant",
                                       Consequence == "splice_polypyrimidine_tract_variant" ~ "Splice variant",
                                       Consequence == "start_lost" ~ "Start/Stop lost/gained variant",
                                       Consequence == "stop_lost" ~ "Start/Stop lost/gained variant",
                                       Consequence == "stop_gained" ~ "Start/Stop lost/gained variant",
                                       Consequence == "stop_retained_variant" ~ "Start/Stop lost/gained variant",
                                       Consequence == "3_prime_UTR_variant" ~ "5'/3' UTR variant",
                                       Consequence == "5_prime_UTR_variant" ~ "5'/3' UTR variant",
                                       Consequence == "upstream_gene_variant" ~ "Upstream/Downstream gene variant",
                                       Consequence == "downstream_gene_variant" ~ "Upstream/Downstream gene variant",
                                       Consequence == "inframe_deletion" ~ "Inframe indel",
                                       Consequence == "inframe_insertion" ~ "Inframe indel",
                                       Consequence == "mature_miRNA_variant" ~ "Non-coding transcript/miRNA variant",
                                       Consequence == "non_coding_transcript_exon_variant" ~ "Non-coding transcript/miRNA variant",
                                       Consequence == "frameshift_variant" ~ "Frameshift variant",
                                       Consequence == "intergenic_variant" ~ "Intergenic variant",
                                       Consequence == "intron_variant" ~ "Intronic variant",
                                       Consequence == "missense_variant" ~ "Missense/Synonymous variant",
                                       Consequence == "synonymous_variant" ~ "Missense/Synonymous variant",
                                       .default = as.character(Consequence)))

test = data %>%
  group_by(Imputed, New_ID) %>%
  summarise(n = n())  %>%
  mutate(freq = n / sum(n), .groups = "drop") %>% arrange(desc(freq)) %>% mutate(Levs = factor(New_ID),
                                                               Levs = fct_reorder(New_ID, freq, .desc = TRUE))
test$Levs
View(test)

test$New_ID <- factor(test$New_ID)
test$New_ID
library(RColorBrewer)
ggplot(test, aes(x = Imputed, y = (freq * 100), fill = Levs, label = round(freq * 100, 2))) +
  geom_bar(stat = "identity",width = 0.7) +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) + theme_classic() + 
  scale_fill_manual(name = NULL, values = c(brewer.pal(8, "Dark2"), "darkblue", "darkred")) +
  scale_x_discrete(breaks = c("CALLED", "IMPUTED"),
    labels = c("Called variants\n(n = 428,027)", "Imputed variants\n(n = 1,078,921)")) +
  labs(x = NULL, y = "Frequency of functional consequence (%)") +
    theme(axis.text.x = element_text(size = 15, colour = "black"),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, colour = "black", face = "bold"))
ggsave("/home/workspace/jogrady/heQTL/results/eQTL/results/Variant_location_called_V_Imputed.pdf", width = 12, height = 15, dpi = 600)



head(test)



df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb)

View(df)

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Car features")
data datadata <- data %>% filter(CANONICAL == "YES" | Consequence == "intergenic_variant")


data$Consequence <- gsub(",.*", "", data$Consequence)


library(networkD3)
