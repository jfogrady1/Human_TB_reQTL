# load tidyverse package
library(tidyverse)
library(dplyr)
library(ggforce)


args = commandArgs(trailingOnly=TRUE)


# read in data
#eigen vectors
pca <- read.table(args[1], sep = " ", header = F)

# eigenval
eigenval <- scan(args[2])
eigenval

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "Ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca$Ind <- paste0("P_", pca$Ind)

Ancestry = read.table(args[3], sep = "\t", header = T)
pca <- left_join(pca, Ancestry)
colnames(Ancestry)
colnames(pca)
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + ylab("Percentage (%) variance explained") + theme_bw() + coord_equal() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

#Define color and shape values for the plot
color_values <- c("#ca0020", "#0571b0", "#7b3294", "#f4a582")
shape_values <- 0:12

# Create the ggplot object
b <- ggplot(pca, aes(PC1, PC2, col = Main_Ancestry, shape = Sub_Ancestry)) +
  geom_point(size = 4) +
  scale_colour_manual(values = color_values) +
  scale_shape_manual(values = shape_values) +
  theme_bw() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  guides(
    col = guide_legend(title = "Main ethnicity", size = 4),
    shape = guide_legend(title = "Sub ethnicity")
  ) +
  theme(axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  coord_equal() +
  stat_ellipse(level = .95, linetype = 1)

        # Display the plot
cowplot::plot_grid(b, a, rel_widths = c(1.5, 1), nrow = 1)
ggsave(args[4], device = "pdf", width = 12, height = 8, dpi = 600)
