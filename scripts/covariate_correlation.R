
# Load in packages
library(data.table)
library(tidyverse)
library(ggplot2)



groupPredict<-function(dataResponse,dataPredictors,R2Type=c("unadjusted","adjusted")){
  # set.seed(1)
  # n<-100
  # dataResponse<-matrix(data=rnorm(n=n*3),nrow=n) #100*3.
  # dataPredictors<-matrix(data=rnorm(n=n*10),nrow=n) #100*10.
  # R2Type<-"unadjusted"

  R2s<-rep(NA,ncol(dataResponse))
  for(j in 1:ncol(dataResponse)){
    # j<-1
    mod<-lm(as.numeric(dataResponse[,j])~as.matrix(dataPredictors))
    # summary(mod)
    if(R2Type=="unadjusted"){
      R2s[j]<-summary(mod)$r.squared
    }else{
      R2s[j]<-summary(mod)$adj.r.squared
    }
  }

  return(R2s)
}

# Function to order by sampleID
custom_sort <- function(s) {
  parts <- strsplit(s, "_")[[1]]
  as.integer(parts[2])
}

calc_R2_individual <- function(sample) {
  # Read in the counts
  # Extact relevant columns (sample expression)
  bed <- fread(paste0('/home/workspace/jogrady/heQTL/work/RNA_seq/quantification/TMM_INT_counts_',sample,'.bed.gz'))
  expr <- t(bed[,-(1:4)]) # CHR , START, END, GeneID - all removed before transposing


  # This is PCA path
  # Run PCA
  prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE) #This should take less than a minute.
  PCs<-prcompResult$x #368*368. The columns are PCs.

  # Choose K
  resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)


  K_elbow<-resultRunElbow #12.

  # Load in the known covariates and extract relevant ones that we think are important
  cov <- read.table(paste0('/home/workspace/jogrady/heQTL/data/covariate/',sample,'_metadata.txt'), sep = "\t", header = T) %>% select(-c(Group, subgroup)) # uninformative

  # convert to factor and numeric
  # Should also be sex but not for discussion here
  cov$gender <- factor(cov$gender, levels = c("male", "female"), labels = c(0,1))

  # Specify the levels for the "smear_results" column
  cov$smear_results <- factor(cov$smear_results, levels = c("Positive", "Negative"), labels = c(1,0))

  # Same for subgroup_att
  cov$subgroup_att <- factor(cov$subgroup_att, levels = c("Outbreak TB strain", "Short_ATT", "Long_ATT", "Difficult TB Cases", "TB Drug Resistance"), labels = c(1,2,3,4,5))

  # Same of days from ATT
  cov$days_from_att <- scale(as.numeric(cov$days_from_att), center = TRUE)

  # Age converted to Numeric
  cov$Age <- scale(as.numeric(cov$Age), center = TRUE)
  # Response group
  cov$response_group <- factor(cov$response_group, levels = c('weaker', 'similar', 'stronger_initial', 'stronger_delayed'), labels = c(0,1,2,3))


  # Read in the genotype PCs and join with the covariates - here read in top 10 as it will be interesting for reviewer
  # Join and sort
  pca <- read.table('/home/workspace/jogrady/heQTL/data/covariate/Genotypes.PCA_eigenvect.txt', sep = "\t", header = T) %>% select(1:3) # top 3 genotype PCS
  known <- left_join(pca, cov, by = c("SampleID" = "Patient_ID")) # join
  known_sorted <- known %>% arrange(sapply(SampleID, custom_sort)) # sort so in same order as PC
  rownames(known_sorted) <- known_sorted$SampleID
  known_sorted <- known_sorted %>% select(-1)

  known_sorted <- known_sorted[rownames(expr),] # get in same order to avoid any doubt.

  dim(known_sorted) # check dimensions

  identical(rownames(known_sorted),rownames(expr)) #TRUE is good.


  colnames(known_sorted)[1:3]<-paste0("genotypePC",1:3) 
  known_sorted$response_group <- factor(known_sorted$response_group, labels = c(0,1,2,3)) # make sure this is a factor for later


  # ============================================================================
  # Calculate R2 explained by each PC for each covariate
  # ============================================================================

  PCsTop <- PCs[, 1:K_elbow]

  # Function to calculate R2 for each PC against a covariate
  calculatePC_R2 <- function(covariate_vec, PC_matrix) {
    # covariate_vec: numeric vector or factor (will be converted to numeric)
    # PC_matrix: matrix of PCs (samples x PCs)
    
    # Convert factors to numeric while preserving the underlying values
    if (is.factor(covariate_vec)) {
      covariate_vec <- as.numeric(covariate_vec)
    }
    
    R2_values <- rep(NA, ncol(PC_matrix))
    names(R2_values) <- colnames(PC_matrix)
    
    for (i in 1:ncol(PC_matrix)) {
      # Fit model: covariate ~ PC_i
      mod <- lm(covariate_vec ~ PC_matrix[, i])
      R2_values[i] <- summary(mod)$r.squared
    }
    
    return(R2_values)
  }

  # Calculate R2 for each PC for each covariate (keeping factors as factors)
  PC_R2_list <- list()
  for (cov_name in colnames(known_sorted)) {
    PC_R2_list[[cov_name]] <- calculatePC_R2(known_sorted[[cov_name]], PCsTop)
  }

  # Convert to data frame for plotting
  PC_R2_df <- as.data.frame(PC_R2_list) %>%
    rownames_to_column(var = "PC") %>%
    pivot_longer(cols = -PC, names_to = "Covariate", values_to = "R2") %>% mutate(Timepoint = sample)

  # Print summary
  print("R2 explained by each PC for each covariate:")
  print(PC_R2_df %>% pivot_wider(names_from = PC, values_from = R2) %>% head(10))

  # Create visualization 1: Heatmap of R2 values
  plot_heatmap_r2 <- ggplot(PC_R2_df, aes(x = PC, y = Covariate, fill = R2)) +
    geom_tile() +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "R²") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "R² Explained by Each PC for Each Covariate",
        x = "Principal Component", y = "Covariate")
  return(PC_R2_df)
}

T0 <- calc_R2_individual("T0")
T1 <- calc_R2_individual("T1")
T2 <- calc_R2_individual("T2")
T3 <- calc_R2_individual("T3")
T4 <- calc_R2_individual("T4")

plot_individual <- rbind(T0,T1,T2,T3,T4)

ggplot(plot_individual, aes(x = PC, y = Covariate, fill = R2)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "R² Explained by Each PC for Each Covariate",
      x = "Principal Component", y = "Covariate") +
  facet_wrap(~ Timepoint, scales = "free_x")

# Create visualization 1: Heatmap of R2 values
plot_heatmap_r2 <- ggplot(PC_R2_df, aes(x = PC, y = Covariate, fill = R2)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "R² Explained by Each PC for Each Covariate",
      x = "Principal Component", y = "Covariate")

# Check correlation between known covariates and PCs
ggsave('/home/workspace/jogrady/heQTL/results/Response_1/Figures/individual_correlation_heatmap.pdf', width = 12, height = 12, dpi = 600)
