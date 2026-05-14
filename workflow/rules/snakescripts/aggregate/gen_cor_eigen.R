library(tidyverse)

INPUT <- snakemake@input
OUTPUT <- snakemake@output

# INPUT <- list("results/aggregate/tenk10k_phase1.gen_cor.ldak.tsv")
# OUTPUT <- list("results/aggregate/tenk10k_phase1.gen_cor_eigen.tsv")

df_all <- read_tsv(INPUT[[1]])


# Remove traits with NA value
# neversmokers has a negative heritability estimate -> not sure why?
exclude_traits <- c("pancreas", "neversmokers")

# get heritability estimates, taking maximum
df_her <- df_all %>% 
    filter(Component == "Her1_All") %>% 
    distinct(trait1, Value) %>% 
    group_by(trait1) %>%
    slice_max(Value)

filter_criteria <- expression(
    !trait1 %in% exclude_traits & 
    !trait2 %in% exclude_traits & 
    Component == "Cor_All" &
    !is.na(Value)
)

# create a square symettrical matrix of correlations
df1 <- df_all %>% 
    filter(eval(filter_criteria))  %>% 
    select(trait1, trait2, value = Value)

df2 <- df_all %>% 
    filter(eval(filter_criteria))  %>% 
    select(trait1 = trait2, trait2 = trait1, value = Value)

traits <- unique(c(df1$trait1, df1$trait2))
mat_cor  <- bind_rows(df1, df2) %>% 
    pivot_wider(names_from = trait2, values_fill = 1) %>% 
    column_to_rownames("trait1") %>% 
    data.matrix() %>% 
    .[traits, traits]

# eigenv <- eigen(mat_cor)$values

# Li and Ji (2005) method
calculate_meff <- function(cor_mat) {
  # Compute eigenvalues of the correlation matrix
  eigenvalues <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  # Apply the function: I(x ≥ 1) + (x - floor(x)), for x > 0
  f <- function(x) {
    if (x <= 0) return(0)
    indicator <- as.numeric(x >= 1)
    fractional <- x - floor(x)
    return(indicator + fractional)
  }
  # Sum over all eigenvalues
  meff <- sum(sapply(eigenvalues, f))
  return(meff)
}

n_effective_test <- calculate_meff(mat_cor)

# Perform hierarchical clustering on the signed correlation matrix
dist_signed_cor <- as.dist(1 - abs(mat_cor))  # Convert signed correlation matrix to distance matrix
hc_signed <- hclust(dist_signed_cor, method = "ward.D2")  # Hierarchical clustering using Ward's method

# Cut the dendrogram to determine the number of independent clusters
num_clusters <- length(unique(cutree(hc_signed, h = 0.5)))  # Adjust 'h' as needed for your threshold

# Save dendrogram plot for signed correlation matrix
pdf(sub(".pdf", "_signed.pdf", OUTPUT[[1]]))
plot(hc_signed, main = "Hierarchical Clustering Dendrogram (Signed Correlation)", xlab = "Traits", sub = "", cex = 0.6)
dev.off()