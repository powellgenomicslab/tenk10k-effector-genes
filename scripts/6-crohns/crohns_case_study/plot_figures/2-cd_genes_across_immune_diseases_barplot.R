# how many diseases share genes with crohn's sig genes?

library(tidyverse)
tenk_all <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_sig.RDS")
crohns <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS")

df <- tenk_all %>% filter(supercategory == "disease", pheno_cat == "Immune", probeID %in% crohns$probeID) %>% 
  select(pheno_label, phenotype, Gene) %>% 
  distinct() %>% 
  group_by(phenotype, pheno_label) %>% 
  tally() %>% 
  arrange(desc(n))


# 
# df <- tenk_all %>% filter(supercategory == "disease", pheno_cat == "Immune", probeID %in% crohns$probeID) %>% 
#   group_by(Gene) %>% 
#   summarise(unique_diseases = length(unique(phenotype)))
#   
# 
# max(df$unique_diseases)
# 
# df %>% filter(unique_diseases == 3)

# # A tibble: 9 × 2
# phenotype     n
# <chr>     <int>
#   1 crohns     1791
# 2 ibd        1508
# 3 uc          405
# 4 sle         388
# 5 t1dm        343
# 6 ra          317
# 7 psoriasis   300
# 8 ms          224
# 9 eczema      135

#  df
# # A tibble: 9 x 2
#   phenotype     n
#   <chr>     <int>
# 1 crohns     2014
# 2 ibd        1707
# 3 ra         1069
# 4 t1dm       1050
# 5 psoriasis   935
# 6 uc          832
# 7 sle         755
# 8 ms          325
# 9 eczema      268


cd_shared_ds <- ggplot(df, aes(y=fct_reorder(pheno_label, n), x=n)) +
  #geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_bar(stat="identity") +
  labs(x = "Number of MR + Sensitivity + Coloc Genes", y = NULL) +
  #scale_x_continuous(limits = c(0, 200), breaks = c(0, 50, 100, 150, 180)) +
  scale_x_continuous(limits = c(0, 180), breaks = c(0, 50, 100, 150, 180), expand = expansion(c(0, 0.01))) +  #geom_segment( aes(xend=phenotype, yend=0)) +
  #geom_point( size=4, color="orange") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(face = "bold", size = 8, vjust = 1, hjust =0.5),
        axis.text.y = element_text(face = "bold", size = 8, vjust = 1, hjust =1))

cd_shared_ds
# df2 <- tenk_all %>% filter(supercategory == "disease", probeID %in% crohns$probeID) %>% group_by(phenotype) %>% tally() %>% arrange(desc(n))

saveRDS(cd_shared_ds, "resources/crohns_case_study/figures/shared_with_other_immuneds_barplot.RDS")

ggsave("resources/crohns_case_study/figures/immune_ds_barplot.png", cd_shared_ds, device = ragg::agg_png(),
       width = 5.0, height = 2, bg = "white", scaling = 1.0, dpi = 300)
