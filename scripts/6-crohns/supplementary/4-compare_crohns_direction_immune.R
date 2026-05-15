library(data.table)
library(tidyverse)

source("scripts/util/helper.R")
df_msmsr <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_sig.RDS") %>% ungroup()

crohns_genes <- df_msmr %>% 
  filter(phenotype == "crohns") %>% 
  select(probeID, cell_type)


df_crohns_immune <- df_msmr %>% 
  inner_join(crohns_genes, by = c("probeID", "cell_type")) %>%
  filter(pheno_cat == "Immune") %>%
  mutate(signed_minlog10p = -log10(p_SMR_multi) * sign(b_SMR),
         dir = ifelse(b_SMR < 0, "Negative", "Positive"))

phenos <- c("crohns", "uc", "ibd", "ra", "t1dm", "sle", "eczema", "psoriasis", "ms")

pheno_label <- df_trait_map %>%
  filter(trait_id %in% phenos) %>%
  arrange(match(trait_id, phenos)) %>% 
  select(trait_id, label) %>%
  deframe()

prefix <- "signed -log10 P MR for"
pheno_label_prefix <- paste(prefix, pheno_label[phenos])

tbl_crohns_immune <- df_crohns_immune %>%  
  pivot_wider(id_cols = c(cell_type, Gene),
            names_from = "phenotype",
            values_from = signed_minlog10p) %>% 
  mutate(count_pos = rowSums(across(all_of(phenos)) > 0, na.rm = TRUE),
         count_neg = rowSums(across(all_of(phenos)) < 0, na.rm = TRUE))

write_gs(tbl_crohns_immune, "crohns_immune")
# df_crohns_immune %>% 
#   filter(Gene %in% c("NCF4", "ZBTB38")) %>%
#   ggplot(aes(x = Gene, y = cell_type, fill = dir)) +
#   facet_grid(rows = vars(pheno_label),
#              scales = "free", space = "fixed") +
#   geom_tile(color = "black") +
#   coord_fixed() +
#   theme_bw() +
#   geom_bar()
#   paletteer::scale_fill_paletteer_c(
#     "ggthemes::Red-Blue-White Diverging",
#     na.value = "grey90", direction = -1,
#     guide = guide_colorbar(theme = theme(legend.key.width = unit(7.5, "lines"),
#                                          legend.key.height = unit(0.75, "lines")))) +
#   theme(legend.position = "bottom")

############################################################################################################################################################################################################

# Not in publication but for referece. 

library(tidyverse)
tenk_all <- readRDS("resources/crohns_case_study/postprocess/tenk_alltraits_sig.RDS")
crohns <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS")

df <- tenk_all %>% filter(supercategory == "disease", pheno_cat == "Immune", probeID %in% crohns$probeID) %>% group_by(phenotype) %>% tally() %>% arrange(desc(n))

df 
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

# bar plot 
df %>%
  mutate(phenotype = fct_reorder(phenotype, desc(n))) %>%
  ggplot( aes(x=phenotype, y=n)) +
  # geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  geom_segment( aes(xend=phenotype, yend=0)) +
  geom_point( size=4, color="orange") +
  theme_bw()


