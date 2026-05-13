library(tidyverse)
library(patchwork)
library(dplyr)
library(forcats)

# Get DEG results for annotation bar plot 
deg_genes <- readRDS("resources/crohns_case_study/deg/mr_max_evidence_and_deg_innerjoin_by_majct.RDS") %>% 
  select(Gene, major_cell_type, deg_maj_ct) %>% 
  mutate(across(where(is.logical), function(x) replace_na(x, FALSE)))

# Left join to tenk results - NOTE need to ensure major cell type annotation is updated 
# Create a no annotation col 
tenk_sig <- readRDS("resources/crohns_case_study/postprocess/tenk_crohns_sig.RDS") %>% 
  ungroup() %>% 
  mutate(across(where(is.logical), function(x) replace_na(x, FALSE))) %>% 
  mutate(p_transform = -log10(p_SMR_multi)*sign(b_SMR)) %>% 
  mutate(major_cell_type = as.character(major_cell_type) %>% dplyr::replace_when(cell_type == "ILC" ~ "NK", cell_type == "Plasmablast" ~ "Plasma B")) %>% 
  left_join(deg_genes, by = join_by(Gene, major_cell_type)) %>%
  mutate(no_cd_known = !cd_known) %>% 
  mutate(across(where(is.logical), function(x) replace_na(x, FALSE))) 

# Each bar plot should show the proportion of significant genes that are annotated by any of those categories
annot_cols <- c("eqtlgen_sig", "sys_review", "open_targets_cd", "liuetal", "drug_target", "drug_pathway_tnf", 
                "drug_pathway_il23", "drug_pathway_integrin", "cd_known", "deg_maj_ct", "no_cd_known")

# Rename the columns for figure, add star for 'non literature'
rename_cols <- c("Gene", "*Bulk eQTLGen MR", "Systematic Review", "95th Pctl Open Targets", "IBD metaGWAS Nearest Gene", "Drug Target", "TNF Pathway", 
                "IL23 Pathway", "Integrin Pathway", "Any Literature Annotation", "*DEG in Matched Major Cell Type", "No Literature Annotation")

# Attempt 2 
plot_data <- tenk_sig %>% 
  select(Gene, all_of(annot_cols)) %>% 
  distinct() %>% 
  rename_with(~ rename_cols) %>% 
  pivot_longer(cols = where(is.logical), names_to = "annotation", values_to = "status") %>% 
  group_by(annotation) %>% 
  summarise(total = sum(status), 
            prop = round(sum(status)/n_distinct(Gene)*100))

annot_df <- readRDS("resources/crohns_case_study/crohns_relevant_gene_lists/external_source_gene_annot_df.RDS")

annot_df <- annot_df %>% 
  #mutate(`No Literature Annotation` = !cd_known) %>% 
  pivot_longer(cols = where(is.logical), names_to = "annotation", values_to = "status") %>% 
  mutate(annotation = case_when(
    annotation == "sys_review" ~ "Systematic Review",
    annotation == "open_targets_cd"  ~ "95th Pctl Open Targets",
    annotation == "liuetal" ~ "IBD metaGWAS Nearest Gene",
    annotation == "drug_target" ~ "Drug Target",
    annotation == "drug_pathway_tnf" ~ "TNF Pathway",
    annotation == "drug_pathway_il23" ~ "IL23 Pathway",
    annotation == "drug_pathway_integrin" ~ "Integrin Pathway",
    annotation == "cd_known" ~ "Any Literature Annotation")) %>% 
  group_by(annotation) %>% 
  summarise(total_lit = sum(status), 
            prop = round(sum(status)/n_distinct(Gene)*100))


plot_data2 <- plot_data %>% 
  left_join(annot_df %>% select(-prop), by = "annotation")

# Fill up total lit with values 
deg <- readRDS("resources/crohns_case_study/deg/crohns_deg_pre-processed_revision.RDS")
deg <- deg %>% filter(major_cell_type %in% tenk_sig$major_cell_type)
length(unique(deg$Gene)) # 5667

bulkeqtlgen <- readRDS("resources/crohns_case_study/postprocess/eqtlgen_alltraits_raw.RDS") %>% 
  filter(phenotype == "crohns") %>% 
  filter(sig == TRUE) 
length(unique(bulkeqtlgen$Gene)) # 578

plot_data2[plot_data2$annotation == "*Bulk eQTLGen MR", ]$total_lit = length(unique(bulkeqtlgen$Gene))
plot_data2[plot_data2$annotation == "*DEG in Matched Major Cell Type", ]$total_lit = length(unique(deg$Gene))
plot_data2[plot_data2$annotation == "No Literature Annotation", ]$total_lit = as.numeric(plot_data2[plot_data2$annotation == "Any Literature Annotation",]$total_lit)

# Make factor 
plot_data2$annotation <- factor(plot_data2$annotation, levels = c("Any Literature Annotation", "No Literature Annotation", "Systematic Review", "95th Pctl Open Targets", "IBD metaGWAS Nearest Gene", "Drug Target", "TNF Pathway", 
                               "IL23 Pathway", "Integrin Pathway", "*Bulk eQTLGen MR", "*DEG in Matched Major Cell Type"))

#plot_data2$annotation_new <- paste0(plot_data2$annotation, " (", plot_data2$total, " out of ", plot_data2$total_lit, ")")
plot_data2$annotation_new <- paste0(plot_data2$annotation, " (", plot_data2$total_lit, ")")
plot_data2 <- plot_data2 %>%  arrange(annotation)
plot_data2$annotation_new <- factor(plot_data2$annotation_new, levels = plot_data2$annotation_new)

# New plot 
p1 <- ggplot(plot_data2 %>% arrange(annotation), aes(y = annotation_new, x = total)) +
  geom_col() +
  geom_text(size = 2, aes(label = paste0(total)), hjust = -0.3, vjust = 0.6, colour = "black", fontface = "bold") +
  #coord_flip() +
  #facet_wrap(~higher_annotation, scales = "free", space = "free_x") + 
  #scale_fill_manual(values = c("TRUE" = "#4682B4", "FALSE" = "#D3D3D3"), guide = guide_legend(reverse = TRUE)) +
  #scale_x_continuous(limits = c(0, 200), breaks = c(0, 50, 100, 150, 180)) +
  scale_x_continuous(limits = c(0, 180), breaks = c(0, 50, 100, 150, 180), expand = expansion(c(0, 0.01))) +  #geom_segment( aes(xend=phenotype, yend=0)) +
  scale_y_discrete(limits = rev) +
  theme_minimal() + 
  labs(x = "Number of MR + Sensitivity + Coloc Genes", y = NULL) +
  theme(axis.title.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(face = "bold", size = 8, vjust = 1, hjust =0.5),
    axis.text.y = element_text(face = "bold", size = 8, vjust = 1, hjust =1))
p1

saveRDS(p1, "resources/crohns_case_study/figures/annotation_fig.RDS")

ggsave("resources/crohns_case_study/figures/ext_source_annotation_barplot.png", p1, device = ragg::agg_png(),
       width = 5.0, height = 2, bg = "white", scaling = 1.0, dpi = 300)

