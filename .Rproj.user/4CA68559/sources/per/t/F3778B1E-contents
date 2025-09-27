# ============================
# FATE1 vs steroid phenotype
# ============================

# Packages
library(tidyverse)
library(rstatix)   # shapiro_test, levene_test, t_test, wilcox_test, cohens_d, wilcox_effsize, get_summary_stats
library(ggpubr)    # stat_pvalue_manual

# ---------------------------------------------------------
# 0) Start from your tibble: `dados_filt` (78 x 6 as shown)
#    Required columns:
#      - sample_id (chr)
#      - FATE1_log2TPM (num)
#      - mRNA_K4 (chr)  e.g., "steroid-phenotype-high+proliferation"
# ---------------------------------------------------------




# Pacotes
library(data.table)
library(tidyverse)

# 1) Ler o arquivo (genes x amostras) ---------------------------------------
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-ACC.star_tpm.tsv.gz"
mat_dt <- fread(url)   # 1ª coluna = Ensembl gene ID (com ou sem versão)

# Detecta a coluna de IDs de gene
gene_col <- names(mat_dt)[1]
ens <- mat_dt[[gene_col]]

# Remove a coluna de IDs e vira matriz numérica
mat <- as.matrix(mat_dt[, -1, with = FALSE])
mode(mat) <- "numeric"
rownames(mat) <- ens

# 2) (Opcional) remover versão dos IDs Ensembl (ENSG00000... .xx) ----------
ens_nov <- sub("\\.\\d+$", "", rownames(mat))
rownames(mat) <- ens_nov

# 3) (Opcional) mapear para símbolo (GENCODE/Ensembl) -----------------------
# Dica: para evitar ambiguidade, mantenha Ensembl como chave principal.
# Se quiser símbolo, use biomaRt/AnnotationDbi; exemplo rápido com biomaRt:
# BiocManager::install("biomaRt")
# library(biomaRt)
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# map <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
#              filters="ensembl_gene_id", values=rownames(mat), mart=mart)
# sym <- map$hgnc_symbol[match(rownames(mat), map$ensembl_gene_id)]

# 4) Extrair FATE1 (por Ensembl) -------------------------------------------
# Ensembl de FATE1 (GRCh38): ENSG00000147378
stopifnot("ENSG00000147378" %in% rownames(mat))
fate1 <- mat["ENSG00000147378", ]                # já está em log2(TPM+1)

# 5) Padronizar barcodes e filtrar tumor primário (01) ----------------------
# As colunas aqui já parecem no nível 'sample' (ex.: TCGA-OR-A5K2-01A).
sample_barcode <- colnames(mat)

get_sample_type <- function(bc) substr(bc, 14, 15)   # posições 14-15 (ex.: "01")
sample_type <- get_sample_type(sample_barcode)

fate1_tbl <- tibble(
  sample_id = sample_barcode,              # use assim mesmo (TCGA-XX-XXXX-01A)
  sample_type = sample_type,
  FATE1_log2TPM = as.numeric(fate1)
) %>%
  filter(sample_type == "01") %>%
  select(sample_id, FATE1_log2TPM)

# 6) Integrar com suas métricas TCR/BCR -------------------------------------
# repertoire <- read.csv("repertoire_metrics.csv")  # deve ter 'sample_id' igual aos do Xena
# dat <- repertoire %>% inner_join(fate1_tbl, by = "sample_id")



# Build steroid_phenotype from mRNA_K4 (case-insensitive)
df <- dados_filt %>%
  mutate(
    steroid_phenotype = case_when(
      str_detect(tolower(mRNA_K4), "high") ~ "HSP",
      str_detect(tolower(mRNA_K4), "low")  ~ "LSP",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sample_id, FATE1_log2TPM, steroid_phenotype) %>%
  filter(!is.na(FATE1_log2TPM), !is.na(steroid_phenotype)) %>%
  mutate(steroid_phenotype = factor(steroid_phenotype, levels = c("LSP","HSP")))

# Quick check of group sizes and basic stats
cat("\nGroup counts:\n")
print(count(df, steroid_phenotype))

cat("\nSummary stats per group (mean/SD + median/IQR):\n")
sum_stats <- df %>%
  group_by(steroid_phenotype) %>%
  get_summary_stats(FATE1_log2TPM, type = "full")
print(sum_stats)

# ---------------------------------------------------------
# 1) Assumption checks
# ---------------------------------------------------------
normality <- df %>%
  group_by(steroid_phenotype) %>%
  shapiro_test(FATE1_log2TPM)

homogeneity <- df %>%
  levene_test(FATE1_log2TPM ~ steroid_phenotype)

cat("\nShapiro-Wilk normality (per group):\n")
print(normality)

cat("\nLevene's test for homogeneity of variances:\n")
print(homogeneity)

all_normal <- all(normality$p > 0.05)
equal_var  <- homogeneity$p > 0.05

# ---------------------------------------------------------
# 2) Pick the test + effect size
# ---------------------------------------------------------
if (all_normal && equal_var) {
  chosen <- "Student t-test (equal variances)"
  stat <- df %>% t_test(FATE1_log2TPM ~ steroid_phenotype, var.equal = TRUE) %>%
    add_significance()
  eff  <- df %>% cohens_d(FATE1_log2TPM ~ steroid_phenotype, var.equal = TRUE)
} else if (all_normal && !equal_var) {
  chosen <- "Welch t-test (unequal variances)"
  stat <- df %>% t_test(FATE1_log2TPM ~ steroid_phenotype, var.equal = FALSE) %>%
    add_significance()
  eff  <- df %>% cohens_d(FATE1_log2TPM ~ steroid_phenotype, var.equal = FALSE)
} else {
  chosen <- "Wilcoxon rank-sum test"
  stat <- df %>% wilcox_test(FATE1_log2TPM ~ steroid_phenotype, detailed = TRUE) %>%
    add_significance()
  eff  <- df %>% wilcox_effsize(FATE1_log2TPM ~ steroid_phenotype, ci = TRUE)
  # Optional: Cliff's delta
  # eff_cliff <- effectsize::cliffs_delta(FATE1_log2TPM ~ steroid_phenotype, data = df)
}

cat("\nChosen test:\n"); print(chosen)
cat("\nTest result:\n"); print(stat)
cat("\nEffect size:\n"); print(eff)

# ---------------------------------------------------------
# 3) Plot (violin + box + jitter) with p-value annotation
# ---------------------------------------------------------
# Position of p-value annotation
ypos <- max(df$FATE1_log2TPM, na.rm = TRUE) * 1.05

# Build label for plot
stat_plot <- stat %>%
  mutate(group1 = "LSP", group2 = "HSP",
         y.position = ypos,
         label = paste0("p = ", signif(p, 3)))

p <- ggplot(df, aes(x = steroid_phenotype, y = FATE1_log2TPM, fill = steroid_phenotype)) +
  geom_violin(trim = FALSE, alpha = 0.25, linewidth = 0.3) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.7, size = 2) +
  stat_pvalue_manual(stat_plot, tip.length = 0.01, label = "label") +
  labs(
    title = "FATE1 expression by steroid phenotype (LSP vs HSP)",
    subtitle = paste0("Test: ", chosen),
    x = NULL,
    y = "FATE1 (log2 TPM + 1)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(p)

# (Optional) Save the figure
# ggsave("FATE1_LSP_vs_HSP_boxplot.png", p, width = 6, height = 4, dpi = 300)
# ggsave("FATE1_LSP_vs_HSP_boxplot.pdf", p, width = 6, height = 4)
