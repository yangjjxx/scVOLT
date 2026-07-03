library(data.table)
library(uwot)
library(scales)

# 读取 gene × cell 表达矩阵
# 第一列是 gene 名称，后面的列是 cell
expr_df <- read.csv("../data/02_time_series_analysis/GSE114461_Bortezomib.csv.gz", header = TRUE)

gene_names <- rownames(expr_df)

expr_mat <- as.matrix(expr_df)
mode(expr_mat) <- "numeric"

# gene × cell 转成 cell × gene
expr_cell_gene <- t(expr_mat)

cell_names <- rownames(expr_cell_gene)

# 读取预测值
drug_score <- read.csv("../results/02_time_series_analysis/Bortezomib_score.csv", header = FALSE)[[1]]

# 检查长度是否和 cell 数一致
stopifnot(length(drug_score) == nrow(expr_cell_gene))


umap_res <- uwot::umap(
  expr_cell_gene,
  n_neighbors = 40,
  min_dist = 1.5,
  n_components = 2,
  metric = "cosine"
)

#构建umap_df
umap_df <- data.frame(
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2],
  cell = cell_names,
  drug_score = drug_score,
  stringsAsFactors = FALSE
)

# 从 cell 名称中提取时间
umap_df$time <- ifelse(grepl("t0", umap_df$cell), "0",
                       ifelse(grepl("t12", umap_df$cell), "12",
                              ifelse(grepl("t48", umap_df$cell), "48",
                                     ifelse(grepl("t96", umap_df$cell), "96", NA))))

umap_df$time <- factor(umap_df$time, levels = c("0", "12", "48", "96"))


# 计算 drug_score 百分位数
score <- umap_df$drug_score

umap_df$drug_percentile <- 100 * 
  (rank(score, ties.method = "average", na.last = "keep") - 1) /
  (sum(!is.na(score)) - 1)

# 查看结果
head(umap_df)
