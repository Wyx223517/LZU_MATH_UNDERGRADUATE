library(GGally)
library(HSAUR2)
library(tidyr)
library(dplyr)
library(forcats)

# 设置保存路径
save_dir <- "E:/719465160/课程作业/多元统计分析/HW04/聚类分析/Figure"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

data("pottery", package = "HSAUR2")

## 散点图矩阵
p1 <- ggscatmat(
    data = pottery,
    columns = 1:9, color = "kiln"
)
ggsave(file.path(save_dir, "scatter_matrix.png"), plot = p1, width = 12, height = 10, dpi = 300)

## 各个变量的方差
var_tbl <- pottery |>
    summarise(across(
        -kiln, var
    )) |>
    pivot_longer(everything(),
        names_to = "Chemical",
        values_to = "Var"
    )
knitr::kable(var_tbl, format = "latex", booktabs = TRUE, digits = 4)

## 标准化处理
pottery_scale <- pottery |>
    mutate(across(-kiln, \(x) (x - mean(x)) / sd(x)))

## 计算主成分
pot_pca <- princomp(pottery[, 1:9], cor = TRUE, scores = TRUE) |>
    predict()

# PCA 按 kiln 分类绘图并保存
as_tibble(pot_pca) |>
    mutate(kiln = pottery[["kiln"]]) |>
    ggplot(aes(
        x = Comp.1, y = Comp.2, color = kiln
    )) +
    geom_point() -> p2
ggsave(file.path(save_dir, "PCA_by_kiln.png"), plot = p2, width = 8, height = 6, dpi = 300)

## 对不同的类数计算解释百分比
elbow_plot <- function(
    data,
    max_k = min(ncol(data), 10)) {
    rat <- numeric(max_k)
    rat[1] <- 0
    for (k in 2:max_k) {
        rkm <- kmeans(
            data,
            centers = k,
            nstart = 20
        )
        rat[k] <- (1 - rkm$tot.withinss / rkm$totss) * 100
    }
    plot(rat,
        xlab = "类数", ylab = "离差平方和解释百分比",
        type = "b"
    )
}
# 保存 Elbow Plot
png(file.path(save_dir, "elbow_plot.png"), width = 800, height = 600)
elbow_plot(pottery_scale[, 1:9], max_k = 6)
dev.off()

## 作k均值聚类
pot_km <- kmeans(pottery_scale[, 1:9],
    centers = 3, nstart = 20
)
pot_km

## 按各类中心的变量值对类别重新编号
cluster_reorder <- function(
    data, km_res,
    label = rownames(data)) {
    xrmean <- rowMeans(data)
    clus_old <- km_res$cluster
    names(clus_old) <- label
    clus_new <- clus_old
    clus_new[] <- as.integer(fct_reorder(
        factor(clus_old), xrmean, mean,
        .desc = TRUE
    ))
    clus_new
}
clus_new <- cluster_reorder(
    pottery_scale[, 1:9],
    pot_km,
    label = pottery_scale[["kiln"]]
)

## 对比kiln和clus_new的值：
table(clus_new, pottery_scale[["kiln"]])

## 按分类作主成分图：
as_tibble(pot_pca) |>
    mutate(cluster = factor(clus_new)) |>
    ggplot(aes(
        x = Comp.1, y = Comp.2, color = cluster
    )) +
    geom_point() -> p3
ggsave(file.path(save_dir, "PCA_by_cluster.png"), plot = p3, width = 8, height = 6, dpi = 300)
