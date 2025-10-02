library(MASS)
library(ncvreg)
library(ggplot2)
library(reshape2)

# 加载 Boston 数据
data("Boston")
x <- as.matrix(Boston[, -14])
y <- Boston$medv

# SCAD 回归拟合
fit_SCAD <- ncvreg(x, y, family = "gaussian", penalty = "SCAD", nlambda = 100)

# 提取系数矩阵
beta <- fit_SCAD$beta[-1, ]
lambda <- fit_SCAD$lambda
log_lambda <- log(lambda)

# 构建数据框用于绘图
coef_df <- as.data.frame(beta)
coef_df$Variable <- rownames(beta)
coef_df <- melt(coef_df, id.vars = "Variable", variable.name = "Step", value.name = "Coefficient")

coef_df$log_lambda <- rep(log_lambda, each = nrow(beta))

# 绘图
SCAD_plot <- ggplot(coef_df, aes(x = log_lambda, y = Coefficient, color = Variable)) +
    geom_line(linewidth = 1) +
    labs(
        title = "SCAD Regression Coefficient Paths",
        x = "log(lambda)", y = "Coefficient"
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "Variables"))
ggsave("SCAD_coefficients_loglambda.png", plot = SCAD_plot, dpi = 300, width = 8, height = 6, bg = "white")

set.seed(2025)
cv.SCAD <- cv.ncvreg(x, y, family = "gaussian", penalty = "SCAD")
png("SCAD_cv_plot.png", width = 800, height = 600)
plot(cv.SCAD, xlab = expression(log(lambda)))
dev.off()
