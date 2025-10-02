library(MASS)
library(glmnet)
library(ggplot2)
library(broom)
library(dplyr)

# 加载 Boston 房价数据集
data("Boston")
x <- as.matrix(Boston[, -14])
y <- Boston$medv

# 岭回归
fit_ridge <- glmnet(x, y, alpha = 0, nlambda = 100)

# 整理数据
ridge_tidy <- tidy(fit_ridge) %>%
    filter(term != "(Intercept)") %>%
    mutate(log_lambda = log(lambda))

# 绘图
ridge_plot <- ggplot(ridge_tidy, aes(x = log_lambda, y = estimate, color = term)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Ridge Regression Coefficient Paths",
        x = "log(lambda)", y = "Coefficient"
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "Variables"))
ggsave("ridge_coefficients_loglambda.png", plot = ridge_plot, dpi = 300, width = 8, height = 6, bg = "white")

# 选择最优 lambda
set.seed(2025)
cv.ridge <- cv.glmnet(x, y, alpha = 0)
png("ridge_cv_plot.png", width = 800, height = 600)
plot(cv.ridge, xlab = expression(log(lambda)))
dev.off()

# 得到参数
cv.ridge$lambda.min
cv.ridge$lambda.1se
coef(cv.ridge, s = "lambda.min")
coef(cv.ridge, s = "lambda.1se")
