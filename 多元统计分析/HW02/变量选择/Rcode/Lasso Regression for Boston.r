library(MASS)
library(glmnet)
library(ggplot2)
library(broom)
library(dplyr)

# 加载 Boston 房价数据集
data("Boston")
x <- as.matrix(Boston[, -14])
y <- Boston$medv

# Lasso回归
fit_lasso <- glmnet(x, y, alpha = 1, nlambda = 100)

# 整理数据
lasso_tidy <- tidy(fit_lasso) %>%
    filter(term != "(Intercept)") %>%
    mutate(log_lambda = log(lambda))

# 绘图
lasso_plot <- ggplot(lasso_tidy, aes(x = log_lambda, y = estimate, color = term)) +
    geom_line(linewidth = 1) +
    labs(
        title = "lasso Regression Coefficient Paths",
        x = "log(lambda)", y = "Coefficient"
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "Variables"))
ggsave("lasso_coefficients_loglambda.png", plot = lasso_plot, dpi = 300, width = 8, height = 6, bg = "white")

# 选择最优 lambda
set.seed(2025)
cv.lasso <- cv.glmnet(x, y, alpha = 1)
png("lasso_cv_plot.png", width = 800, height = 600)
plot(cv.lasso, xlab = expression(log(lambda)))
dev.off()

# 得到参数
cv.lasso$lambda.min
cv.lasso$lambda.1se
coef(cv.lasso, s = "lambda.min")
coef(cv.lasso, s = "lambda.1se")
