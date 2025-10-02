library(msgps)
library(MASS)

# 加载 Boston 数据
data("Boston")
x <- as.matrix(Boston[, -14])
y <- Boston$medv

# 自适应Lasso回归
alasso_fit <- msgps(x, y, penalty = "alasso", gamma = 1, lambda = 0)

# 绘图
png("alasso_plot.png", width = 800, height = 400)
par(mfrow = c(1, 2))
plot(alasso_fit, criterion = "gcv", xvar = "t", main = "GCV")
plot(alasso_fit, criterion = "bic", xvar = "t", main = "BIC")
dev.off()

summary(alasso_fit)
