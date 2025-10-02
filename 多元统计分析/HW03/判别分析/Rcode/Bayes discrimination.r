library(MASS)
library(rrcov)
library(ggplot2)
library(xtable)

## LDA函数
data(salmon)
lda.fit <- lda(
    Origin ~ Freshwater + Marine,
    prior = c(1, 1) / 2,
    data = salmon
)
print(lda.fit)

## 输出混淆矩阵
lda.predict <- predict(lda.fit, newdata = salmon)
table(salmon[, 4], lda.predict$class)

## 绘制线性判别函数决策边界
g.mean <- lda.fit$prior %*% lda.fit$means
const <- as.numeric(g.mean %*% lda.fit$scaling)
slope <- -lda.fit$scaling[1] / lda.fit$scaling[2]
intercept <- const / lda.fit$scaling[2]

png("lda_boundary_plot.png", width = 800, height = 800, res = 150)
plot(salmon[, 2:3],
    pch = rep(c(15, 20), each = 50),
    col = rep(c(2, 4), each = 50)
)
abline(intercept, slope)
legend("topleft", legend = c("Alaskan", "Canadian"), pch = c(15, 20), col = c(2, 4))
dev.off()


## QDA函数
qda.fit <- qda(
    Origin ~ Freshwater + Marine,
    prior = c(1, 1) / 2,
    data = salmon
)
print(qda.fit)

## 输出混淆矩阵
qda.predict <- predict(qda.fit, newdata = salmon)
table(salmon[, 4], qda.predict$class)

## 绘制二次判别函数决策边界
x <- seq(50, 200, 0.2)
y <- seq(300, 600, 0.2)
z <- as.matrix(expand.grid(x, y), 0)
colnames(z) <- c("Freshwater", "Marine")

z <- as.data.frame(z)
m <- length(x)
n <- length(y)

z.predict <- as.numeric(predict(object = qda.fit, newdata = z)$class)

png("qda_boundary_plot.png", width = 800, height = 800, res = 150)
plot(salmon[, 2:3],
    pch = rep(c(15, 20), each = 50),
    col = rep(c(2, 4), each = 50)
)

contour(x, y, matrix(z.predict, m, n),
    add = TRUE, drawlabels = FALSE, lty = 1
)

legend("topleft",
    legend = c("Alaskan", "Canadian"),
    pch = c(15, 20),
    col = c(2, 4)
)
dev.off()
