library(MASS)
library(knitr)

# 拟合LDA模型
lda1 <- lda(
    Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
    prior = rep(1 / 3, 3), data = iris
)
print(lda1)

# 预测
res1 <- predict(lda1)
print(res1)

# 混淆矩阵
newG <- res1$class
tab <- table(iris$Species, newG)

# 输出混淆矩阵
kable(tab, format = "latex", caption = "LDA混淆矩阵")

# 输出分类错误的样本
kable(cbind(iris, Predicted = newG)[iris$Species != newG, ],
    format = "latex", caption = "分类错误样本"
)
