# 定义线性判别函数
W <- function(x, mu1, mu2, Sigma) {
    a <- solve(Sigma, mu1 - mu2)
    mu3 <- (mu1 + mu2) / 2
    sum(a * (x - mu3))
}

# 设置组均值和协方差矩阵
mu1 <- c(6, 5) # 雌虫平均值
mu2 <- c(8, 6) # 雄虫平均值
Sigma <- rbind(c(9, 2), c(2, 4)) # 协方差矩阵

# 输入待分类样本
x1 <- c(7.2, 5.6)
y1 <- W(x1, mu1, mu2, Sigma)

x2 <- c(6.3, 4.9)
y2 <- W(x2, mu1, mu2, Sigma)

# 输出结果
cat("当测量值为：", x1, "时，线性判别函数数值 =", y1, "\n")

cat("当测量值为：", x2, "时，线性判别函数数值 =", y2, "\n")
